#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Curvature.hpp"
#include "../../Algorithms/MeshFeatures/DistanceHistogram.hpp"
#include "../../Algorithms/MeshFeatures/LocalPCA.hpp"
#include "../../Algorithms/MeshFeatures/ShapeDiameter.hpp"
#include "../../Algorithms/BestFitSphere3.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Array.hpp"
#include "../../IOStream.hpp"
#include "../../Matrix.hpp"
#include "../../Vector3.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Graphics;
using namespace Algorithms;

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

typedef MeshKDTree<Mesh> KDTree;

template <typename MeshTriangleT>
Vector3
smoothNormal(MeshTriangleT const & tri, Vector3 const & p)
{
  Vector3 b = tri.barycentricCoordinates(p);

  Vector3 n0 = tri.getVertices().getVertexNormal(0);
  Vector3 n1 = tri.getVertices().getVertexNormal(1);
  Vector3 n2 = tri.getVertices().getVertexNormal(2);

  return b[0] * n0 + b[1] * n1 + b[2] * n2;
}

struct MeshScaleType
{
  enum Value
  {
    BSPHERE,
    BBOX,
    AVG_DIST,
  };

  THEA_ENUM_CLASS_BODY(MeshScaleType)

  string toString() const
  {
    switch (*this)
    {
      case BSPHERE: return "bounding sphere";
      case BBOX: return "bounding box";
      case AVG_DIST: return "average distance from centroid";
      default: return "unknown";
    }
  }
};

MeshScaleType mesh_scale_type = MeshScaleType::BSPHERE;
bool normalize_by_mesh_scale = false;
double mesh_scale = 1;

int usage(int argc, char * argv[]);
double meshScale(MG & mg, MeshScaleType mesh_scale_type);
bool computeSDF(KDTree const & kdtree, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                TheaArray<double> & values);
bool computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                                TheaArray<double> & values);
bool computeDistanceHistograms(MG const & mg, TheaArray<Vector3> const & positions, long num_bins, double max_distance,
                               long num_samples, Matrix<double, MatrixLayout::ROW_MAJOR> & values);
bool computeLocalPCA(MG const & mg, TheaArray<Vector3> const & positions, bool pca_full, TheaArray<double> & values);

int
main(int argc, char * argv[])
{
  string mesh_path;
  string pts_path;
  string out_path;

  bool shift_to_01 = false;
  bool abs_values = false;
  bool feat_scale = false;
  double feat_scale_factor = 1;
  bool binary = false;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (!beginsWith(arg, "--"))
    {
      switch (curr_opt)
      {
        case 0: mesh_path = arg; break;
        case 1: pts_path = arg; break;
        case 2: out_path = arg; break;
        default:
        {
          THEA_ERROR << "Too many positional arguments";
          return -1;
        }
      }

      curr_opt++;
    }
    else if (arg == "--shift01")
    {
      shift_to_01 = true;
    }
    else if (arg == "--abs")
    {
      abs_values = true;
    }
    else if (beginsWith(arg, "--featscale="))
    {
      if (sscanf(arg.c_str(), "--featscale=%lf", &feat_scale_factor) == 1)
        feat_scale = true;
      else
      {
        THEA_ERROR << "Couldn't parse feat_scale factor";
        return -1;
      }
    }
    else if (beginsWith(arg, "--meshscale="))
    {
      string type = arg.substr(strlen("--meshscale="));
      if (type == "bsphere")
        mesh_scale_type = MeshScaleType::BSPHERE;
      else if (type == "bbox")
        mesh_scale_type = MeshScaleType::BBOX;
      else if (type == "avgdist")
        mesh_scale_type = MeshScaleType::AVG_DIST;
      else
      {
        THEA_ERROR << "Unsupported mesh scale type";
        return -1;
      }
    }
    else if (arg == "--normalize")
    {
      normalize_by_mesh_scale = true;
    }
    else if (arg == "--binary")
    {
      binary = true;
    }
    else
      continue;

    argv[i][0] = 0;  // zero out the argument so it won't be considered as a feature later
  }

  if (curr_opt < 3)
    return usage(argc, argv);

  MG mg;
  try
  {
    mg.load(mesh_path);
    mesh_scale = meshScale(mg, mesh_scale_type);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load mesh %s", mesh_path.c_str())

  THEA_CONSOLE << "Loaded mesh from " << mesh_path << " with scale " << mesh_scale << " (based on "
               << mesh_scale_type.toString() << ')';

  // Load points
  TheaArray<Vector3> pts;
  {
    ifstream in(pts_path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not load points from file " << pts_path;
      return -1;
    }

    string line;
    long line_num = 0;
    while (getline(in, line))
    {
      line_num++;

      boost::trim(line);
      if (line.empty())
        continue;

      istringstream line_in(line);
      Vector3 p;
      if (!(line_in >> p[0] >> p[1] >> p[2]))
      {
        THEA_ERROR << "Could not read point on line " << line_num << " of file " << pts_path;
        return -1;
      }

      pts.push_back(p);
    }
  }

  THEA_CONSOLE << "Loaded " << pts.size() << " point(s) from " << pts_path;

  // Snap to surface points and normals
  KDTree kdtree;
  kdtree.add(mg);
  kdtree.init();

  THEA_CONSOLE << "Created mesh kd-tree";

  TheaArray<Vector3> positions(pts.size());
  TheaArray<Vector3> face_normals(pts.size());
  TheaArray<Vector3> smooth_normals(pts.size());

  for (array_size_t i = 0; i < pts.size(); ++i)
  {
    Vector3 cp(0);
    long elem = kdtree.closestElement<MetricL2>(pts[i], -1, NULL, &cp);
    if (elem < 0)
    {
      THEA_ERROR << "Could not find nearest neighbor of query point " << pts[i] << " on mesh";
      return false;
    }

    Vector3 cn = kdtree.getElements()[(array_size_t)elem].getNormal();
    Vector3 csn = smoothNormal(kdtree.getElements()[(array_size_t)elem], cp);

    positions[i] = cp;
    face_normals[i] = cn;
    smooth_normals[i] = csn;
  }

  THEA_CONSOLE << "Snapped query points to mesh";

  // Compute features
  TheaArray< TheaArray<double> > features(positions.size());
  TheaArray<string> feat_names;

  for (int i = 1; i < argc; ++i)
  {
    string feat = string(argv[i]);
    if (feat == "--sdf")
    {
      TheaArray<double> values;
      if (!computeSDF(kdtree, positions, face_normals, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of SDF values doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    else if (feat == "--projcurv")
    {
      TheaArray<double> values;
      if (!computeProjectedCurvatures(mg, positions, smooth_normals, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of projected curvatures doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    else if (beginsWith(feat, "--dh="))
    {
      long num_bins, num_samples;
      double max_distance;

      long num_params = sscanf(feat.c_str(), "--dh=%ld,%lf,%ld", &num_bins, &max_distance, &num_samples);
      if (num_params < 1)
      {
        THEA_ERROR << "Couldn't parse distance histogram parameters";
        return -1;
      }
      else if (num_params < 3)
      {
        THEA_WARNING << "Approximate number of samples not specified for distance histogram, using default value";
        num_samples = -1;

        if (num_params < 2)
        {
          THEA_WARNING << "Distance limit for not specified for distance histogram, using default of mesh scale";
          max_distance = -1;
        }
      }

      if (num_bins <= 0)
      {
        THEA_ERROR << "Number of histogram bins must be > 0";
        return -1;
      }

      Matrix<double, MatrixLayout::ROW_MAJOR> values;  // each row is a histogram
      if (!computeDistanceHistograms(mg, positions, num_bins, max_distance, num_samples, values))
        return -1;

      alwaysAssertM(values.numRows() == (long)positions.size(), "Number of distance histograms doesn't match number of points");
      alwaysAssertM(positions.empty() || values.numColumns() == num_bins,
                    "Number of distance histogram bins doesn't match input parameter");

      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        double const * row_start = &values((long)j, 0);
        features[j].insert(features[j].end(), row_start, row_start + num_bins);
      }
    }
    else if (beginsWith(feat, "--pca"))
    {
      bool pca_full = (feat == "--pca=full");

      TheaArray<double> values;
      if (!computeLocalPCA(mg, positions, pca_full, values))
        return -1;

      array_size_t step = (pca_full ? 12 : 3);
      alwaysAssertM(values.size() == step * positions.size(), "Number of PCA feature vectors doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        array_size_t base = step * j;

        for (array_size_t k = 0; k < step; ++k)
          features[j].push_back(values[base + k]);
      }
    }
    else
    {
      if (!feat.empty())
        THEA_WARNING << "Ignoring unsupported option: " << feat;

      continue;
    }

    feat_names.push_back(feat);
  }

  ostringstream feat_str;
  for (array_size_t i = 0; i < feat_names.size(); ++i)
  {
    if (i > 0) feat_str << ", ";
    feat_str << feat_names[i];
  }

  THEA_CONSOLE << "Computed " << feat_names.size() << " feature(s): " << feat_str.str();

  // Write features to file
  if (binary)
  {
    BinaryOutputStream out(out_path, Endianness::LITTLE);
    if (!out.ok())
    {
      THEA_ERROR << "Could not open output file " << out_path << " for writing features";
      return -1;
    }

    out.writeInt64((int64)features.size());
    out.writeInt64((int64)(features.empty() ? 0 : features[0].size()));

    for (array_size_t i = 0; i < features.size(); ++i)
    {
      out.writeFloat32(pts[i][0]);
      out.writeFloat32(pts[i][1]);
      out.writeFloat32(pts[i][2]);

      for (array_size_t j = 0; j < features[i].size(); ++j)
      {
        double f = features[i][j];

        if (feat_scale) f *= feat_scale_factor;
        if (shift_to_01) f = 0.5 * (1.0 + f);
        if (abs_values) f = fabs(f);

        out.writeFloat32((float32)f);
      }
    }
  }
  else
  {
    ofstream out(out_path.c_str());
    if (!out)
    {
      THEA_ERROR << "Could not open output file " << out_path << " for writing features";
      return -1;
    }

    for (array_size_t i = 0; i < features.size(); ++i)
    {
      out << pts[i][0] << ' ' << pts[i][1] << ' ' << pts[i][2];

      for (array_size_t j = 0; j < features[i].size(); ++j)
      {
        double f = features[i][j];

        if (feat_scale) f *= feat_scale_factor;
        if (shift_to_01) f = 0.5 * (1.0 + f);
        if (abs_values) f = fabs(f);

        out << ' ' << f;
      }

      out << '\n';
    }
  }

  THEA_CONSOLE << "Wrote " << features.size() << " feature vectors to " << out_path;

  return 0;
}

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <points> <outfile> [<feature0> <feature1> ...]";
  THEA_CONSOLE << "    <featureN> must be one of:";
  THEA_CONSOLE << "        --sdf";
  THEA_CONSOLE << "        --projcurv";
  THEA_CONSOLE << "        --dh=<num-bins>[,<max_distance>[,<num-samples>]]";
  THEA_CONSOLE << "        --pca[=full] (eigenvalues [+ eigenvectors] in decreasing order)";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "    The following options may also be specified:";
  THEA_CONSOLE << "        --meshscale={bsphere|bbox|avgdist} (used to set neighborhood scales)";
  THEA_CONSOLE << "        --normalize (rescale mesh so --meshscale == 1)";
  THEA_CONSOLE << "        --shift01 (maps features in [-1, 1] to [0, 1])";
  THEA_CONSOLE << "        --featscale=<factor> (scales feature values by the factor)";
  THEA_CONSOLE << "        --binary (outputs features in binary format)";
  THEA_CONSOLE << "";

  return -1;
}

double
meshScale(MG & mg, MeshScaleType mesh_scale_type)
{
  switch (mesh_scale_type)
  {
    case MeshScaleType::BBOX:
    {
      mg.updateBounds();
      return mg.getBounds().getExtent().length();
    }

    case MeshScaleType::AVG_DIST:
    {
      MeshSampler<Mesh> sampler(mg);
      TheaArray<Vector3> samples;
      sampler.sampleEvenlyByArea(50000, samples);

      if (samples.size() <= 1)
        return 0.0;

      Vector3 centroid = CentroidN<Vector3, 3>::compute(samples.begin(), samples.end());

      double sum_dists = 0;
      for (array_size_t i = 0; i < samples.size(); ++i)
        sum_dists += (samples[i] - centroid).length();

      return (3.0 * sum_dists) / samples.size();  // 3 instead of 2 since the avg is in general smaller than the max
    }

    default:
    {
      BestFitSphere3 bsphere;
      PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mg);
      return bsphere.getDiameter();
    }
  }
}

bool
computeSDF(KDTree const & kdtree, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
           TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing SDF features";

  values.resize(positions.size());
  MeshFeatures::ShapeDiameter<Mesh> sdf(&kdtree, (Real)mesh_scale);
  double scaling = (normalize_by_mesh_scale ? 1 : mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    double v0 = sdf.compute(positions[i],  normals[i], true);
    if (v0 < 0)
      v0 = sdf.compute(positions[i],  normals[i], false);

    double v1 = sdf.compute(positions[i], -normals[i], true);
    if (v1 < 0)
      v1 = sdf.compute(positions[i], -normals[i], false);

    double vmin = (v1 < 0 || (v0 >= 0 && v0 < v1)) ? v0 : v1;
    values[i] = (vmin < 0 ? 0 : vmin) * scaling;
  }

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                           TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing projected curvatures";

  values.resize(positions.size());
  MeshFeatures::Curvature<Mesh> projcurv(mg, -1, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
    values[i] = projcurv.computeProjectedCurvature(positions[i], normals[i]);

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeDistanceHistograms(MG const & mg, TheaArray<Vector3> const & positions, long num_bins, double max_distance,
                          long num_samples, Matrix<double, MatrixLayout::ROW_MAJOR> & values)
{
  THEA_CONSOLE << "Computing distance histograms";

  values.resize((long)positions.size(), num_bins);
  MeshFeatures::DistanceHistogram<Mesh> dh(mg, num_samples, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
    dh.compute(positions[i], num_bins, &values((long)i, 0), max_distance);

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeLocalPCA(MG const & mg, TheaArray<Vector3> const & positions, bool pca_full, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing local PCA features";

  values.clear();
  MeshFeatures::LocalPCA<Mesh> pca(mg, -1, (Real)mesh_scale);

  Vector3 evecs[3];
  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    Vector3 evals = pca.compute(positions[i], evecs);
    if (normalize_by_mesh_scale)
      evals /= mesh_scale;

    values.push_back(evals[0]);
    values.push_back(evals[1]);
    values.push_back(evals[2]);

    if (pca_full)
    {
      for (array_size_t j = 0; j < 3; ++j)
        for (array_size_t k = 0; k < 3; ++k)
          values.push_back(evecs[j][k]);
    }
  }

  THEA_CONSOLE << "  -- done";

  return true;
}
