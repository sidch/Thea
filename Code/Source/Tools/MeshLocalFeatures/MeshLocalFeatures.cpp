#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Local/AverageDistance.hpp"
#include "../../Algorithms/MeshFeatures/Local/Curvature.hpp"
#include "../../Algorithms/MeshFeatures/Local/LocalDistanceHistogram.hpp"
#include "../../Algorithms/MeshFeatures/Local/LocalPCA.hpp"
#include "../../Algorithms/MeshFeatures/Local/RandomWalks.hpp"
#include "../../Algorithms/MeshFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/MeshFeatures/Local/SpinImage.hpp"
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
bool is_oriented = false;  // all normals point outwards

int usage(int argc, char * argv[]);
double meshScale(MG & mg, MeshScaleType mesh_scale_type);
bool computeSDF(KDTree const & kdtree, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                TheaArray<double> & values);
bool computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                                long num_samples, double nbd_radius, TheaArray<double> & values);
bool computeAverageDistances(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, DistanceType dist_type,
                             double max_distance, TheaArray<double> & values);
bool computeLocalDistanceHistograms(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, long num_bins,
                                    DistanceType dist_type, double max_distance, double reduction_ratio,
                                    Matrix<double, MatrixLayout::ROW_MAJOR> & values);
bool computeLocalPCA(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, double nbd_radius, bool pca_full,
                     TheaArray<double> & values);
bool computeLocalPCARatios(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, double nbd_radius,
                           TheaArray<double> & values);
bool computeSpinImages(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, int num_radial_bins,
                       int num_height_bins, Matrix<double, MatrixLayout::ROW_MAJOR> & values);
bool computeRandomWalks(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, long num_steps, long num_walks,
                        Matrix<double, MatrixLayout::ROW_MAJOR> & values);

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
    else if (arg == "--abs")
    {
      abs_values = true;
    }
    else if (arg == "--binary")
    {
      binary = true;
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
    else if (arg == "--is-oriented")
    {
      is_oriented = true;
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
    else if (arg == "--shift01")
    {
      shift_to_01 = true;
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

    //=========================================================================================================================
    // Average Distance
    //=========================================================================================================================
    if (beginsWith(feat, "--avgd="))
    {
      char dist_str[260];
      DistanceType dist_type;
      long num_samples;
      double max_distance;

      int num_params = sscanf(feat.c_str(), "--avgd=%256[^,],%ld,%lf", dist_str, &num_samples, &max_distance);
      if (num_params < 1)
      {
        THEA_ERROR << "Couldn't parse average distance parameters";
        return -1;
      }

      if (!dist_type.fromString(dist_str))
      {
        THEA_ERROR << "Unknown distance metric: " << dist_str;
        return -1;
      }

      if (num_params < 3)
      {
        THEA_WARNING << "Distance limit not specified for average distance, using default of mesh scale";
        max_distance = -1;

        if (num_params < 2)
        {
          THEA_WARNING << "Number of samples not specified for average distance, using default value";
          num_samples = -1;
        }
      }

      TheaArray<double> values;
      if (!computeAverageDistances(mg, positions, num_samples, dist_type, max_distance, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of average distances doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    //=========================================================================================================================
    // Distance Histogram
    //=========================================================================================================================
    else if (beginsWith(feat, "--dh="))
    {
      char dist_str[260];
      DistanceType dist_type;
      long num_bins, num_samples;
      double max_distance;
      double reduction_ratio;

      int num_params = sscanf(feat.c_str(), "--dh=%256[^,],%ld,%ld,%lf,%lf",
                              dist_str, &num_bins, &num_samples, &max_distance, &reduction_ratio);
      if (num_params < 2)
      {
        THEA_ERROR << "Couldn't parse distance histogram parameters";
        return -1;
      }

      if (!dist_type.fromString(dist_str))
      {
        THEA_ERROR << "Unknown distance metric: " << dist_str;
        return -1;
      }

      if (num_params < 5)
      {
        THEA_WARNING << "Sample reduction ratio not specified for distance histogram, using default value";
        reduction_ratio = -1;

        if (num_params < 4)
        {
          THEA_WARNING << "Distance limit for not specified for distance histogram, using default of mesh scale";
          max_distance = -1;

          if (num_params < 3)
          {
            THEA_WARNING << "Number of samples not specified for distance histogram, using default value";
            num_samples = -1;
          }
        }
      }

      Matrix<double, MatrixLayout::ROW_MAJOR> values;  // each row is a histogram
      if (!computeLocalDistanceHistograms(mg, positions, num_samples, num_bins, dist_type, max_distance, reduction_ratio,
                                          values))
      {
        return -1;
      }

      alwaysAssertM(values.numRows() == (long)positions.size(), "Number of distance histograms doesn't match number of points");
      alwaysAssertM(positions.empty() || values.numColumns() == num_bins,
                    "Number of distance histogram bins doesn't match input parameter");

      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        double const * row_start = &values((long)j, 0);
        features[j].insert(features[j].end(), row_start, row_start + num_bins);
      }
    }
    //=========================================================================================================================
    // PCA
    //=========================================================================================================================
    else if (beginsWith(feat, "--pca"))
    {
      enum { PCA_DEFAULT, PCA_FULL, PCA_RATIO } pca_type = PCA_DEFAULT;
      array_size_t num_pca_features = 3;

      string opts;
      if (beginsWith(feat, "--pca-full"))
      {
        pca_type = PCA_FULL;
        num_pca_features = 12;
        opts = feat.substr(strlen("--pca-full"));
      }
      else if (beginsWith(feat, "--pca-ratio"))
      {
        pca_type = PCA_RATIO;
        num_pca_features = 2;
        opts = feat.substr(strlen("--pca-ratio"));
      }
      else
        opts = feat.substr(strlen("--pca"));

      double nbd_radius;
      long num_samples;
      int num_params = sscanf(opts.c_str(), "=%lf,%ld", &nbd_radius, &num_samples);
      if (num_params < 2)
      {
        THEA_WARNING << "Number of samples not specified for PCA, using default value";
        num_samples = -1;

        if (num_params < 1)
        {
          THEA_WARNING << "Neighborhood radius not specified for PCA, using default value";
          nbd_radius = -1;
        }
      }

      TheaArray<double> values;
      if (pca_type == PCA_RATIO)
      {
        if (!computeLocalPCARatios(mg, positions, num_samples, nbd_radius, values))
          return -1;
      }
      else
      {
        if (!computeLocalPCA(mg, positions, num_samples, nbd_radius, (pca_type == PCA_FULL), values))
          return -1;
      }

      alwaysAssertM(values.size() == num_pca_features * positions.size(),
                    "Number of PCA feature vectors doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        array_size_t base = num_pca_features * j;

        for (array_size_t k = 0; k < num_pca_features; ++k)
          features[j].push_back(values[base + k]);
      }
    }
    //=========================================================================================================================
    // Projected Curvature
    //=========================================================================================================================
    else if (beginsWith(feat, "--projcurv"))
    {
      double nbd_radius;
      long num_samples;
      int num_params = sscanf(feat.c_str(), "--projcurv=%lf,%ld", &nbd_radius, &num_samples);
      if (num_params < 2)
      {
        THEA_WARNING << "Number of samples not specified for projected curvature, using default value";
        num_samples = -1;

        if (num_params < 1)
        {
          THEA_WARNING << "Neighborhood radius not specified for projected curvature, using default value";
          nbd_radius = -1;
        }
      }

      TheaArray<double> values;
      if (!computeProjectedCurvatures(mg, positions, smooth_normals, num_samples, nbd_radius, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of projected curvatures doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    //=========================================================================================================================
    // Shape Diameter
    //=========================================================================================================================
    else if (feat == "--sdf")
    {
      TheaArray<double> values;
      if (!computeSDF(kdtree, positions, face_normals, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of SDF values doesn't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    //=========================================================================================================================
    // Spin Image
    //=========================================================================================================================
    else if (beginsWith(feat, "--spin="))
    {
      int num_radial_bins, num_height_bins;
      long num_samples;

      int num_params = sscanf(feat.c_str(), "--spin=%d,%d,%ld", &num_radial_bins, &num_height_bins, &num_samples);
      if (num_params < 2)
      {
        THEA_ERROR << "Couldn't parse spin image parameters";
        return -1;
      }

      if (num_params < 3)
      {
        THEA_WARNING << "Number of samples not specified for spin image, using default value";
        num_samples = -1;
      }

      Matrix<double, MatrixLayout::ROW_MAJOR> values;  // each row is the spin image at a query point
      if (!computeSpinImages(mg, positions, num_samples, num_radial_bins, num_height_bins, values))
      {
        return -1;
      }

      alwaysAssertM(values.numRows() == (long)positions.size(), "Number of spin images doesn't match number of points");
      alwaysAssertM(positions.empty() || values.numColumns() == num_radial_bins * num_height_bins,
                    "Number of spin image bins doesn't match input parameter");

      long num_features = num_radial_bins * num_height_bins;
      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        double const * row_start = &values((long)j, 0);
        features[j].insert(features[j].end(), row_start, row_start + num_features);
      }
    }
    //=========================================================================================================================
    // Random Walks
    //=========================================================================================================================
    else if (beginsWith(feat, "--rndwalk="))
    {
      long num_steps, num_walks, num_samples;

      int num_params = sscanf(feat.c_str(), "--rndwalk=%ld,%ld,%ld", &num_steps, &num_walks, &num_samples);
      if (num_params < 1)
      {
        THEA_ERROR << "Couldn't parse random walks parameters";
        return -1;
      }

      if (num_params < 3)
      {
        THEA_WARNING << "Number of samples not specified for random walks, using default value";
        num_samples = -1;

        if (num_params < 2)
        {
          THEA_WARNING << "Number of walks per point not specified for random walks, using default value";
          num_walks = -1;
        }
      }

      Matrix<double, MatrixLayout::ROW_MAJOR> values;  // each row is the feature vector of a query point
      if (!computeRandomWalks(mg, positions, num_samples, num_steps, num_walks, values))
      {
        return -1;
      }

      alwaysAssertM(values.numRows() == (long)positions.size(),
                    "Number of random walk descriptors doesn't match number of points");
      alwaysAssertM(positions.empty() || values.numColumns() == 3 * num_steps,
                    "Length of random walk descriptor != 3 * number of steps");

      long num_features = values.numColumns();
      for (array_size_t j = 0; j < positions.size(); ++j)
      {
        double const * row_start = &values((long)j, 0);
        features[j].insert(features[j].end(), row_start, row_start + num_features);
      }
    }
    else
    {
      if (!feat.empty())
        THEA_WARNING << "Ignoring unsupported option: " << feat;

      continue;
    }

    feat_names.push_back(feat.substr(2));
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

  THEA_CONSOLE << "Wrote " << features.size() << " feature vector(s) to " << out_path;

  return 0;
}

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <points> <outfile> [<feature0> <feature1> ...]";
  THEA_CONSOLE << "    <featureN> must be one of:";
  THEA_CONSOLE << "        --avgd=<metric>[,<#samples>[,<max-distance>]]";
  THEA_CONSOLE << "        --dh=<metric>,<#bins>[,<#samples>[,<max-distance>[,<reduction-ratio>]]]";
  THEA_CONSOLE << "        --pca[=<nbd-radius>[,<#samples>]]";
  THEA_CONSOLE << "              (eigenvalues in decreasing order)";
  THEA_CONSOLE << "        --pca-full[=<nbd-radius>[,<#samples>]]";
  THEA_CONSOLE << "              (eigenvalues + eigenvectors in decreasing order)";
  THEA_CONSOLE << "        --pca-ratio[=<nbd-radius>[,<#samples>]]";
  THEA_CONSOLE << "              (ratios of 2nd and 3rd eigenvalues to max eigenvalue)";
  THEA_CONSOLE << "        --projcurv[=<nbd-radius>[,<#samples>]]";
  THEA_CONSOLE << "        --rndwalk=<#steps>[,<#walks>[,<#samples>]]";
  THEA_CONSOLE << "        --sdf";
  THEA_CONSOLE << "        --spin=<#radial-bins>,<#height-bins>[,<#samples>]";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "    The following options may also be specified:";
  THEA_CONSOLE << "        --abs (uses the absolute value of every feature)";
  THEA_CONSOLE << "        --binary (outputs features in binary format)";
  THEA_CONSOLE << "        --featscale=<factor> (scales feature values by the factor)";
  THEA_CONSOLE << "        --is-oriented (assumes mesh normals consistently point outward)";
  THEA_CONSOLE << "        --meshscale={bsphere|bbox|avgdist} (used to set neighborhood scales)";
  THEA_CONSOLE << "        --normalize (rescale mesh so --meshscale == 1)";
  THEA_CONSOLE << "        --shift01 (maps features in [-1, 1] to [0, 1])";
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
  MeshFeatures::Local::ShapeDiameter<Mesh> sdf(&kdtree, (Real)mesh_scale);
  double scaling = (normalize_by_mesh_scale ? 1 : mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    double v0 = sdf.compute(positions[i], normals[i], true);
    if (v0 < 0)
      v0 = sdf.compute(positions[i], normals[i], false);

    if (is_oriented)
      values[i] = v0 * scaling;
    else
    {
      double v1 = sdf.compute(positions[i], -normals[i], true);
      if (v1 < 0)
        v1 = sdf.compute(positions[i], -normals[i], false);

      double vmin = (v1 < 0 || (v0 >= 0 && v0 < v1)) ? v0 : v1;
      values[i] = (vmin < 0 ? 0 : vmin) * scaling;
    }
  }

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                           long num_samples, double nbd_radius, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing projected curvatures";

  values.resize(positions.size());
  MeshFeatures::Local::Curvature<> projcurv(mg, num_samples, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
    values[i] = projcurv.computeProjectedCurvature(positions[i], normals[i], (Real)nbd_radius);

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeAverageDistances(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, DistanceType dist_type,
                        double max_distance, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing average " << dist_type.toString() << " distances";

  values.resize(positions.size());
  MeshFeatures::Local::AverageDistance<> avgd(mg, num_samples, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
    values[i] = avgd.compute(positions[i], dist_type, (Real)max_distance);

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeLocalDistanceHistograms(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, long num_bins,
                               DistanceType dist_type, double max_distance, double reduction_ratio,
                               Matrix<double, MatrixLayout::ROW_MAJOR> & values)
{
  THEA_CONSOLE << "Computing " << dist_type.toString() << " distance histograms";

  if (num_bins <= 0)
  {
    THEA_ERROR << "Number of histogram bins must be > 0";
    return false;
  }

  values.resize((long)positions.size(), num_bins);
  MeshFeatures::Local::LocalDistanceHistogram<> dh(mg, num_samples, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    Histogram histogram(num_bins, &values((long)i, 0));
    dh.compute(positions[i], histogram, dist_type, (Real)max_distance, (Real)reduction_ratio);
    histogram.normalize();
  }

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeLocalPCA(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, double nbd_radius, bool pca_full,
                TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing local PCA features";

  values.clear();
  MeshFeatures::Local::LocalPCA<> pca(mg, num_samples, (Real)mesh_scale);

  Vector3 evecs[3];
  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    Vector3 evals = pca.compute(positions[i], evecs, (Real)nbd_radius);
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

bool
computeLocalPCARatios(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, double nbd_radius,
                      TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing local PCA ratios";

  values.clear();
  MeshFeatures::Local::LocalPCA<> pca(mg, num_samples, (Real)mesh_scale);

  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    Vector3 evals = pca.compute(positions[i], NULL, (Real)nbd_radius);
    if (evals[0] > 0)
    {
      values.push_back(evals[1] / evals[0]);
      values.push_back(evals[2] / evals[0]);
    }
    else
    {
      values.push_back(0);
      values.push_back(0);
    }
  }

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeSpinImages(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, int num_radial_bins,
                  int num_height_bins, Matrix<double, MatrixLayout::ROW_MAJOR> & values)
{
  THEA_CONSOLE << "Computing spin images";

  long num_features = num_radial_bins * num_height_bins;
  values.resize((long)positions.size(), num_features);

  MeshFeatures::Local::SpinImage<> spin_image(mg, num_samples, (Real)mesh_scale);

  Matrix<double> f;
  for (array_size_t i = 0; i < positions.size(); ++i)
  {
    spin_image.compute(positions[i], num_radial_bins, num_height_bins, f);
    int j = 0;
    for (int r = 0; r < num_radial_bins; ++r)
      for (int h = 0; h < num_height_bins; ++h, ++j)
        values((long)i, j) = f(r, h);
  }

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeRandomWalks(MG const & mg, TheaArray<Vector3> const & positions, long num_samples, long num_steps, long num_walks,
                   Matrix<double, MatrixLayout::ROW_MAJOR> & values)
{
  THEA_CONSOLE << "Computing random walks";

  values.resize((long)positions.size(), 3 * (array_size_t)num_steps);
  MeshFeatures::Local::RandomWalks<> rw(mg, num_samples);

  for (array_size_t i = 0; i < positions.size(); ++i)
    rw.compute(positions[i], num_steps, &values((long)i, 0), num_walks);

  THEA_CONSOLE << "  -- done";

  return true;
}
