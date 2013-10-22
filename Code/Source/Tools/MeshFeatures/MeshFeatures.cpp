#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Curvature.hpp"
#include "../../Algorithms/MeshFeatures/ShapeDiameter.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Array.hpp"
#include "../../Vector3.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <cstdio>
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

bool computeSDF(KDTree const & kdtree, TheaArray<Vector3> const & pts, TheaArray<Vector3> const & normals,
                TheaArray<Real> & values);
bool computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & pts, TheaArray<Vector3> const & normals,
                                TheaArray<Real> & values);

int
main(int argc, char * argv[])
{
  string mesh_path;
  string pts_path;
  string out_path;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
    {
      switch (curr_opt)
      {
        case 0: mesh_path = argv[i]; break;
        case 1: pts_path = argv[i]; break;
        case 2: out_path = argv[i]; break;
      }

      curr_opt++;
      if (curr_opt >= 3)
        break;
    }
  }

  if (curr_opt < 3)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <points> <outfile> [<feature0> <feature1> ...]";
    THEA_CONSOLE << "    <featureN> must be one of:";
    THEA_CONSOLE << "        --sdf";
    THEA_CONSOLE << "        --projcurv";
    THEA_CONSOLE << "        --shift01 (not a feature, maps features in [-1, 1] to [0, 1])";

    return 0;
  }

  MG mg;
  try
  {
    mg.load(mesh_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load mesh %s", mesh_path.c_str())

  // Load points
  TheaArray<Vector3> pts;
  {
    ifstream in(pts_path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not load points from file " << pts_path;
      return -1;
    }

    THEA_CONSOLE << "Loaded mesh from " << mesh_path;

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
  TheaArray< TheaArray<Real> > features(positions.size());
  TheaArray<string> feat_names;

  bool shift_to_01 = false;
  bool abs_values = false;
  bool scale = false;
  Real scale_factor = 1;

  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
      continue;

    string feat = string(argv[i]).substr(2);
    if (feat == "sdf")
    {
      TheaArray<Real> values;
      if (!computeSDF(kdtree, positions, face_normals, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of SDF values don't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    else if (feat == "projcurv")
    {
      TheaArray<Real> values;
      if (!computeProjectedCurvatures(mg, positions, smooth_normals, values))
        return -1;

      alwaysAssertM(values.size() == positions.size(), "Number of projected curvatures don't match number of points");

      for (array_size_t j = 0; j < positions.size(); ++j)
        features[j].push_back(values[j]);
    }
    else if (feat == "shift01")
    {
      shift_to_01 = true;
    }
    else if (feat == "abs")
    {
      abs_values = true;
    }
    else if (beginsWith(feat, "scale"))
    {
      if (sscanf(feat.c_str(), "scale=%f", &scale_factor) == 1)
        scale = true;
    }
    else
    {
      THEA_WARNING << "Ignoring unsupported feature type: " << feat;
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

  THEA_CONSOLE << "Computed " << feat_names.size() << " features: " << feat_str.str();

  // Write features to file
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
      Real f = features[i][j];

      if (scale) f *= scale_factor;
      if (shift_to_01) f = 0.5 * (1.0 + f);
      if (abs_values) f = fabs(f);

      out << ' ' << f;
    }

    out << '\n';
  }

  THEA_CONSOLE << "Wrote " << features.size() << " feature vectors to " << out_path;

  return 0;
}

bool
computeSDF(KDTree const & kdtree, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
           TheaArray<Real> & values)
{
  THEA_CONSOLE << "Computing SDF features";

  values.clear();
  MeshFeatures::ShapeDiameter<Mesh> sdf(&kdtree);
  sdf.compute(positions, normals, values);

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeProjectedCurvatures(MG const & mg, TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                           TheaArray<Real> & values)
{
  THEA_CONSOLE << "Computing projected curvatures";

  values.clear();
  MeshFeatures::Curvature<Mesh> projcurv(mg);
  projcurv.computeProjectedCurvatures(positions, normals, values);

  THEA_CONSOLE << "  -- done";

  return true;
}
