#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/ShapeDiameter.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Array.hpp"
#include "../../Vector3.hpp"
#include <boost/algorithm/string/trim.hpp>
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

bool computeSDF(MG const & mg, TheaArray<Vector3> const & pts, TheaArray<Real> & values);

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

  // Compute features
  TheaArray< TheaArray<Real> > features(pts.size());
  TheaArray<string> feat_names;

  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
      continue;

    string feat = string(argv[i]).substr(2);
    if (feat == "sdf")
    {
      TheaArray<Real> values;
      if (!computeSDF(mg, pts, values))
        return -1;

      alwaysAssertM(values.size() == pts.size(), "Number of feature values don't match number of points");

      for (array_size_t j = 0; j < pts.size(); ++j)
        features[j].push_back(values[j]);
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
    for (array_size_t j = 0; j < features[i].size(); ++j)
    {
      if (j > 0) out << ' ';
      out << features[i][j];
    }

    out << '\n';
  }

  THEA_CONSOLE << "Wrote " << features.size() << " feature vectors to " << out_path;

  return 0;
}

typedef MeshKDTree<Mesh> KDTree;

bool
computeSDF(MG const & mg, TheaArray<Vector3> const & pts, TheaArray<Real> & values)
{
  THEA_CONSOLE << "Computing SDF features";

  // Initialize kdtree
  KDTree kdtree;
  kdtree.add(const_cast<MG &>(mg));
  kdtree.init();

  THEA_CONSOLE << "  -- initialized mesh kdtree";

  // Find closest points and normals
  TheaArray<Vector3> positions;
  TheaArray<Vector3> normals;
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
    positions.push_back(cp);
    normals.push_back(cn);
  }

  THEA_CONSOLE << "  -- mapped query points to mesh";

  // Compute SDF values
  values.clear();
  MeshFeatures::ShapeDiameter<Mesh> sdf(&kdtree);
  sdf.compute(positions, normals, values);

  THEA_CONSOLE << "  -- computed SDF values";

  return true;
}
