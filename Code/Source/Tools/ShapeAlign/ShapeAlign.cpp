#include "../../Common.hpp"
#include "../../Algorithms/ICP3.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " <from> <to> [<output>]";
  return -1;
}

bool
sampleMesh(string const & mesh_path, TheaArray<Vector3> & samples)
{
  typedef DisplayMesh Mesh;
  typedef MeshGroup<Mesh> MG;

  static int const NUM_SAMPLES = 50000;

  try
  {
    MG mg("MeshGroup");
    mg.load(mesh_path);

    MeshSampler<Mesh> sampler(mg);
    samples.clear();
    sampler.sampleEvenlyByArea(NUM_SAMPLES, samples);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "An error occurred")

  THEA_CONSOLE << "Sampled " << samples.size() << " points from mesh " << mesh_path;

  return true;
}

bool
loadPoints(string const & points_path, TheaArray<Vector3> & points)
{
  string path_lc = toLower(points_path);
  if (endsWith(path_lc, ".obj") || endsWith(path_lc, ".3ds") || endsWith(path_lc, ".ply"))  // this is a mesh file!
    return sampleMesh(points_path, points);

  ifstream in(points_path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open points file " << points_path;
    return false;
  }

  points.clear();

  if (endsWith(path_lc, ".off"))
  {
    string magic;
    in >> magic;
    if (magic != "OFF")
    {
      THEA_ERROR << "OFF file '" << points_path << "' does not start with 'OFF'";
      return false;
    }

    long nv, nf, ne;
    if (!(in >> nv >> nf >> ne))
    {
      THEA_ERROR << "Could not read vertex/face/edge counts from " << points_path;
      return false;
    }

    if (nv < 0)
    {
      THEA_ERROR << "Invalid vertex count " << nv << " in " << points_path;
      return false;
    }

    if (nf > 0)  // this is a mesh file! close, reopen as mesh, and sample
    {
      return sampleMesh(points_path, points);
    }

    points.resize((array_size_t)nv);

    for (long i = 0; i < nv; ++i)
      if (!(in >> points[i][0] >> points[i][1] >> points[i][2]))
      {
        THEA_ERROR << "Could not read vertex " << i << " from " << points_path;
        return false;
      }
  }
  else if (endsWith(path_lc, ".arff"))
  {
    bool seen_data = false;
    string line;
    while (getline(in, line))
    {
      if (beginsWith(toUpper(trimWhitespace(line)), "@DATA"))
      {
        seen_data = true;
        break;
      }
    }

    if (!seen_data)
    {
      THEA_ERROR << "No points found in " << points_path;
      return false;
    }

    Vector3 p;
    string field;
    while (getline(in, line))
    {
      line = trimWhitespace(line);
      if (line.empty())
        continue;

      istringstream line_in(line);
      for (int i = 0; i < 3; ++i)
      {
        if (!getline(line_in, field, ','))
        {
          THEA_ERROR << "Could not read coordinate of " << i << " of point " << points.size() << " from " << points_path;
          return false;
        }

        istringstream field_in(field);
        if (!(field_in >> p[i]))
        {
          THEA_ERROR << "Could not parse coordinate of " << i << " of point " << points.size() << " from " << points_path;
          return false;
        }
      }

      points.push_back(p);
    }
  }
  else
  {
    Vector3 p;
    string line;
    while (getline(in, line))
    {
      line = trimWhitespace(line);
      if (line.empty())
        continue;

      istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2]))
      {
        THEA_ERROR << "Could not read point " << points.size() << " from " << points_path;
        return false;
      }

      points.push_back(p);
    }
  }

  THEA_CONSOLE << "Read " << points.size() << " from " << points_path;

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
    return usage(argc, argv);

  string from_path = argv[1];
  string to_path = argv[2];
  string out_path = (argc > 3 ? argv[3] : "");

  TheaArray<Vector3> from_pts;
  if (!loadPoints(from_path, from_pts))
    return -1;

  TheaArray<Vector3> to_pts;
  if (!loadPoints(to_path, to_pts))
    return -1;

  AffineTransformN<3, double> tr = AffineTransformN<3, double>::identity();
  double error = 0;
  if (from_pts.empty())
  {
    THEA_WARNING << "Source shape is empty: returning identity transform";
  }
  else if (to_pts.empty())
  {
    THEA_WARNING << "Target shape is empty: returning identity transform";
  }
  else
  {
    ICP3<double> icp(-1, -1, true);
    tr = icp.align((long)from_pts.size(), &from_pts[0], (long)to_pts.size(), &to_pts[0], &error);
  }

  THEA_CONSOLE << "Alignment error = " << error;
  THEA_CONSOLE << "Alignment = " << tr.toString();

  return 0;
}
