#include "../../Common.hpp"
#include "../../Algorithms/BestFitSphere3.hpp"
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

long num_mesh_samples = 5000;
bool match_scale = false;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [<options>] <from> <to> [<output>]";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -n <num-mesh-samples>  :  Approximate #points to sample from mesh";
  THEA_CONSOLE << "  -s                     :  Match shape scales";
  return -1;
}

bool
sampleMesh(string const & mesh_path, TheaArray<Vector3> & samples)
{
  typedef DisplayMesh Mesh;
  typedef MeshGroup<Mesh> MG;

  try
  {
    MG mg("MeshGroup");
    mg.load(mesh_path);

    MeshSampler<Mesh> sampler(mg);
    samples.clear();
    sampler.sampleEvenlyByArea(num_mesh_samples, samples);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "An error occurred")

  THEA_CONSOLE << "Sampled " << samples.size() << " points from mesh " << mesh_path;

  return true;
}

bool
loadPoints(string const & points_path, TheaArray<Vector3> & points)
{
  string path_lc = toLower(points_path);
  if (endsWith(path_lc, ".obj")
   || endsWith(path_lc, ".3ds")
   || endsWith(path_lc, ".ply")
   || endsWith(path_lc, ".off.bin"))  // this is a mesh file!
  {
    return sampleMesh(points_path, points);
  }

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

  THEA_CONSOLE << "Read " << points.size() << " points from " << points_path;

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
    return usage(argc, argv);

  int pos_arg = 0;
  string from_path;
  string to_path;
  string out_path;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (arg.empty())
      continue;

    if (arg[0] == '-')
    {
      if (arg == "-n")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        istringstream iss(argv[++i]);
        if (!(iss >> num_mesh_samples))
        {
          THEA_ERROR << "Could not parse number of mesh samples";
          return -1;
        }

        if (num_mesh_samples < 1)
        {
          THEA_ERROR << "Invalid number of samples: " << num_mesh_samples;
          return -1;
        }

        THEA_CONSOLE << "Number of mesh samples: " << num_mesh_samples;
      }
      else if (arg == "-s")
      {
        match_scale = true;
        THEA_CONSOLE << "Scales will be matched";
      }
      else
        return usage(argc, argv);
    }
    else
    {
      switch (pos_arg)
      {
        case 0: from_path = arg; break;
        case 1: to_path = arg; break;
        case 2: out_path = arg; break;
        default: return usage(argc, argv);
      }

      pos_arg++;
    }
  }

  if (pos_arg < 2)
    return usage(argc, argv);

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
    double rescaling = 1.0;
    if (match_scale)
    {
      if (from_pts.size() > 1 && to_pts.size() > 1)
      {
        BestFitSphere3 from_bsphere;
        from_bsphere.addPoints(from_pts.begin(), from_pts.end());

        BestFitSphere3 to_bsphere;
        to_bsphere.addPoints(to_pts.begin(), to_pts.end());

        if (from_bsphere.getRadius() > 0)
          rescaling = (double)to_bsphere.getRadius() / (double)from_bsphere.getRadius();
      }

      THEA_CONSOLE << "Rescaling = " << rescaling;

      for (array_size_t i = 0; i < from_pts.size(); ++i)
        from_pts[i] *= rescaling;
    }

    ICP3<double> icp(-1, -1, true);
    tr = icp.align((long)from_pts.size(), &from_pts[0], (long)to_pts.size(), &to_pts[0], &error);

    if (match_scale)
      tr = tr * AffineTransformN<3, double>::scaling(rescaling);
  }

  THEA_CONSOLE << "Alignment error = " << error;
  THEA_CONSOLE << "Alignment = " << tr.toString();

  if (!out_path.empty())
  {
    ofstream out(out_path.c_str());
    if (!out)
    {
      THEA_ERROR << "Couldn't open output file " << out_path;
      return -1;
    }

    out << tr.getLinear()(0, 0) << ' ' << tr.getLinear()(0, 1) << ' ' << tr.getLinear()(0, 2) << tr.getTranslation()[0] << '\n'
        << tr.getLinear()(1, 0) << ' ' << tr.getLinear()(1, 1) << ' ' << tr.getLinear()(1, 2) << tr.getTranslation()[1] << '\n'
        << tr.getLinear()(2, 0) << ' ' << tr.getLinear()(2, 1) << ' ' << tr.getLinear()(2, 2) << tr.getTranslation()[2] << '\n';
  }

  return 0;
}
