#include "../../Common.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/ICP3.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../AffineTransformN.hpp"
#include "../../MatrixMN.hpp"
#include "../../VectorN.hpp"
#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

long num_mesh_samples = 5000;
bool normalize_scales = false;
bool prealign_centroids = false;
bool rotate_axis_aligned = false;

typedef VectorN<3, double> DVector3;
typedef MatrixMN<3, 3, double> DMatrix3;
typedef AffineTransformN<3, double> DAffineTransform3;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [<options>] <from> <to> [<output>]";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -n <#mesh-samples>  :  Approximate #points to sample from mesh";
  THEA_CONSOLE << "  -s                  :  Normalize shape scales";
  THEA_CONSOLE << "  -c                  :  Align shape centroids before starting ICP";
  THEA_CONSOLE << "  -x                  :  Test over various axis-aligned rotations";
  // THEA_CONSOLE << "  -r                  :  Test over various rotations";
  return -1;
}

bool
sampleMesh(string const & mesh_path, TheaArray<DVector3> & samples)
{
  typedef DisplayMesh Mesh;
  typedef MeshGroup<Mesh> MG;

  try
  {
    MG mg("MeshGroup");
    mg.load(mesh_path);

    MeshSampler<Mesh> sampler(mg);
    TheaArray<Vector3> pts;
    sampler.sampleEvenlyByArea(num_mesh_samples, pts);

    samples.resize(pts.size());
    for (array_size_t i = 0; i < samples.size(); ++i)
      samples[i] = DVector3(pts[i]);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "An error occurred")

  THEA_CONSOLE << "Sampled " << samples.size() << " points from mesh " << mesh_path;

  return true;
}

bool
loadPoints(string const & points_path, TheaArray<DVector3> & points)
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

    DVector3 p;
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
    DVector3 p;
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

DAffineTransform3
normalization(TheaArray<DVector3> const & from_pts, TheaArray<DVector3> const & to_pts, double & from_scale, double & to_scale)
{
  from_scale = to_scale = 1.0;

  if (from_pts.empty() || to_pts.empty())
    return DAffineTransform3::identity();

  DVector3 from_center  =  CentroidN<DVector3, 3, double>::compute(from_pts.begin(), from_pts.end());
  DVector3 to_center    =  CentroidN<DVector3, 3, double>::compute(to_pts.begin(), to_pts.end());

  double rescaling = 1.0;
  if (normalize_scales)
  {
    // Measure average distance of from_pts to from_center
    double from_avg_dist = 0;
    for (array_size_t i = 0; i < from_pts.size(); ++i)
      from_avg_dist += (from_pts[i] - from_center).length();

    from_avg_dist /= from_pts.size();

    // Measure average distance of to_pts to to_center
    double to_avg_dist = 0;
    for (array_size_t i = 0; i < to_pts.size(); ++i)
      to_avg_dist += (to_pts[i] - to_center).length();

    to_avg_dist /= to_pts.size();

    if (from_avg_dist > 0 && to_avg_dist > 0)
      rescaling = to_avg_dist / from_avg_dist;

    from_scale  =  from_avg_dist;
    to_scale    =  to_avg_dist;

    THEA_CONSOLE << "Rescaling = " << rescaling;
  }

  if (normalize_scales)
  {
    DAffineTransform3 tr = DAffineTransform3::scaling(rescaling)
                         * DAffineTransform3::translation(-from_center);

    if (prealign_centroids)
      return DAffineTransform3::translation(to_center) * tr;
    else
      return DAffineTransform3::translation(from_center) * tr;  // scale in place
  }
  else
  {
    if (prealign_centroids)
      return DAffineTransform3::translation(to_center - from_center);
    else
      return DAffineTransform3::identity();
  }
}

double
normalizeError(double err, double scale, long num_points)
{
  if (normalize_scales && scale > 0)
    return err / (4 * scale * scale * num_points);  // normalized scale makes average distance from centroid == 0.5
  else
    return err / num_points;
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
        normalize_scales = true;
        THEA_CONSOLE << "Scales will be matched";
      }
      else if (arg == "-c")
      {
        prealign_centroids = true;
        THEA_CONSOLE << "Centroids will be pre-aligned";
      }
      else if (arg == "-x")
      {
        rotate_axis_aligned = true;
        THEA_CONSOLE << "Alignment will search over axis-aligned rotations";
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

  TheaArray<DVector3> from_pts;
  if (!loadPoints(from_path, from_pts))
    return -1;

  TheaArray<DVector3> to_pts;
  if (!loadPoints(to_path, to_pts))
    return -1;

  DAffineTransform3 tr = DAffineTransform3::identity();
  double from_scale = 1.0, to_scale = 1.0;
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
    DAffineTransform3 init_tr = normalization(from_pts, to_pts, from_scale, to_scale);
    for (array_size_t i = 0; i < from_pts.size(); ++i)
      from_pts[i] = init_tr * from_pts[i];

    KDTreeN<DVector3, 3, double> to_kdtree(to_pts.begin(), to_pts.end());

    if (rotate_axis_aligned)
    {
      ICP3<double> icp(-1, -1, false);
      DVector3 from_center = CentroidN<DVector3, 3, double>::compute(from_pts.begin(), from_pts.end());
      TheaArray<DVector3> rot_from_pts(from_pts.size());
      double rot_error;
      bool first = true;

      for (int u = 0; u < 2; ++u)
        for (int su = -1; su <= 1; su += 2)
          for (int v = 1; v < 2; ++v)
            for (int sv = -1; sv <= 1; sv += 2)
            {
              DVector3 du(0, 0, 0); du[u] = su;
              DVector3 dv(0, 0, 0); dv[(u + v) % 3] = sv;
              DVector3 dw = du.cross(dv);
              DMatrix3 rot(du[0], dv[0], dw[0],
                           du[1], dv[1], dw[1],
                           du[2], dv[2], dw[2]);

              DAffineTransform3 rot_tr = DAffineTransform3::translation(from_center)
                                       * DAffineTransform3(rot)
                                       * DAffineTransform3::translation(-from_center);

              for (array_size_t i = 0; i < from_pts.size(); ++i)
                rot_from_pts[i] = rot_tr * from_pts[i];

              DAffineTransform3 icp_tr = icp.align((long)rot_from_pts.size(), &rot_from_pts[0], to_kdtree, &rot_error);
              rot_error = normalizeError(rot_error, to_scale, (long)from_pts.size());
              if (first || rot_error < error)
              {
                error = rot_error;
                tr = icp_tr * rot_tr;
                first = false;

                THEA_CONSOLE << "-- rotation " << rot.toString() << " reduced error to " << error;
              }
              else
                THEA_CONSOLE << "-- rotation " << rot.toString() << ", no reduction in error";
            }

      alwaysAssertM(!first, "No alignment found after trying axis-aligned rotations");

      tr = tr * init_tr;
    }
    else
    {
      ICP3<double> icp(-1, -1, true);
      tr = icp.align((long)from_pts.size(), &from_pts[0], to_kdtree, &error);
      tr = tr * init_tr;

      error = normalizeError(error, to_scale, (long)from_pts.size());
    }
  }

  THEA_CONSOLE << "Alignment error = " << std::setprecision(10) << error;
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
