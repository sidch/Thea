#include "../../Common.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/Icp3.hpp"
#include "../../Algorithms/KdTreeN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../AffineTransformN.hpp"
#include "../../FilePath.hpp"
#include "../../MatVec.hpp"
#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef AffineTransformN<3, double> AffineTransform3d;

intx num_mesh_samples = 5000;
bool normalize_scales = false;
bool do_initial_placement = false;
Vector3d initial_from_position(0, 0, 0);
bool prealign_centroids = false;
bool rotate_axis_aligned = false;
bool rotate_arbitrary = false;
intx num_rotations = -1;
bool has_up_vector = false;
Vector3d up_vector;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [<options>] <from> <to> [<output>]";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -n <#mesh-samples>  :  Approximate #points to sample from mesh";
  THEA_CONSOLE << "  -s                  :  Normalize shape scales";
  THEA_CONSOLE << "  -t x y z            :  Initial location for centroid of <from>";
  THEA_CONSOLE << "  -c                  :  Align shape centroids before starting ICP";
  THEA_CONSOLE << "  -x                  :  Test over various axis-aligned rotations";
  THEA_CONSOLE << "  -r <#rotations>     :  Test over rotations around an up vector (requires -u).";
  THEA_CONSOLE << "                         <#rotations> <= 0 picks a default number.";
  THEA_CONSOLE << "  -u x y z            :  Up vector";
  return -1;
}

bool
sampleMesh(string const & mesh_path, Array<Vector3d> & samples)
{
  typedef DisplayMesh Mesh;
  typedef MeshGroup<Mesh> MG;

  try
  {
    MG mg("MeshGroup");
    mg.load(mesh_path);

    MeshSampler<Mesh> sampler(mg);
    Array<Vector3> pts;
    sampler.sampleEvenlyByArea(num_mesh_samples, pts);

    samples.resize(pts.size());
    for (size_t i = 0; i < samples.size(); ++i)
      samples[i] = pts[i].cast<double>();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "An error occurred")

  THEA_CONSOLE << "Sampled " << samples.size() << " points from mesh " << mesh_path;

  return true;
}

bool
loadPoints(string const & points_path, Array<Vector3d> & points)
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

    intx nv, nf, ne;
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

    points.resize((size_t)nv);

    for (intx i = 0; i < nv; ++i)
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

    Vector3d p;
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
    Vector3d p;
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

AffineTransform3d
normalization(Array<Vector3d> const & from_pts, Array<Vector3d> const & to_pts, double & from_scale, double & to_scale)
{
  from_scale = to_scale = 1.0;

  if (from_pts.empty() || to_pts.empty())
    return AffineTransform3d::identity();

  Vector3d from_center  =  CentroidN<Vector3d, 3, double>::compute(from_pts.begin(), from_pts.end());
  Vector3d to_center    =  CentroidN<Vector3d, 3, double>::compute(to_pts.begin(), to_pts.end());

  double rescaling = 1.0;
  if (normalize_scales)
  {
    // Measure average distance of from_pts to from_center
    double from_avg_dist = 0;
    for (size_t i = 0; i < from_pts.size(); ++i)
      from_avg_dist += (from_pts[i] - from_center).norm();

    from_avg_dist /= from_pts.size();

    // Measure average distance of to_pts to to_center
    double to_avg_dist = 0;
    for (size_t i = 0; i < to_pts.size(); ++i)
      to_avg_dist += (to_pts[i] - to_center).norm();

    to_avg_dist /= to_pts.size();

    if (from_avg_dist > 0 && to_avg_dist > 0)
      rescaling = to_avg_dist / from_avg_dist;

    from_scale  =  from_avg_dist;
    to_scale    =  to_avg_dist;

    THEA_CONSOLE << "Rescaling = " << rescaling;
  }

  if (normalize_scales)
  {
    AffineTransform3d tr = AffineTransform3d::scaling(rescaling)
                         * AffineTransform3d::translation(-from_center);

    if (prealign_centroids)
      return AffineTransform3d::translation(to_center) * tr;
    else
      return AffineTransform3d::translation(from_center) * tr;  // scale in place
  }
  else
  {
    if (prealign_centroids)
      return AffineTransform3d::translation(to_center - from_center);
    else
      return AffineTransform3d::identity();
  }
}

double
normalizeError(double err, double scale, intx num_points)
{
  if (normalize_scales && scale > 0)
    return err / (4 * scale * scale * num_points);  // normalized scale makes average distance from centroid == 0.5
  else
    return err / num_points;
}

bool
loadShapePaths(string const & path, Array<string> & shape_paths)
{
  shape_paths.clear();

  if (endsWith(toLower(path), ".txt"))
  {
    ifstream in(path.c_str());
    if (!in)
    {
      THEA_ERROR << "Couldn't open list of shape paths: " << path;
      return false;
    }

    string dir = FilePath::parent(path);
    string line;
    while (getline(in, line))
    {
      line = trimWhitespace(line);
      if (!line.empty())
        shape_paths.push_back(FilePath::concat(dir, line));
    }
  }
  else
    shape_paths.push_back(path);

  return true;
}

int
alignShapes(string const & from_path, string const & to_path, std::ostream * out)
{
  Array<Vector3d> from_pts;
  if (!loadPoints(from_path, from_pts))
    return -1;

  Array<Vector3d> to_pts;
  if (!loadPoints(to_path, to_pts))
    return -1;

  AffineTransform3d tr = AffineTransform3d::identity();
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
    AffineTransform3d init_tr = normalization(from_pts, to_pts, from_scale, to_scale);
    for (size_t i = 0; i < from_pts.size(); ++i)
      from_pts[i] = init_tr * from_pts[i];

    Vector3d from_center = CentroidN<Vector3d, 3, double>::compute(from_pts.begin(), from_pts.end());
    if (do_initial_placement)
    {
      Vector3d offset = initial_from_position - from_center;
      init_tr = AffineTransform3d::translation(offset) * init_tr;

      for (size_t i = 0; i < from_pts.size(); ++i)
        from_pts[i] += offset;
    }

    KdTreeN<Vector3d, 3, double> to_kdtree(to_pts.begin(), to_pts.end());

    if (rotate_axis_aligned)
    {
      Icp3<double> icp(-1, -1, false);
      if (has_up_vector)
        icp.setUpVector(up_vector);

      Array<Vector3d> rot_from_pts(from_pts.size());
      double rot_error;
      bool first = true;

      for (int u = 0; u < 2; ++u)
        for (int su = -1; su <= 1; su += 2)
          for (int v = 1; v < 2; ++v)
            for (int sv = -1; sv <= 1; sv += 2)
            {
              Vector3d du(0, 0, 0); du[u] = su;
              Vector3d dv(0, 0, 0); dv[(u + v) % 3] = sv;
              Vector3d dw = du.cross(dv);
              Matrix3d rot; rot << du, dv, dw;

              AffineTransform3d rot_tr = AffineTransform3d::translation(from_center)
                                       * AffineTransform3d(rot)
                                       * AffineTransform3d::translation(-from_center);

              for (size_t i = 0; i < from_pts.size(); ++i)
                rot_from_pts[i] = rot_tr * from_pts[i];

              AffineTransform3d icp_tr = icp.align((intx)rot_from_pts.size(), &rot_from_pts[0], to_kdtree, &rot_error);
              rot_error = normalizeError(rot_error, to_scale, (intx)from_pts.size());
              if (first || rot_error < error)
              {
                error = rot_error;
                tr = icp_tr * rot_tr;
                first = false;

                THEA_CONSOLE << "-- rotation " << toString(rot) << " reduced error to " << error;
              }
              else
                THEA_CONSOLE << "-- rotation " << toString(rot) << ", no reduction in error";
            }

      alwaysAssertM(!first, "No alignment found after trying axis-aligned rotations");

      tr = tr * init_tr;
    }
    else if (rotate_arbitrary)
    {
      if (!has_up_vector)
      {
        THEA_ERROR << "Testing rotations without a valid up vector not currently supported";
        return -1;
      }

      Icp3<double> icp(-1, -1, false);
      if (has_up_vector)
        icp.setUpVector(up_vector);

      Array<Vector3d> rot_from_pts(from_pts.size());
      double rot_error;
      bool first = true;

      static intx const DEFAULT_NUM_ROTATIONS = 16;
      intx nr = (num_rotations > 0 ? num_rotations : DEFAULT_NUM_ROTATIONS);
      for (intx i = 0; i < nr; ++i)
      {
        double angle = i * (Math::twoPi() / nr);
        Matrix3d rot = Math::rotationAxisAngle(up_vector, angle);
        AffineTransform3d rot_tr = AffineTransform3d::translation(from_center)
                                 * AffineTransform3d(rot)
                                 * AffineTransform3d::translation(-from_center);

        for (size_t j = 0; j < from_pts.size(); ++j)
          rot_from_pts[j] = rot_tr * from_pts[j];

        AffineTransform3d icp_tr = icp.align((intx)rot_from_pts.size(), &rot_from_pts[0], to_kdtree, &rot_error);
        rot_error = normalizeError(rot_error, to_scale, (intx)from_pts.size());
        if (first || rot_error < error)
        {
          error = rot_error;
          tr = icp_tr * rot_tr;
          first = false;

          THEA_CONSOLE << "-- rotation by " << Math::radiansToDegrees(angle) << " degrees reduced error to " << error;
        }
        else
          THEA_CONSOLE << "-- rotation by " << Math::radiansToDegrees(angle) << " degrees, no reduction in error";
      }

      alwaysAssertM(!first, "No alignment found after trying axis-aligned rotations");

      tr = tr * init_tr;
    }
    else
    {
      Icp3<double> icp(-1, -1, true);
      if (has_up_vector)
        icp.setUpVector(up_vector);

      tr = icp.align((intx)from_pts.size(), &from_pts[0], to_kdtree, &error);
      tr = tr * init_tr;

      error = normalizeError(error, to_scale, (intx)from_pts.size());
    }
  }

  THEA_CONSOLE << "Alignment error = " << std::setprecision(10) << error;
  THEA_CONSOLE << "Alignment = " << tr.toString();

  if (out)
  {
    Matrix3d const & lin = tr.getLinear();
    Vector3d const & trn = tr.getTranslation();

    *out << from_path << '\n'
         << to_path << '\n'
         << error << '\n'
         << lin(0, 0) << ' ' << lin(0, 1) << ' ' << lin(0, 2) << ' ' << trn[0] << '\n'
         << lin(1, 0) << ' ' << lin(1, 1) << ' ' << lin(1, 2) << ' ' << trn[1] << '\n'
         << lin(2, 0) << ' ' << lin(2, 1) << ' ' << lin(2, 2) << ' ' << trn[2] << endl;
  }

  return 0;
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
      else if (arg == "-t")
      {
        if (i >= argc - 3)
          return usage(argc, argv);

        do_initial_placement = true;
        for (int j = 0; j < 3; ++j)
        {
          istringstream iss(argv[++i]);
          if (!(iss >> initial_from_position[j]))
          {
            THEA_ERROR << "Could not parse coordinate " << j << " of initial location";
            return -1;
          }
        }

        THEA_CONSOLE << "Source shape will be initially placed at " << initial_from_position;
      }
      else if (arg == "-u")
      {
        if (i >= argc - 3)
          return usage(argc, argv);

        has_up_vector = true;
        for (int j = 0; j < 3; ++j)
        {
          istringstream iss(argv[++i]);
          if (!(iss >> up_vector[j]))
          {
            THEA_ERROR << "Could not parse coordinate " << j << " of up vector";
            return -1;
          }
        }

        THEA_CONSOLE << "Up vector is " << toString(up_vector);
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
      else if (arg == "-r")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        istringstream iss(argv[++i]);
        if (!(iss >> num_rotations))
        {
          THEA_ERROR << "Could not parse number of rotations sampled around up vector";
          return -1;
        }

        rotate_arbitrary = true;
        if (num_rotations <= 0)
          THEA_CONSOLE << "Alignment will search over rotations around up vector";
        else
          THEA_CONSOLE << "Alignment will search over " << num_rotations << " rotation(s) around up vector";
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

  if (prealign_centroids && do_initial_placement)
  {
    THEA_ERROR << "Cannot prealign centroids AND do initial placement";
    return -1;
  }

  ostream * out = nullptr;
  ofstream fout;
  if (!out_path.empty())
  {
    fout.open(out_path.c_str());
    if (!fout)
    {
      THEA_ERROR << "Couldn't open output file " << out_path;
      return -1;
    }

    out = &fout;
  }

  Array<string> from_shape_paths, to_shape_paths;
  if (!loadShapePaths(from_path, from_shape_paths) || !loadShapePaths(to_path, to_shape_paths))
    return -1;

  for (size_t i = 0; i < from_shape_paths.size(); ++i)
    for (size_t j = 0; j < to_shape_paths.size(); ++j)
    {
      if (!(i == 0 && j == 0) && out)
        *out << '\n';

      THEA_CONSOLE << "*** Aligning " << from_shape_paths[i] << " to " << to_shape_paths[j] << " ***";
      int status = alignShapes(from_shape_paths[i], to_shape_paths[j], out);
      if (status != 0)
        return -1;
    }

  return 0;
}
