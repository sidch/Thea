#include "../../Common.hpp"
#include "../../Algorithms/MeshBvh.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/PointSet3.hpp"
#include "../../Algorithms/ShortestPaths.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../AffineTransform3.hpp"
#include "../../Array.hpp"
#include "../../FilePath.hpp"
#include "../../MatVec.hpp"
#include "../../Noncopyable.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Graphics;
using namespace Algorithms;

int max_nbrs = 8;
long min_samples = 50000;
bool consistent_normals = false;
bool reachability = false;
bool pairwise_distances = false;

enum { LOAD_ERROR = 1, PARSE_ERROR, UNSUPPORTED_FORMAT };

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;
typedef MeshBvh<Mesh> Bvh;

int loadSamples(string const & samples_path, Array<Vector3> & positions, Array<Vector3> & normals);

struct MeshTransformer
{
  MeshTransformer(AffineTransform3 const & tr_) : tr(tr_) {}
  bool operator()(Mesh & mesh) const
  {
    for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      vi->setPosition(tr * vi->getPosition());

    return false;
  }

  AffineTransform3 tr;
};

struct DistanceCallback : public Noncopyable
{
  DistanceCallback(intx n) : m(n, n), current_source(-1) { m.fill(-1); }

  bool operator()(PointSet3::SampleGraph::VertexHandle vertex, double distance, bool has_pred,
                  PointSet3::SampleGraph::VertexHandle pred)
  {
    if (!vertex)
      return false;

    intx i = vertex->getIndex();
    m(current_source, i) = m(i, current_source) = distance;

    return false;
  }

  MatrixX<double> m;
  intx current_source;
};

int
main(int argc, char * argv[])
{
  string mesh_path;
  string samples_path;
  string out_path;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (!beginsWith(arg, "-"))
    {
      switch (curr_opt)
      {
        case 0: mesh_path = arg; break;
        case 1: samples_path = arg; break;
        case 2: out_path = arg; break;
      }

      curr_opt++;
      if (curr_opt > 3)
        break;
    }
    else
    {
      if (beginsWith(arg, "--max-nbrs="))
      {
        if (sscanf(arg.c_str(), "--max-nbrs=%d", &max_nbrs) != 1 || max_nbrs <= 0)
        {
          THEA_ERROR << "Invalid --max-nbrs parameter";
          return -1;
        }
      }
      else if (beginsWith(arg, "--min-samples="))
      {
        if (sscanf(arg.c_str(), "--min-samples=%ld", &min_samples) != 1 || min_samples <= 0)
        {
          THEA_ERROR << "Invalid --min-samples parameter";
          return -1;
        }
      }
      else if (arg == "-n" || arg == "--normals")
      {
        consistent_normals = true;
      }
      else if (arg == "-r" || arg == "--reachability")
      {
        reachability = true;
      }
      else if (arg == "-d" || arg == "--distances")
      {
        pairwise_distances = true;
      }
      else
      {
        THEA_ERROR << "Unknown parameter: " << arg;
        return -1;
      }
    }
  }

  if (curr_opt != 3)
  {
    THEA_CONSOLE << "";
    THEA_CONSOLE << "Usage: " << argv[0] << " [OPTIONS] <mesh|dense-points> <points> <surf-file>";
    THEA_CONSOLE << "";
    THEA_CONSOLE << "Options:";
    THEA_CONSOLE << "  --max-nbrs=N          Maximum degree of proximity surf";
    THEA_CONSOLE << "  --min-samples=N       Minimum number of original plus generated samples";
    THEA_CONSOLE << "  --normals | -n        Run extra tests assuming consistently oriented mesh normals";
    THEA_CONSOLE << "  --reachability | -r   Reachability test for adjacency (requires -n)";
    THEA_CONSOLE << "  --distances | -d      Output distances b/w all pairs of points as a .dist file";
    THEA_CONSOLE << "                        Requires surf to already exist";
    THEA_CONSOLE << "";

    return -1;
  }

  if (reachability && !consistent_normals)
  {
    THEA_ERROR << "Reachability test requires consistent normals";
    return -1;
  }

  //===========================================================================================================================
  // Load sample graph if it already exists, compute and write all pairwise distances, and quit
  //===========================================================================================================================

  if (pairwise_distances)
  {
    PointSet3 surf;
    if (!surf.load(samples_path, out_path))
    {
      THEA_CONSOLE << "Could not load graph from file " << out_path;
      return -1;
    }

    auto const & samples = surf.getSamples();
    DistanceCallback distance_callback((intx)samples.size());
    ShortestPaths<PointSet3::SampleGraph> shortest_paths;

    distance_callback.m(0, 0) = 0;
    for (size_t i = 1; i < samples.size(); ++i)  // matrix is symmetric so no need to have 0 as source
    {
      distance_callback.current_source = (intx)i;
      shortest_paths.dijkstraWithCallback(const_cast<PointSet3::SampleGraph &>(surf.getGraph()),
                                          const_cast<PointSet3::SampleGraph::VertexHandle>(&samples[i]),
                                          std::ref(distance_callback));
    }

    ofstream d_out(FilePath::changeExtension(out_path, "dist").c_str());
    for (intx r = 0; r < distance_callback.m.rows(); ++r)
    {
      for (intx c = 0; c < distance_callback.m.cols(); ++c)
      {
        if (c > 0) d_out << ' ';
        d_out << distance_callback.m(r, c);
      }

      d_out << endl;
    }

    return 0;
  }

  //===========================================================================================================================
  // Load points
  //===========================================================================================================================

  Array<Vector3> sample_positions;
  Array<Vector3> sample_normals;
  if (loadSamples(samples_path, sample_positions, sample_normals))
    return -1;

  bool has_normals = (!sample_positions.empty() && sample_normals.size() == sample_positions.size());

  THEA_CONSOLE << "Loaded " << sample_positions.size() << " sample(s) from " << samples_path;

  //===========================================================================================================================
  // Load mesh or dense samples
  //===========================================================================================================================

  Array<Vector3> dense_positions;
  Array<Vector3> dense_normals;
  bool dense_has_normals = false;
  Bvh bvh;
  if (mesh_path != "-")
  {
    // First try to load the file as a set of points
    int dense_load_status = loadSamples(mesh_path, dense_positions, dense_normals);
    dense_has_normals = (!dense_positions.empty() && dense_normals.size() == dense_positions.size());

    if (dense_load_status != 0 && dense_load_status != UNSUPPORTED_FORMAT)
    {
      return -1;
    }
    else if (dense_load_status == 0)
    {
      THEA_CONSOLE << dense_positions.size() << " extra samples added from: " << mesh_path;
    }
    else  // Now try to load the file as a mesh, if the above failed
    {
      MG mg;
      try
      {
        mg.load(mesh_path);
      }
      THEA_CATCH(return -1;, ERROR, "Could not load mesh %s", mesh_path.c_str())

      THEA_CONSOLE << "Loaded mesh from " << mesh_path;

      // Make sure the mesh is properly scaled
      AxisAlignedBox3 mesh_bounds = mg.getBounds();
      AxisAlignedBox3 samples_bounds;
      for (size_t i = 0; i < sample_positions.size(); ++i)
        samples_bounds.merge(sample_positions[i]);

      Real scale_error = (mesh_bounds.getLow()  - samples_bounds.getLow()).norm()
                       + (mesh_bounds.getHigh() - samples_bounds.getHigh()).norm();
      if (scale_error > 0.01 * mesh_bounds.getExtent().norm())
      {
        // Rescale the mesh.
        Real scale = (samples_bounds.getExtent().cwiseQuotient(mesh_bounds.getExtent())).maxCoeff();  // samples give a smaller
                                                                                     // bound than the  true bound, so take the
                                                                                     // axis in which the approximation is best

        AffineTransform3 tr = AffineTransform3::translation(samples_bounds.getCenter())
                            * AffineTransform3::scaling(scale)
                            * AffineTransform3::translation(-mesh_bounds.getCenter());
        MeshTransformer func(tr);
        mg.forEachMeshUntil(std::cref(func));
        mg.updateBounds();

        THEA_CONSOLE << "Matched scale of source mesh and original samples";
      }

      bvh.add(mg);
      bvh.init();

      // If the samples have no normals, compute them
      if (consistent_normals && !has_normals)
      {
        sample_normals.resize(sample_positions.size());

        for (size_t i = 0; i < sample_positions.size(); ++i)
        {
          intx elem = bvh.closestElement<MetricL2>(sample_positions[i]);
          if (elem < 0)
          {
            THEA_ERROR << "Could not find nearest neighbor of sample " << i << " on mesh";
            return -1;
          }

          sample_normals[i] = bvh.getElements()[(size_t)elem].getNormal();
        }

        has_normals = true;

        THEA_CONSOLE << "Computed sample normals";
      }

      // Augment the number of samples if necessary
      if ((intx)sample_positions.size() < min_samples)
      {
        dense_positions.clear();
        dense_normals.clear();

        MeshSampler<Mesh> sampler(mg);
        sampler.sampleEvenlyByArea((intx)(min_samples - sample_positions.size()), dense_positions,
                                   (consistent_normals ? &dense_normals : nullptr));

        if (!dense_positions.empty())
          THEA_CONSOLE << dense_positions.size() << " extra samples added to original set, for density";

        dense_has_normals = (!dense_positions.empty() && dense_normals.size() == dense_positions.size());
      }
    }
  }

  if (consistent_normals && !has_normals)
  {
    THEA_ERROR << "Samples loaded from " << samples_path << " do not have normals for consistency test";
    return -1;
  }

  if (consistent_normals && !dense_positions.empty() && !dense_has_normals)
  {
    THEA_ERROR << "Dense samples loaded from " << mesh_path << " do not have normals for consistency test";
    return -1;
  }

  //===========================================================================================================================
  // Create adjacency graph on samples
  //===========================================================================================================================

  PointSet3::Options opts;
  opts.setMaxDegree(max_nbrs);
  PointSet3 surf(opts);
  surf.addSamples((intx)sample_positions.size(), &sample_positions[0], (consistent_normals ? &sample_normals[0] : nullptr));
  surf.addOversampling((intx)dense_positions.size(), &dense_positions[0], (consistent_normals ? &dense_normals[0] : nullptr));
  surf.updateGraph(reachability && !bvh.empty() ? &bvh : nullptr);

  THEA_CONSOLE << "Computed sample graph";

  //===========================================================================================================================
  // Write surf to file
  //===========================================================================================================================

  if (!surf.save(/* samples_path = */ "", out_path, true))
    return -1;

  double sum_degrees = 0;
  PointSet3::SampleArray const & samples = surf.getSamples();
  for (size_t i = 0; i < samples.size(); ++i)
    sum_degrees += samples[i].getNeighbors().size();

  THEA_CONSOLE << "Wrote sample graph of average degree " << sum_degrees / samples.size() << " to " << out_path;

  return 0;
}

int
loadSamples(string const & samples_path, Array<Vector3> & positions, Array<Vector3> & normals)
{
  ifstream in(samples_path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not load points from file " << samples_path;
    return LOAD_ERROR;
  }

  positions.clear();
  normals.clear();

  Vector3 p, n;
  string path_lc = toLower(samples_path);
  bool has_normals = false;
  if (endsWith(path_lc, ".pts"))
  {
    string line;
    intx line_num = 0;
    while (getline(in, line))
    {
      line_num++;

      line = trimWhitespace(line);
      if (line.empty())
        continue;

      istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2]))
      {
        THEA_ERROR << "Could not read point on line " << line_num << " of file " << samples_path;
        return PARSE_ERROR;
      }

      positions.push_back(p);

      if (positions.size() == 1)
      {
        if (line_in >> n[0] >> n[1] >> n[2])
        {
          has_normals = true;
          normals.push_back(n);
        }
      }
      else if (has_normals)
      {
        if (!(line_in >> n[0] >> n[1] >> n[2]))
        {
          THEA_ERROR << "Could not read normal on line " << line_num << " of file " << samples_path;
          return PARSE_ERROR;
        }

        normals.push_back(n);
      }
    }
  }
  else if (endsWith(path_lc, ".off"))
  {
    string magic;
    in >> magic;
    if (magic != "OFF")
    {
      THEA_ERROR << "Could not read OFF header from file " << samples_path;
      return PARSE_ERROR;
    }

    intx nv, nf, ne;
    if (!(in >> nv >> nf >> ne))
    {
      THEA_ERROR << "Could not read element counts from OFF file " << samples_path;
      return PARSE_ERROR;
    }

    if (nf == 0)
    {
      for (intx i = 0; i < nv; ++i)
      {
        if (!(in >> p[0] >> p[1] >> p[2]))
        {
          THEA_ERROR << "Could not parse point " << i << " from OFF file " << samples_path;
          return PARSE_ERROR;
        }

        positions.push_back(p);
      }
    }
    else
      return UNSUPPORTED_FORMAT;
  }
  else
    return UNSUPPORTED_FORMAT;

  return 0;
}
