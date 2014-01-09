#include "../../Common.hpp"
#include "../../Algorithms/IntersectionTester.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/RayIntersectionTester.hpp"
#include "../../Algorithms/ShortestPaths.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Array.hpp"
#include "../../Ball3.hpp"
#include "../../BoundedSortedArray.hpp"
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

int max_nbrs = 8;

class SurfaceSample;

struct Neighbor
{
  SurfaceSample * sample;
  Real separation;

  Neighbor(SurfaceSample * sample_ = NULL, Real separation_ = 0) : sample(sample_), separation(separation_) {}

  bool operator<(Neighbor const & rhs) const
  {
    return separation < rhs.separation || (separation == rhs.separation && std::less<SurfaceSample *>()(sample, rhs.sample));
  }

}; // struct Neighbor

namespace boost {

template <> struct has_trivial_assign<Neighbor> : public boost::true_type {};

} // namespace boost

class SurfaceSample
{
  public:
    typedef BoundedSortedArray<Neighbor> NeighborSet;

    SurfaceSample() : nbrs(max_nbrs) {}

    Vector3 p, n;
    NeighborSet nbrs;
};

namespace Thea {
namespace Algorithms {

template <>
class IsPointN< ::SurfaceSample, 3 >
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN< ::SurfaceSample, 3 >::getPosition(::SurfaceSample const & sample)
{
  return sample.p;
}

} // namespace Algorithms
} // namespace Thea

template <typename RefIterator>
struct RefToPtrIterator : public RefIterator
{
  RefToPtrIterator(RefIterator i) : RefIterator(i) {}

  typename RefIterator::value_type * operator*() const { return &(this->RefIterator::operator*()); }

}; // struct RefToPtrIterator

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;
typedef MeshKDTree<Mesh> KDTree;
typedef KDTreeN<SurfaceSample *, 3> SampleKDTree;

TheaArray<SurfaceSample> samples;
array_size_t orig_num_samples = 0;
KDTree kdtree;
SampleKDTree sample_kdtree;
Real sample_nbd_radius;

void computeSampleNeighborhoodRadius();
void computeSampleNeighbors(SurfaceSample & sample);
void extractOriginalAdjacencies();

int
main(int argc, char * argv[])
{
  string mesh_path;
  string samples_path;
  string out_path;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
    {
      switch (curr_opt)
      {
        case 0: mesh_path = argv[i]; break;
        case 1: samples_path = argv[i]; break;
        case 2: out_path = argv[i]; break;
      }

      curr_opt++;
      if (curr_opt >= 3)
        break;
    }
    else
    {
      if (beginsWith(argv[i], "--max-nbrs="))
      {
        if (!sscanf(argv[i], "--max-nbrs=%d", &max_nbrs) == 1 || max_nbrs <= 0)
        {
          THEA_ERROR << "Invalid --max-nbrs parameter";
          return -1;
        }
      }
      else
      {
        THEA_ERROR << "Unknown parameter: " << argv[i];
        return -1;
      }
    }
  }

  if (curr_opt < 3)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " [--max-nbrs=N] <mesh> <points> <graph-outfile>";
    return 0;
  }

  //===========================================================================================================================
  // Load mesh
  //===========================================================================================================================

  MG mg;
  try
  {
    mg.load(mesh_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load mesh %s", mesh_path.c_str())

  //===========================================================================================================================
  // Load points
  //===========================================================================================================================

  ifstream in(samples_path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not load points from file " << samples_path;
    return -1;
  }

  THEA_CONSOLE << "Loaded mesh from " << mesh_path;

  SurfaceSample s;
  string line;
  long line_num = 0;
  bool has_normals = false;
  while (getline(in, line))
  {
    line_num++;

    boost::trim(line);
    if (line.empty())
      continue;

    istringstream line_in(line);
    if (!(line_in >> s.p[0] >> s.p[1] >> s.p[2]))
    {
      THEA_ERROR << "Could not read point on line " << line_num << " of file " << samples_path;
      return -1;
    }

    if (samples.empty() && (line_in >> s.n[0] >> s.n[1] >> s.n[2]))
      has_normals = true;
    else if (has_normals)
    {
      if (!(line_in >> s.n[0] >> s.n[1] >> s.n[2]))
      {
        THEA_ERROR << "Could not read normal on line " << line_num << " of file " << samples_path;
        return -1;
      }
    }

    samples.push_back(s);
  }

  THEA_CONSOLE << "Loaded " << samples.size() << " samples(s) from " << samples_path;

  //===========================================================================================================================
  // KD-tree on mesh
  //===========================================================================================================================

  kdtree.add(mg);
  kdtree.init();

  //===========================================================================================================================
  // If the samples have no normals, compute them
  //===========================================================================================================================

  if (!has_normals)
  {
    for (array_size_t i = 0; i < samples.size(); ++i)
    {
      long elem = kdtree.closestElement<MetricL2>(samples[i].p);
      if (elem < 0)
      {
        THEA_ERROR << "Could not find nearest neighbor of sample " << i << " on mesh";
        return false;
      }

      samples[i].n = kdtree.getElements()[(array_size_t)elem].getNormal();
    }

    THEA_CONSOLE << "Computed sample normals";
  }

  //===========================================================================================================================
  // Augment the number of samples if necessary
  //===========================================================================================================================

  static array_size_t const MIN_SAMPLES = 10000;

  orig_num_samples = samples.size();
  bool added_extra_samples = false;
  if (orig_num_samples < MIN_SAMPLES)
  {
    TheaArray<Vector3> positions;
    TheaArray<Vector3> face_normals;

    MeshSampler<Mesh> sampler(mg);
    sampler.sampleEvenlyByArea((long)(MIN_SAMPLES - orig_num_samples), positions, &face_normals);

    for (array_size_t i = 0; i < positions.size(); ++i)
    {
      s.p = positions[i];
      s.n = face_normals[i];
      samples.push_back(s);
    }

    added_extra_samples = (samples.size() > orig_num_samples);
    if (added_extra_samples)
      THEA_CONSOLE << samples.size() - orig_num_samples << " extra samples added to original set, for density";
  }

  //===========================================================================================================================
  // KD-tree on samples
  //===========================================================================================================================

  sample_kdtree.init(RefToPtrIterator< TheaArray<SurfaceSample>::iterator >(samples.begin()),
                     RefToPtrIterator< TheaArray<SurfaceSample>::iterator >(samples.end()));

  //===========================================================================================================================
  // Compute neighbors
  //===========================================================================================================================

  computeSampleNeighborhoodRadius();

  for (array_size_t i = 0; i < samples.size(); ++i)
    computeSampleNeighbors(samples[i]);

  THEA_CONSOLE << "Computed sample adjacencies";

  //===========================================================================================================================
  // If we added extra samples, we must recompute neighbors only among the original samples
  //===========================================================================================================================

  if (added_extra_samples)
  {
    extractOriginalAdjacencies();

    THEA_CONSOLE << "Adjacencies computed on original samples using extra samples as intermediates";
  }

  //===========================================================================================================================
  // Write graph to file
  //===========================================================================================================================

  ofstream out(out_path.c_str());
  if (!out)
  {
    THEA_ERROR << "Could not open output file " << out_path << " for writing features";
    return -1;
  }

  double sum_degrees = 0;
  for (array_size_t i = 0; i < orig_num_samples; ++i)
  {
    out << samples[i].nbrs.size();

    for (int j = 0; j < samples[i].nbrs.size(); ++j)
    {
      size_t index = samples[i].nbrs[j].sample - &samples[0];
      alwaysAssertM(index < orig_num_samples, "Neighbor index out of bounds");

      out << ' ' << index;
    }

    out << '\n';

    sum_degrees += samples[i].nbrs.size();
  }

  THEA_CONSOLE << "Wrote sample graph of average degree " << sum_degrees / orig_num_samples << " to " << out_path;

  return 0;
}

void
computeSampleNeighborhoodRadius()
{
  static Real const NBD_RADIUS_FACTOR = 1.0f;

  if (samples.empty())
    sample_nbd_radius = 0;
  else
  {
    double total_area = 0;
    for (long i = 0; i < kdtree.numElements(); ++i)
      total_area += kdtree.getElements()[i].getArea();

    sample_nbd_radius = NBD_RADIUS_FACTOR * Math::fastSqrt(total_area / samples.size());

    // THEA_CONSOLE << "sample_nbd_radius = " << sample_nbd_radius << ", bounding box diag = "
    //              << kdtree.getBounds().getExtent().length();
  }
}

bool
validNeighbors(SurfaceSample const * sample1, SurfaceSample const * sample2, Real & sep)
{
  // Tests are (hopefully) in decreasing order of speed

  // Identity test
  if (sample1 == sample2)
    return false;

  // Angle test
  static Real const MIN_DOT = -0.5f;
  if (sample1->n.dot(sample2->n) < MIN_DOT)
    return false;

  // Duplication test
  for (int i = 0; i < sample1->nbrs.size(); ++i)
    if (sample2 == sample1->nbrs[i].sample)
      return false;

  // Compute the separation of the samples
  Vector3 diff = sample2->p - sample1->p;
  Real sqsep = diff.squaredLength();
  sep = Math::fastSqrt(sqsep);

#if 1
  // Reachability test: lift points off the surface by an amount proportional (currently, equal) to their separation, and
  // see if the line connecting them is blocked by the surface
  static Real const LIFT_FACTOR = 5;
  Vector3 lift_dir = (sample1->n + sample2->n).fastUnit();
  Ray3 ray(sample1->p + LIFT_FACTOR * sep * lift_dir, diff);
  if (kdtree.rayIntersects<RayIntersectionTester>(ray, 1))
    return false;
#endif

  return true;
}

class NeighborFunctor
{
  public:
    NeighborFunctor(SurfaceSample * sample_)
    : sample(sample_)
    {}

    bool operator()(long index, SurfaceSample * nbr)
    {
      Real sep = 0;
      if (validNeighbors(sample, nbr, sep))
        sample->nbrs.insert(Neighbor(nbr, sep));

      return false;
    }

  private:
    SurfaceSample * sample;

}; // struct NeighborFunctor

// #define TIME_UPDATE_SAMPLE_NEIGHBORS

void
computeSampleNeighbors(SurfaceSample & sample)
{
#ifdef TIME_UPDATE_SAMPLE_NEIGHBORS
  StopWatch timer;
  timer.tick();
#endif

  static int const MAX_ITERS = 3;
  static Real const RADIUS_EXPANSION_FACTOR = 2.0f;
  static int min_nbrs = max(4, (int)(0.25 * max_nbrs));

  Real radius = sample_nbd_radius;
  for (int i = 0; i < MAX_ITERS; ++i)
  {
    sample.nbrs.clear();
    NeighborFunctor func(&sample);

    Ball3 nbd(sample.p, radius);
    sample_kdtree.processRangeUntil<IntersectionTester>(nbd, &func);

    if (sample.nbrs.size() >= min_nbrs)
      break;
    else
      radius *= RADIUS_EXPANSION_FACTOR;
  }

#ifdef TIME_UPDATE_SAMPLE_NEIGHBORS
  timer.tock();
  THEA_CONSOLE << getName() << ": Neighbors of sample " << &sample << " updated in " << 1000 * timer.elapsedTime()
               << "ms, sample has " << sample.nbrs.size() << " neighbors";
#endif
}

struct Graph
{
  typedef TheaArray<SurfaceSample> NodeArray;
  typedef SurfaceSample::NeighborSet NeighborSet;

  NodeArray * nodes;

  typedef SurfaceSample * VertexHandle;
  typedef SurfaceSample const * VertexConstHandle;
  typedef NodeArray::iterator VertexIterator;
  typedef NodeArray::const_iterator VertexConstIterator;

  typedef Neighbor * NeighborIterator;
  typedef Neighbor const * NeighborConstIterator;

  Graph(NodeArray * nodes_) : nodes(nodes_) {}

  long numVertices() const { return (long)nodes->size(); }

  VertexIterator verticesBegin() { return nodes->begin(); }
  VertexConstIterator verticesBegin() const { return nodes->begin(); }

  VertexIterator verticesEnd() { return nodes->end(); }
  VertexConstIterator verticesEnd() const { return nodes->end(); }

  VertexHandle getVertex(VertexIterator vi) { return &(*vi); }
  VertexConstHandle getVertex(VertexConstIterator vi) const { return &(*vi); }

  long numNeighbors(VertexConstHandle vertex) const { return vertex->nbrs.size(); }

  NeighborIterator neighborsBegin(VertexHandle vertex)
  {
    return vertex->nbrs.isEmpty() ? NULL : const_cast<Neighbor *>(&vertex->nbrs[0]);
  }

  NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
  {
    return vertex->nbrs.isEmpty() ? NULL : &vertex->nbrs[0];
  }

  NeighborIterator neighborsEnd(VertexHandle vertex) { return neighborsBegin(vertex) + numNeighbors(vertex); }
  NeighborConstIterator neighborsEnd(VertexConstHandle vertex) const { return neighborsBegin(vertex) + numNeighbors(vertex); }

  VertexHandle getVertex(NeighborIterator ni) { return ni->sample; }
  VertexConstHandle getVertex(NeighborConstIterator ni) const { return ni->sample; }

  double distance(VertexConstHandle v, NeighborConstIterator ni) const { return ni->separation; }
};

struct DijkstraCallback
{
  DijkstraCallback(SurfaceSample * sample_, array_size_t src_index_) : sample(sample_), src_index(src_index_)
  {
    sample->nbrs.clear();
  }

  bool operator()(Graph::VertexHandle vertex, double distance, bool has_pred, Graph::VertexHandle pred)
  {
    array_size_t index = vertex - &samples[0];
    if (index < orig_num_samples && index != src_index)
      sample->nbrs.insert(Neighbor(vertex, (Real)distance));

    return false;
  }

  SurfaceSample * sample;
  array_size_t src_index;
};

void
extractOriginalAdjacencies()
{
  Graph graph(&samples);
  ShortestPaths<Graph> shortest_paths;

  TheaArray<SurfaceSample> samples_with_new_nbrs(orig_num_samples);
  for (array_size_t i = 0; i < orig_num_samples; ++i)
  {
    samples_with_new_nbrs[i] = samples[i];
    DijkstraCallback callback(&samples_with_new_nbrs[i], i);
    shortest_paths.dijkstra(graph, &samples[i], &callback);
  }

  for (array_size_t i = 0; i < orig_num_samples; ++i)
    samples[i] = samples_with_new_nbrs[i];
}
