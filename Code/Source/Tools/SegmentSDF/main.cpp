#include "Graph.hpp"

#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/Clustering.hpp"
#include "../../Algorithms/ConvexHull3.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Colors.hpp"
#include "../../Math.hpp"
#include "../../System.hpp"
#include "../../UnionFind.hpp"
#include "../../UnorderedMap.hpp"
#include "../../VectorN.hpp"
#include <boost/functional/hash.hpp>
#include <boost/thread/thread.hpp>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

// If defined, map the label to a color and print it as a fake normal after the position, instead of the SDF and label
#define DEBUG_PTS

using namespace std;
using namespace Thea;
using namespace Graphics;
using namespace Algorithms;

int segmentSDF(int argc, char * argv[]);

int
main(int argc, char * argv[])
{
  int status = 0;
  try
  {
    status = segmentSDF(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return status;
}

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;
typedef MeshKDTree<Mesh> KDTree;
typedef VectorN<4, double> ClusterablePoint;

ClusterablePoint
toClusterablePoint(Vector3 const & pos, Real sdf)
{
  ClusterablePoint cp;
  cp[0] = pos[0];
  cp[1] = pos[1];
  cp[2] = pos[2];

#define SDF_WEIGHT 3
  cp[3] = SDF_WEIGHT * sdf;

  return cp;
}

int countSDFModes(TheaArray<Real> const & sdf_values);
ColorRGB const & getPaletteColor(int i);
int combineClustersByConvexity(TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals, KDTree const & kdtree,
                               TheaArray<int> & labels, double concavity_threshold, double score_threshold);

int
segmentSDF(int argc, char * argv[])
{
  if (argc < 3)
  {
    cout << "Usage: " << argv[0] << " <mesh-file> <out-file> [<approx-num-samples>]" << endl;
    return 0;
  }

#ifdef THEA_OSX
  THEA_WARNING << "!!! On OS X, this program will only print the guessed number of segments !!!";
  THEA_WARNING << "!!! No other output will be produced !!!";
#endif

  string inpath = argv[1];
  string outpath = argv[2];
  long approx_num_samples = 5000;
  if (argc > 3)
  {
    approx_num_samples = atoi(argv[3]);
    if (approx_num_samples < 1)
    {
      THEA_ERROR << "Invalid number of samples: " << approx_num_samples;
      return -1;
    }
  }

  // Load mesh
  MG mg;
  mg.load(inpath);

  // Initialize kdtree
  KDTree kdtree;
  kdtree.add(mg);
  kdtree.init();

  // Compute samples
  double total_area = 0;
  for (long i = 0; i < kdtree.numElements(); ++i)
    total_area += kdtree.getElements()[(size_t)i].getArea();

  double density = approx_num_samples / total_area;

  TheaArray<Vector3> positions, normals;
  for (long i = 0; i < kdtree.numElements(); ++i)
  {
    KDTree::Element const & elem = kdtree.getElements()[(size_t)i];
    double num_samples = density * elem.getArea();
    double rem = num_samples;
    for (long j = 1; j < num_samples; ++j)
    {
      positions.push_back(elem.randomPoint());
      normals.push_back(elem.getNormal());
      rem -= 1;
    }

    if (Random::common().uniform01() < rem)
    {
      positions.push_back(elem.randomPoint());
      normals.push_back(elem.getNormal());
    }
  }

  THEA_CONSOLE << "Computed " << positions.size() << " sample points on the mesh";

  // Compute SDF values
  TheaArray<Real> sdf_values(positions.size());
  MeshFeatures::Local::ShapeDiameter<Mesh> sdf(&kdtree);
  for (size_t i = 0; i < positions.size(); ++i)
    sdf_values[i] = (Real)sdf.compute(positions[i], normals[i]);

  // Undo normalization
  Real scale = sdf.getNormalizationScale();
  for (size_t i = 0; i < sdf_values.size(); ++i)
    sdf_values[i] *= scale;

  THEA_CONSOLE << "Computed SDF values";

  // Estimate number of clusters
  int num_clusters_hint = countSDFModes(sdf_values);
  if (num_clusters_hint <= 2)
    num_clusters_hint = 7;
  else if (num_clusters_hint > 20)
    num_clusters_hint = 20;

  THEA_CONSOLE << "Estimated " << num_clusters_hint << " modes in the SDF distribution";

#ifndef THEA_OSX

  TheaArray<ClusterablePoint> clusterable_pts(positions.size());
  for (size_t i = 0; i < positions.size(); ++i)
    clusterable_pts[i] = toClusterablePoint(positions[i], sdf_values[i]);

  // Set clustering options (try G1, SLINK_W, UPGMA, CUT)
  Clustering::FlatOptions options(Clustering::FlatMethod::CLUTO_GRAPH_CLUSTER_RB,
                                  Clustering::SimilarityMeasure::L2_DISTANCE,
                                  Clustering::GraphModel::SYMMETRIC_DIRECT, 40,
                                  Clustering::ClusteringCriterion::CLUTO_CUT,
                                  Clustering::SplitPriority::BEST_FIRST);

  // Do the clustering
  TheaArray<int> labels(positions.size());
  int num_clusters = Clustering::computeFlat(clusterable_pts, labels, options, num_clusters_hint);

  THEA_CONSOLE << "Segmented points into " << num_clusters << " clusters";

  // Merge clusters whose union is approximately convex
  num_clusters = combineClustersByConvexity(positions, normals, kdtree, labels, -1, -1);

  THEA_CONSOLE << "Merged clusters into " << num_clusters << " clusters";

  // Write out the labels
  ofstream out(outpath.c_str());
  out << positions.size() << endl;
  out << num_clusters << endl;
  for (size_t i = 0; i < positions.size(); ++i)
  {
#ifdef DEBUG_PTS
    ColorRGB pseudo_n = getPaletteColor(labels[i]);
    // Real ns = sdf_values[i] / sdf.getNormalizationScale();
    // ColorRGB pseudo_n(ns, sqrt(1 - ns * ns), 0);
    out << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << ' '
        << pseudo_n.r()    << ' ' << pseudo_n.g()    << ' ' << pseudo_n.b() << ' ' << sdf_values[i] << ' ' << labels[i] << endl;
#else
    out << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << ' ' << sdf_values[i] << ' ' << labels[i]
        << endl;
#endif
  }

#endif

  return 0;
}

static int const NUM_PALETTE_COLORS = 24;
ColorRGB COLOR_PALETTE[NUM_PALETTE_COLORS] = {
  ColorRGB::fromARGB(0x298edb),
  ColorRGB::fromARGB(0x982411),
  ColorRGB::fromARGB(0x6d4e25),
  ColorRGB::fromARGB(0x1b5043),
  ColorRGB::fromARGB(0x6e7662),
  ColorRGB::fromARGB(0xa08b00),
  ColorRGB::fromARGB(0x58427b),
  ColorRGB::fromARGB(0x1d2f5b),
  ColorRGB::fromARGB(0xac5e34),
  ColorRGB::fromARGB(0x804055),
  ColorRGB::fromARGB(0x6d7a00),
  ColorRGB::fromARGB(0x572e2c),

  // Invert each color above
  ColorRGB::fromARGB(~0x298edb),
  ColorRGB::fromARGB(~0x982411),
  ColorRGB::fromARGB(~0x6d4e25),
  ColorRGB::fromARGB(~0x1b5043),
  ColorRGB::fromARGB(~0x6e7662),
  ColorRGB::fromARGB(~0xa08b00),
  ColorRGB::fromARGB(~0x58427b),
  ColorRGB::fromARGB(~0x1d2f5b),
  ColorRGB::fromARGB(~0xac5e34),
  ColorRGB::fromARGB(~0x804055),
  ColorRGB::fromARGB(~0x6d7a00),
  ColorRGB::fromARGB(~0x572e2c)
};

int
numPaletteColors()
{
  return NUM_PALETTE_COLORS;
}

ColorRGB const &
getPaletteColor(int i)
{
  return COLOR_PALETTE[i % numPaletteColors()];
}

namespace MathInternal {

// exp(-x) approximation from http://www.xnoiz.co.cc/fast-exp-x/
inline float fastMinuzExp_1(float x) {

        // err <= 3e-3
        return 1
                -x*(0.9664
                -x*(0.3536));

}

} // namespace MathInternal

/** Computes a fast approximation to exp(-x) for 0 <= x <= 1. */
inline float
fastMinusExp01(float x)
{
  if (x > 0.69)
  {
    float y = MathInternal::fastMinuzExp_1(0.5 * x);
    return y * y;
  }
  else
    return MathInternal::fastMinuzExp_1(x);
}

// Returns a value proportional to the (negative of the) derivative of k(s) w.r.t. s, where the kernel is k(x^2).
inline double
dKernel(double x_squared)
{
  return fastMinusExp01(x_squared);
//   return 1;
}

/** A value, plus an index, sortable by value. */
struct IndexedValue
{
  IndexedValue() {}
  IndexedValue(double value_, int index_) : value(value_), index(index_) {}

  bool operator<(IndexedValue const & rhs) const { return value < rhs.value; }
  bool operator>(IndexedValue const & rhs) const { return value > rhs.value; }

  double value;
  int index;
};

double
estimateBandwidth(TheaArray<IndexedValue> const & sorted_pts)
{
  // See Silverman (1986), eq. 3.31, and
  // http://www1.american.edu/academic.depts/cas/econ/gaussres/utilitys/KERNEL/KONING/KNLIB.PS

  double sum_values = 0, sum_squares = 0;
  for (size_t i = 0; i < sorted_pts.size(); ++i)
  {
    double val = sorted_pts[i].value;
    sum_values  += val;
    sum_squares += (val * val);
  }

  double avg = sum_values / sorted_pts.size();
  double variance = (sum_squares / sorted_pts.size() - avg * avg);

  size_t q1 = sorted_pts.size() / 4;
  size_t q3 = (3 * sorted_pts.size()) / 4;
  double iqr = sorted_pts[q3].value - sorted_pts[q1].value;  // inter-quartile range

  double bandwidth = 0.9 * min(sqrt(variance), iqr / 1.34) / pow((double)sorted_pts.size(), 1.0 / 5.0);
  if (bandwidth < 1e-9) bandwidth = 0.005;

  return bandwidth;
}

double
doMeanShift(double start, TheaArray<IndexedValue> const & sorted_pts, double bandwidth, unsigned int num_iters = 100,
            double threshold = 0.001)
{
  double inv_bandwidth = 1.0 / bandwidth;

  for (unsigned int i = 0; i < num_iters; ++i)
  {
    // Speed up things by only looking at points within a window
    double lo = start - 1.2 * bandwidth;
    double hi = start + 1.2 * bandwidth;

    typedef TheaArray<IndexedValue>::const_iterator Iterator;
    Iterator begin = lower_bound(sorted_pts.begin(), sorted_pts.end(), IndexedValue(lo, -1));
    Iterator end   = upper_bound(sorted_pts.begin(), sorted_pts.end(), IndexedValue(hi, -1));

    double weighted_sum = 0;
    double sum_weights = 0, x, w;
    for (Iterator pi = begin; pi != end; ++pi)
    {
      x = (pi->value - start) * inv_bandwidth;
      w = dKernel(x * x);
      weighted_sum += w * (pi->value);
      sum_weights += w;
    }

    if (Math::fuzzyNe(sum_weights, 0.0))
    {
      double new_start = weighted_sum / sum_weights;
      double shift = fabs(new_start - start);
      start = new_start;
      if (shift < threshold) break;
    }
    else
      break;
  }

  return start;
}

void
doMeanShiftBlock(TheaArray<IndexedValue> const * descs, TheaArray<IndexedValue> * modes, size_t begin, size_t end,
                 double bandwidth, int thread_index)
{
  static unsigned int const NUM_ITERS = 100;
  static double       const THRESHOLD = 0.0001;

  size_t actual_end = min(end, descs->size());
  for (size_t i = begin; i < actual_end; ++i)
    (*modes)[i] = IndexedValue(doMeanShift((*descs)[i].value, *descs, bandwidth, NUM_ITERS, THRESHOLD), (*descs)[i].index);
}

typedef UnionFind<> LabelUnionFind;

int
defaultConcurrency()
{
  return (int)System::concurrency();
}

int
countSDFModes(TheaArray<Real> const & sdf_values)
{
  // Associate descriptors with sample indices
  TheaArray<IndexedValue> descs(sdf_values.size());
  for (size_t i = 0; i < sdf_values.size(); ++i)
    descs[i] = IndexedValue(sdf_values[i], (int)i);

  // Sort them, so we can quickly get all points in a range
  sort(descs.begin(), descs.end());

  double bandwidth = estimateBandwidth(descs);  // 0.005
  THEA_CONSOLE << "Estimated bandwidth = " << bandwidth;

  // Do mean shift on each point to find its nearest mode
  TheaArray<IndexedValue> modes(descs.size());

  size_t num_threads = (size_t)defaultConcurrency();  // numCores()
  THEA_CONSOLE << "Using " << num_threads << " threads to find initial modes via mean-shift";

  size_t samples_per_thread = (descs.size() <= num_threads
                                  ? 1 : (size_t)ceil(descs.size() / (double)num_threads));

  boost::thread_group threads;
  int thread_index = 0;
  for (size_t begin = 0; begin < descs.size(); begin += samples_per_thread, ++thread_index)
  {
    threads.add_thread(new boost::thread(doMeanShiftBlock, &descs, &modes, begin, begin + samples_per_thread, bandwidth,
                                         thread_index));
  }

  threads.join_all();

  // Sort the set of modes
  sort(modes.begin(), modes.end());

  // Create a union-find structure that starts with a set for every sample and finishes with a set for every cluster
  LabelUnionFind uf((long)modes.size());

  // Combine all means that are within a small threshold of each other
  static double const THRESHOLD_SCALE = 1.3;
  static long const MAX_CLUSTERS = 10;
  static long const REQUIRED_DIFF = 3;
  static long const QUICK_STOP = 5;
  static unsigned int MAX_ITERS = 100;
  double merge_threshold = 0.001 * bandwidth;

  long last_num_sets = uf.numSets();
  long last_diff = numeric_limits<long>::max() / 4;
  for (unsigned int iter = 1; iter <= MAX_ITERS; ++iter)
  {
    for (size_t i = 1; i < modes.size(); ++i)
      if (fabs(modes[i].value - modes[i - 1].value) < merge_threshold)
      {
        uf.merge((long)i, (long)i - 1);  // indices are into modes array
      }

    THEA_CONSOLE << uf.numSets() << " clusters identified after merge pass " << iter << " with threshold " << merge_threshold;

    long diff = max(last_num_sets - uf.numSets(), 1L);
    if (uf.numSets() <= MAX_CLUSTERS)
    {
      // If there are too few sets, or we just made a large jump to enter the allowed region, or there is a sharp drop in the
      // number of sets, stop
      if (uf.numSets() <= QUICK_STOP
       || last_num_sets - uf.numSets() >= REQUIRED_DIFF
       || diff > 2 * last_diff)
        break;
    }

    merge_threshold *= THRESHOLD_SCALE;
    last_num_sets = uf.numSets();
    last_diff = diff;
  }

  return (int)uf.numSets();
}

struct SampleCluster;  // forward declaration
typedef Graph<SampleCluster *, double> SampleClusterConnectivityGraph;

/** A cluster of samples, typically coresponding to a single segment of the model. */
struct SampleCluster
{
  int label;
  TheaArray<Vector3 const *> positions;
  TheaArray<Vector3 const *> normals;
  KDTreeN<Vector3 const *, 3> * kdtree;
  double concavity;
  SampleClusterConnectivityGraph::VertexIterator conn_vertex;
  bool changed;

  SampleCluster() : kdtree(new KDTreeN<Vector3 const *, 3>), changed(false) {}

  void updateKDTree()
  {
    if (changed)
    {
      kdtree->init(positions.begin(), positions.end());
      changed = false;
    }
  }
};

typedef TheaUnorderedMap<int, SampleCluster> SampleClusterMap;

double
computeConcavity(TheaArray<Vector3 const *> const & positions, TheaArray<Vector3 const *> const & normals,
                 KDTree const & kdtree)
{
  Real skin_width = 0.01 * kdtree.getBounds().getExtent().length();
  THEA_CONSOLE << "Skin width = " << skin_width;

  ConvexHull3::Options(ConvexHull3::Options::Approx(100, skin_width));
  ConvexHull3 hull;
  for (size_t i = 0; i < positions.size(); ++i)
    hull.addPoint(*positions[i]);

  Mesh hull_mesh;
  hull.computeApprox(hull_mesh);

  KDTree hull_kdtree;
  hull_kdtree.add(hull_mesh);
  hull_kdtree.init();

  Real extent = kdtree.getBounds().getExtent().length();
  THEA_CONSOLE << "Extent = " << extent;
  if (Math::fuzzyEq(extent, (Real)0))
    return 0;

  double sum_distances = 0;
  long num_points = 0;
  for (size_t i = 0; i < positions.size(); ++i)
  {
    Ray3 hull_ray(*positions[i], *normals[i]);
    Real hull_isec_time = hull_kdtree.rayIntersectionTime<RayIntersectionTester>(hull_ray);
    // cout << "Intersection time for ray " << i << " = " << isec_time << endl;

    Ray3 model_ray(hull_ray.getOrigin() + 0.001 * hull_ray.getDirection(), hull_ray.getDirection());
    Real model_isec_time = kdtree.rayIntersectionTime<RayIntersectionTester>(model_ray);
    if (hull_isec_time < model_isec_time)
    {
      sum_distances += min(hull_isec_time, extent);
      num_points++;
    }
  }

  THEA_CONSOLE << "Sum of distances for " << num_points << " points = " << sum_distances;

  return num_points <= 0 ? 0 : sum_distances / (num_points * extent);
}

int
combineClustersByConvexity(TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals, KDTree const & kdtree,
                           TheaArray<int> & labels, double concavity_threshold, double score_threshold)
{
  if (concavity_threshold < 0) concavity_threshold  =  0.015;
  if (score_threshold     < 0) score_threshold      =  0.9;

  SampleClusterMap sample_clusters;
  {
    for (size_t i = 0; i < positions.size(); ++i)
    {
      sample_clusters[labels[i]].positions.push_back(&positions[i]);
      sample_clusters[labels[i]].normals.push_back(&normals[i]);
    }

    for (SampleClusterMap::iterator ci = sample_clusters.begin(); ci != sample_clusters.end(); ++ci)
    {
      ci->second.label = ci->first;
      ci->second.changed = true;
    }
  }

  SampleClusterConnectivityGraph cluster_conn_graph;
  {
    double       const INTERSECTION_THRESHOLD  =  0.01 * kdtree.getBounds().getExtent().length();
    size_t const NNBRS_THRESHOLD         =  10;

    for (SampleClusterMap::iterator ci = sample_clusters.begin(); ci != sample_clusters.end(); ++ci)
    {
      ci->second.updateKDTree();
      ci->second.conn_vertex = cluster_conn_graph.addVertex(&ci->second);

      for (SampleClusterMap::iterator cj = sample_clusters.begin(); cj != ci; ++cj)
      {
        size_t nnbrs = 0;
        for (size_t j = 0; j < cj->second.positions.size(); ++j)
        {
          long nn_index = ci->second.kdtree->closestElement<MetricL2>(*cj->second.positions[j], INTERSECTION_THRESHOLD);
          if (nn_index >= 0)
          {
            if (++nnbrs >= min(min(ci->second.positions.size(), cj->second.positions.size()), NNBRS_THRESHOLD))
            {
              THEA_CONSOLE << "Clusters " << cj->second.label << " and " << ci->second.label << " are connected ";
              cluster_conn_graph.addEdge(cj->second.conn_vertex, ci->second.conn_vertex, 0);
              break;
            }
          }
        }
      }
    }
  }

  // Compute cluster concavities
  for (SampleClusterConnectivityGraph::VertexIterator vi = cluster_conn_graph.verticesBegin();
       vi != cluster_conn_graph.verticesEnd(); ++vi)
  {
    vi->attr()->concavity = computeConcavity(vi->attr()->positions, vi->attr()->normals, kdtree);
  }

  // Compute concavities of connected pairs of clusters
  for (SampleClusterConnectivityGraph::EdgeIterator ei = cluster_conn_graph.edgesBegin(); ei != cluster_conn_graph.edgesEnd();
       ++ei)
  {
    TheaArray<Vector3 const *> combined_positions(ei->getOrigin()->attr()->positions);
    combined_positions.insert(combined_positions.end(), ei->getEnd()->attr()->positions.begin(),
                                                        ei->getEnd()->attr()->positions.end());

    TheaArray<Vector3 const *> combined_normals(ei->getOrigin()->attr()->normals);
    combined_normals.insert(combined_normals.end(), ei->getEnd()->attr()->normals.begin(),
                                                    ei->getEnd()->attr()->normals.end());

    ei->setAttr(computeConcavity(combined_positions, combined_normals, kdtree));
  }

  printGraph(cluster_conn_graph);

  while (cluster_conn_graph.numEdges() > 1)
  {
    SampleClusterConnectivityGraph::EdgeIterator max_edge = cluster_conn_graph.edgesEnd();
    double max_score = -1e+30;
    for (SampleClusterConnectivityGraph::EdgeIterator ei = cluster_conn_graph.edgesBegin(); ei != cluster_conn_graph.edgesEnd();
         ++ei)
    {
      double edge_concavity = ei->attr();
      double vertex0_concavity = ei->getOrigin()->attr()->concavity;
      double vertex1_concavity = ei->getEnd()->attr()->concavity;
      double vertex_concavity = min(vertex0_concavity, vertex1_concavity);
      double score = vertex_concavity / max(edge_concavity, 1e-10);

      THEA_CONSOLE << "Score of merging clusters " << ei->getOrigin()->attr()->label << " and " << ei->getEnd()->attr()->label
                   << " = " << score;

      if (score > max_score && (edge_concavity < concavity_threshold || score > score_threshold))
      {
        max_score = score;
        max_edge = ei;
      }
    }

    if (max_edge == cluster_conn_graph.edgesEnd())
      break;

    THEA_CONSOLE << "Merging clusters " << max_edge->getOrigin()->attr()->label << " and " << max_edge->getEnd()->attr()->label;

    // Collapse the edge. The convexity of the vertex formed by merging the endpoints is the convexity stored at the edge.
    SampleClusterConnectivityGraph::Vertex * ov = max_edge->getOrigin();
    ov->attr()->concavity = max_edge->attr();

    // Merge the sample lists of the two clusters
    TheaArray<Vector3 const *> & opositions = ov->attr()->positions;
    TheaArray<Vector3 const *> & epositions = max_edge->getEnd()->attr()->positions;
    opositions.insert(opositions.end(), epositions.begin(), epositions.end());

    TheaArray<Vector3 const *> & onormals = ov->attr()->normals;
    TheaArray<Vector3 const *> & enormals = max_edge->getEnd()->attr()->normals;
    onormals.insert(onormals.end(), enormals.begin(), enormals.end());

    // Mark the KD-tree for an update
    ov->attr()->changed = true;

    cluster_conn_graph.collapseEdge(max_edge);
    cluster_conn_graph.mergeTwinEdges(SampleClusterConnectivityGraph::VertexIterator(ov));

    printGraph(cluster_conn_graph);

    // Recompute the convexity values for all edges incident at this vertex
    for (SampleClusterConnectivityGraph::Vertex::EdgeIterator ei = ov->incomingEdgesBegin(); ei != ov->incomingEdgesEnd(); ++ei)
    {
      TheaArray<Vector3 const *> combined_positions((*ei)->getOrigin()->attr()->positions);
      combined_positions.insert(combined_positions.end(), opositions.begin(), opositions.end());

      TheaArray<Vector3 const *> combined_normals((*ei)->getOrigin()->attr()->normals);
      combined_normals.insert(combined_normals.end(), onormals.begin(), onormals.end());

      double combined_concavity = computeConcavity(combined_positions, combined_normals, kdtree);
      (*ei)->setAttr(combined_concavity);
    }

    for (SampleClusterConnectivityGraph::Vertex::EdgeIterator ei = ov->outgoingEdgesBegin(); ei != ov->outgoingEdgesEnd(); ++ei)
    {
      TheaArray<Vector3 const *> combined_positions((*ei)->getEnd()->attr()->positions);
      combined_positions.insert(combined_positions.end(), opositions.begin(), opositions.end());

      TheaArray<Vector3 const *> combined_normals((*ei)->getEnd()->attr()->normals);
      combined_normals.insert(combined_normals.end(), onormals.begin(), onormals.end());

      double combined_concavity = computeConcavity(combined_positions, combined_normals, kdtree);
      (*ei)->setAttr(combined_concavity);
    }
  }

  // Now relabel the points
  int curr_label = 0;
  for (SampleClusterConnectivityGraph::VertexIterator vi = cluster_conn_graph.verticesBegin();
       vi != cluster_conn_graph.verticesEnd(); ++vi)
  {
    int label = -1;
    if (vi->attr()->label >= 0)
      label = curr_label++;

    for (size_t i = 0; i < vi->attr()->positions.size(); ++i)
    {
      size_t index = (size_t)(vi->attr()->positions[i] - &positions[0]);
      labels[index] = label;
    }
  }

  return (int)cluster_conn_graph.numVertices();
}
