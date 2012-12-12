#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/ShapeDiameter.hpp"

#ifndef THEA_OSX
#  define THEA_ENABLE_CLUTO
#  include "../../Algorithms/Clustering.hpp"
#endif

#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Math.hpp"
#include "../../VectorN.hpp"
#include <CGAL/Union_find.h>
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
Color3 const & getPaletteColor(int i);

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
    total_area += kdtree.getElements()[(array_size_t)i].getArea();

  double density = approx_num_samples / total_area;

  TheaArray<Vector3> positions, normals;
  for (long i = 0; i < kdtree.numElements(); ++i)
  {
    KDTree::Element const & elem = kdtree.getElements()[(array_size_t)i];
    double num_samples = density * elem.getArea();
    double rem = num_samples;
    for (long j = 1; j < num_samples; ++j)
    {
      positions.push_back(elem.randomPoint());
      normals.push_back(elem.getNormal());
      rem -= 1;
    }

    if (Math::rand01() < rem)
    {
      positions.push_back(elem.randomPoint());
      normals.push_back(elem.getNormal());
    }
  }

  THEA_CONSOLE << "Computed " << positions.size() << " sample points on the mesh";

  // Compute SDF values
  TheaArray<Real> sdf_values;
  MeshFeatures::ShapeDiameter<Mesh> sdf(&kdtree);
  sdf.compute(positions, normals, sdf_values);

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
  for (array_size_t i = 0; i < positions.size(); ++i)
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

  // Write out the labels
  ofstream out(outpath.c_str());
  out << positions.size() << endl;
  out << num_clusters << endl;
  for (array_size_t i = 0; i < positions.size(); ++i)
  {
#ifdef DEBUG_PTS
    Color3 pseudo_n = getPaletteColor(labels[i]);
    out << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << ' '
        << pseudo_n.r      << ' ' << pseudo_n.g      << ' ' << pseudo_n.b << endl;
#else
    out << positions[i][0] << ' ' << positions[i][1] << ' ' << positions[i][2] << ' ' << sdf_values[i] << ' ' << labels[i]
        << endl;
#endif
  }

#endif

  return 0;
}

static int const NUM_PALETTE_COLORS = 24;
Color3 COLOR_PALETTE[NUM_PALETTE_COLORS] = {
  Color3::fromARGB(0x298edb),
  Color3::fromARGB(0x982411),
  Color3::fromARGB(0x6d4e25),
  Color3::fromARGB(0x1b5043),
  Color3::fromARGB(0x6e7662),
  Color3::fromARGB(0xa08b00),
  Color3::fromARGB(0x58427b),
  Color3::fromARGB(0x1d2f5b),
  Color3::fromARGB(0xac5e34),
  Color3::fromARGB(0x804055),
  Color3::fromARGB(0x6d7a00),
  Color3::fromARGB(0x572e2c),

  // Invert each color above
  Color3::fromARGB(~0x298edb),
  Color3::fromARGB(~0x982411),
  Color3::fromARGB(~0x6d4e25),
  Color3::fromARGB(~0x1b5043),
  Color3::fromARGB(~0x6e7662),
  Color3::fromARGB(~0xa08b00),
  Color3::fromARGB(~0x58427b),
  Color3::fromARGB(~0x1d2f5b),
  Color3::fromARGB(~0xac5e34),
  Color3::fromARGB(~0x804055),
  Color3::fromARGB(~0x6d7a00),
  Color3::fromARGB(~0x572e2c)
};

int
numPaletteColors()
{
  return NUM_PALETTE_COLORS;
}

Color3 const &
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
  for (array_size_t i = 0; i < sorted_pts.size(); ++i)
  {
    double val = sorted_pts[i].value;
    sum_values  += val;
    sum_squares += (val * val);
  }

  double avg = sum_values / sorted_pts.size();
  double variance = (sum_squares / sorted_pts.size() - avg * avg);

  array_size_t q1 = sorted_pts.size() / 4;
  array_size_t q3 = (3 * sorted_pts.size()) / 4;
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
doMeanShiftBlock(TheaArray<IndexedValue> const * descs, TheaArray<IndexedValue> * modes, array_size_t begin, array_size_t end,
                 double bandwidth, int thread_index)
{
  static unsigned int const NUM_ITERS = 100;
  static double       const THRESHOLD = 0.0001;

  array_size_t actual_end = min(end, descs->size());
  for (array_size_t i = begin; i < actual_end; ++i)
    (*modes)[i] = IndexedValue(doMeanShift((*descs)[i].value, *descs, bandwidth, NUM_ITERS, THRESHOLD), (*descs)[i].index);
}

/** Union-find data structure, wrapper for CGAL::Union_find. */
template < typename T, typename A = CGAL_ALLOCATOR(T) >
class UnionFind : public CGAL::Union_find<T, A>
{
  private:
    typedef CGAL::Union_find<T, A> BaseType;

  public:
    /** Hash functor for mutable handle. */
    struct HandleHash
    {
      size_t operator()(typename BaseType::handle ufh) const { return reinterpret_cast<size_t>(&(*ufh)); }
    };

    /** Hash functor for immutable handle. */
    struct ConstHandleHash
    {
      size_t operator()(typename BaseType::const_handle ufh) const { return reinterpret_cast<size_t>(&(*ufh)); }
    };

}; // class UnionFind

typedef UnionFind<int> LabelUnionFind;

int
defaultConcurrency()
{
  int cc = (int)boost::thread::hardware_concurrency();
  if (cc <= 0)
    return G3D::System::numCores();
  else
    return cc;
}

int
countSDFModes(TheaArray<Real> const & sdf_values)
{
  // Associate descriptors with sample indices
  TheaArray<IndexedValue> descs(sdf_values.size());
  for (array_size_t i = 0; i < sdf_values.size(); ++i)
    descs[i] = IndexedValue(sdf_values[i], (int)i);

  // Sort them, so we can quickly get all points in a range
  sort(descs.begin(), descs.end());

  double bandwidth = estimateBandwidth(descs);  // 0.005
  THEA_CONSOLE << "Estimated bandwidth = " << bandwidth;

  // Do mean shift on each point to find its nearest mode
  TheaArray<IndexedValue> modes(descs.size());

  array_size_t num_threads = (array_size_t)defaultConcurrency();  // numCores()
  THEA_CONSOLE << "Using " << num_threads << " threads to find initial modes via mean-shift";

  array_size_t samples_per_thread = (descs.size() <= num_threads
                                  ? 1 : (array_size_t)ceil(descs.size() / (double)num_threads));

  boost::thread_group threads;
  int thread_index = 0;
  for (array_size_t begin = 0; begin < descs.size(); begin += samples_per_thread, ++thread_index)
  {
    threads.add_thread(new boost::thread(doMeanShiftBlock, &descs, &modes, begin, begin + samples_per_thread, bandwidth,
                                         thread_index));
  }

  threads.join_all();

  // Sort the set of modes
  sort(modes.begin(), modes.end());

  // Create a union-find structure that starts with a set for every sample and finishes with a set for every cluster
  LabelUnionFind uf;

  TheaArray<LabelUnionFind::handle> handles(modes.size());
  for (array_size_t i = 0; i < modes.size(); ++i)
    handles[i] = uf.push_back(i);  // indices into the modes array

  // Combine all means that are within a small threshold of each other
  static double const THRESHOLD_SCALE = 1.3;
  static size_t const MAX_CLUSTERS = 10;
  static size_t const REQUIRED_DIFF = 3;
  static size_t const QUICK_STOP = 5;
  static unsigned int MAX_ITERS = 100;
  double merge_threshold = 0.001 * bandwidth;

  size_t last_num_sets = uf.number_of_sets();
  size_t last_diff = numeric_limits<size_t>::max() / 4;
  for (unsigned int iter = 1; iter <= MAX_ITERS; ++iter)
  {
    for (array_size_t i = 1; i < modes.size(); ++i)
      if (fabs(modes[i].value - modes[i - 1].value) < merge_threshold)
      {
        uf.unify_sets(handles[i], handles[i - 1]);  // indices are into modes array
      }

    THEA_CONSOLE << uf.number_of_sets() << " clusters identified after merge pass " << iter << " with threshold "
                 << merge_threshold;

    size_t diff = max(last_num_sets - uf.number_of_sets(), (size_t)1);
    if (uf.number_of_sets() <= MAX_CLUSTERS)
    {
      // If there are too few sets, or we just made a large jump to enter the allowed region, or there is a sharp drop in the
      // number of sets, stop
      if (uf.number_of_sets() <= QUICK_STOP
       || last_num_sets - uf.number_of_sets() >= REQUIRED_DIFF
       || diff > 2 * last_diff)
        break;
    }

    merge_threshold *= THRESHOLD_SCALE;
    last_num_sets = uf.number_of_sets();
    last_diff = diff;
  }

  return (int)uf.number_of_sets();
}
