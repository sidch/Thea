#include "../../Common.hpp"
#include "../../FilePath.hpp"
#include "../../Algorithms/Filter.hpp"
#include "../../Algorithms/Iterators.hpp"
#include "../../Algorithms/BvhN.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../AffineTransform3.hpp"
#include "../../BoundedArrayN.hpp"
#include "../../BoundedSortedArrayN.hpp"
#include "../../MatVec.hpp"
#include "../../UnorderedMap.hpp"
#include "../../UnionFind.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>

// #define CONNECTED_COMPONENTS

using namespace std;
using namespace Thea;
using namespace Algorithms;

typedef UnorderedMap<int32, string> IndexLabelMap;

static size_t const MAX_NBRS = 8;
int MAX_ROUNDS = 25;
int MAX_SMOOTHING_ROUNDS = 3;
int LAST_SMOOTHING_ROUNDS = 1;
bool USE_LABELS = false;
bool USE_NORMALS = true;
IndexLabelMap LABELS;

struct Sample
{
  Sample() : active(true) {}
  Sample(Vector3 const & p_, Vector3 const & n_ = Vector3::Zero(), int32 label_ = -1, bool active_ = true)
  : p(p_), n(n_), label(label_), active(active_) {}

  Vector3 p;
  Vector3 n;
  int32 label;
  bool active;
};

class Offset
{
  private:
    Vector3 dir;
    bool valid;

  public:
    Offset() : dir(0, 0, 0), valid(false) {}
    Offset(Vector3 const & d_) : dir(d_), valid(true) {}
    Vector3 const & d() const { return dir; }
    bool isValid() const { return valid; }
    void set(Vector3 const & d_) { dir = d_; valid = true; }
    void unset() { valid = false; dir = Vector3::Zero(); }
};

namespace Thea {
namespace Algorithms {

template <>
class IsPointN<Sample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<Sample, 3>::getPosition(Sample const & sample)
{
  return sample.p;
}

} // namespace Algorithms
} // namespace Thea

typedef Array<Sample> SampleArray;
typedef PtrIterator<SampleArray::const_iterator> SamplePtrIterator;
typedef BvhN<Sample const *, 3> Bvh;
typedef BoundedArrayN<MAX_NBRS, size_t> NeighborSet;
typedef Array<NeighborSet> NeighborSets;
typedef Array<Offset> OffsetArray;

struct NNFilter : public Filter<Sample const *>
{
  NNFilter(Vector3 const & n_, int32 label_) : n(n_), label(label_) {}

  bool allows(Sample const * const & sample) const
  { return (!USE_LABELS || sample->label == label) && (!USE_NORMALS || n.dot(sample->n) >= 0); }

  Vector3 n;
  int32 label;
};

struct NNFilterLabelOnly : public Filter<Sample const *>
{
  NNFilterLabelOnly(int32 label_) : label(label_) {}
  bool allows(Sample const * const & sample) const { return (!USE_LABELS || sample->label == label); }

  int32 label;
};

inline Real
kernelEpanechnikovSqDist(Real squared_dist, Real squared_bandwidth)
{
  static Real const SCALE = (Real)(15.0 / (8.0 * Math::pi()));  // normalizing constant for 3D Epanechnikov kernel

  Real nrm_sqdist = squared_dist / squared_bandwidth;
  return nrm_sqdist < 1 ? SCALE * (1 - nrm_sqdist) : 0;
}

string
readQuotedString(istream & in)
{
  char c;
  while (in.get(c))
  {
    if (c == '"')
    {
      ostringstream oss;
      while (in.get(c))
      {
        if (c == '"') break;
        oss << c;
      }

      return oss.str();
    }
  }

  return "";
}

// Guaranteed to return a value between 1 and 2^31 - 1
int32
labelHash(string const & label)
{
  std::hash<string> hasher;
  return max((int32)(hasher(label) & 0x7FFFFFFF), (int32)1);
}

bool
loadSamples(string const & path, SampleArray & samples, bool read_labels)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open samples file " << path;
    return false;
  }

  samples.clear();
  string line;
  Vector3 p, n = Vector3::Zero();
  while (getline(in, line))
  {
    istringstream line_in(line);
    if (USE_NORMALS)
    { if (!(line_in >> p[0] >> p[1] >> p[2] >> n[0] >> n[1] >> n[2])) { continue; } }
    else
    { if (!(line_in >> p[0] >> p[1] >> p[2])) { continue; } }

    int32 label_hash = 0;
    if (read_labels)
    {
      string label = readQuotedString(line_in);
      if (!trimWhitespace(label).empty())
      {
        label_hash = labelHash(label);
        if (LABELS.find(label_hash) == LABELS.end())
          LABELS[label_hash] = label;
      }
    }

    samples.push_back(Sample(p, n, label_hash));
  }

  if (samples.empty())
  {
    THEA_ERROR << "Could not read any samples from file '" << path << '\'';
    return false;
  }

  THEA_CONSOLE << "Read " << samples.size() << " samples from file '" << path << '\'';

  return true;
}

AxisAlignedBox3
computeBounds(SampleArray const & samples, int32 selected_label = -1)
{
  AxisAlignedBox3 aabb;
  for (size_t i = 0; i < samples.size(); ++i)
  {
    if (selected_label < 0 || samples[i].label == selected_label)
      aabb.merge(samples[i].p);
  }

  return aabb;
}

void
findNeighbors(SampleArray const & samples, Bvh const & bvh, NeighborSets & nbrs)
{
  nbrs.clear();
  nbrs.resize(samples.size());

  BoundedSortedArrayN<3 * MAX_NBRS, Bvh::NeighborPair> init_nbrs;
  size_t last_count = 0;
  for (size_t i = 0; i < samples.size(); ++i)
  {
    NNFilter filter(samples[i].n, samples[i].label);
    const_cast<Bvh &>(bvh).pushFilter(&filter);
      bvh.kClosestPairs<MetricL2>(samples[i].p, init_nbrs);
    const_cast<Bvh &>(bvh).popFilter();

    for (size_t j = 0; j < init_nbrs.size() && nbrs[i].size() < MAX_NBRS; ++j)
    {
      intx tgt_index = init_nbrs[j].getTargetIndex();
      if (tgt_index < 0)
        continue;

      nbrs[i].push_back((size_t)tgt_index);
    }

    if (i > 0 && i % 1000 == 0)
    {
      THEA_CONSOLE << "Found neighbors of " << i << " samples";
      last_count = i;
    }
  }

  if (last_count + 1 < samples.size())
    THEA_CONSOLE << "Found neighbors of " << samples.size() << " samples";

  if (!samples.empty())
  {
    double avg_nnbrs = 0;
    for (size_t i = 0; i < nbrs.size(); ++i)
      avg_nnbrs += nbrs[i].size();

    avg_nnbrs /= nbrs.size();

    THEA_CONSOLE << "Samples have " << avg_nnbrs << " neighbors on average";
  }
}

void
updateOffsets(SampleArray const & src_samples, Bvh const & tgt_bvh, OffsetArray & src_offsets, int32 selected_label = -1)
{
  for (size_t i = 0; i < src_samples.size(); ++i)
  {
    if (selected_label >= 0 && src_samples[i].label != selected_label)
      continue;

    if (!src_samples[i].active)
    {
      src_offsets[i].unset();
      continue;
    }

    Vector3 offset_p = src_samples[i].p + src_offsets[i].d();

    NNFilter filter(src_samples[i].n, src_samples[i].label);
    const_cast<Bvh &>(tgt_bvh).pushFilter(&filter);
      intx nn_index = tgt_bvh.closestElement<MetricL2>(offset_p);
    const_cast<Bvh &>(tgt_bvh).popFilter();

    if (nn_index >= 0)
    {
      Sample const * nn_sample = tgt_bvh.getElements()[(size_t)nn_index];
      src_offsets[i].set(nn_sample->p - src_samples[i].p);  // no transform here
    }
    else
      src_offsets[i].unset();
  }

  if (selected_label < 0)
    cout << " update(" << src_offsets.size() << ')' << flush;
}

void
smoothOffsets(SampleArray const & samples, NeighborSets const & nbrs, OffsetArray & offsets, int num_rounds = 1)
{
  if (num_rounds < 1)
    return;

  for (int round = 0; round < num_rounds; ++round)
  {
    OffsetArray smoothed_offsets(samples.size());

    for (size_t i = 0; i < samples.size(); ++i)
    {
      if (nbrs[i].empty())
        continue;

      Real sq_bandwidth = 4 * (samples[nbrs[i].back()].p - samples[i].p).squaredNorm();

// #define SPHERE_PRIOR
#ifdef SPHERE_PRIOR

      // Even if there is no valid offset, we want the point's original position to be favored
      Real sum_weights = kernelEpanechnikovSqDist(0, sq_bandwidth);

      Vector3 sum_dirs(0, 0, 0);
      Real sum_lengths = 0;
      Real length = 0;
      if (offsets[i].isValid() && (length = offsets[i].d().norm()) > 0)
      {
        sum_dirs = (sum_weights / length) * offsets[i].d();
        sum_lengths = (sum_weights / length);
      }

      for (size_t j = 0; j < nbrs[i].size(); ++j)
      {
        Offset const & offset = offsets[nbrs[i][j]];
        if (!offset.isValid())
          continue;

        Real sqdist = (samples[nbrs[i][j]].p - samples[i].p).squaredNorm();
        Real weight = kernelEpanechnikovSqDist(sqdist, sq_bandwidth);

        length = offset.d().norm();
        if (length > 0)
        {
          sum_dirs += ((weight / length) * offset.d());

          if (offsets[i].isValid())
          {
            if (offset.d().dot(offsets[i].d()) < 0)
              sum_lengths -= (weight * length);
            else
              sum_lengths += (weight * length);
          }
          else
            sum_lengths += fabs(weight * length);
        }

        sum_weights += weight;
      }

      smoothed_offsets[i].set(fabs(sum_lengths / sum_weights) * sum_dirs.stableNormalized());

#else // plane prior

      // Even if there is no valid offset, we want the point's original position to be favored
      Real sum_weights = 0;
      Vector3 sum_offsets = Vector3::Zero();
      if (offsets[i].isValid())
      {
        sum_weights = kernelEpanechnikovSqDist(0, sq_bandwidth);
        sum_offsets = sum_weights * offsets[i].d();
      }

      for (size_t j = 0; j < nbrs[i].size(); ++j)
      {
        Offset const & offset = offsets[nbrs[i][j]];
        if (!offset.isValid())
          continue;

        Real sqdist = (samples[nbrs[i][j]].p - samples[i].p).squaredNorm();
        Real weight = kernelEpanechnikovSqDist(sqdist, sq_bandwidth);

        sum_offsets += (weight * offset.d());
        sum_weights += weight;
      }

      if (sum_weights > 0)
        smoothed_offsets[i].set(sum_offsets / sum_weights);
      else
        smoothed_offsets[i].set(Vector3::Zero());

#endif
    }

    offsets = smoothed_offsets;
  }

  cout << " smooth(" << offsets.size() << ", " << num_rounds << ')' << flush;
}

#ifdef CONNECTED_COMPONENTS

Vector3
computeCentroid(SampleArray const & samples, int32 selected_label = -1, bool only_active = false)
{
  Vector3 sum_c = Vector3::Zero();
  double sum_weights = 0;
  for (size_t i = 0; i < samples.size(); ++i)
  {
    if ((selected_label < 0 || samples[i].label == selected_label) && (!only_active || samples[i].active))
      continue;

    sum_c += samples[i].p;
    sum_weights += 1.0;
  }

  return (sum_weights > 0 ? sum_c / sum_weights : Vector3::Zero());
}

Vector3
computeWeightedCentroid(SampleArray const & samples, NeighborSets const & nbrs, int32 selected_label = -1,
                        bool only_active = false)
{
  Array<size_t> selected_samples;
  for (size_t i = 0; i < samples.size(); ++i)
    if ((selected_label < 0 || samples[i].label == selected_label) && (!only_active || samples[i].active))
      selected_samples.push_back(i);

  if (selected_samples.empty())
    return Vector3::Zero();

  UnionFind<size_t> uf(selected_samples.begin(), selected_samples.end());
  for (size_t i = 0; i < selected_samples.size(); ++i)
  {
    size_t index = selected_samples[i];
    intx src_id = uf.getObjectId(index);
    for (size_t j = 0; j < nbrs[index].size(); ++j)
      uf.merge(src_id, uf.getObjectId(nbrs[index][j]));
  }

  Vector3 sum_c = Vector3::Zero();
  double sum_weights = 0;
  for (size_t i = 0; i < selected_samples.size(); ++i)
  {
    size_t index = selected_samples[i];
    double weight = uf.sizeOfSet(uf.getObjectId(index)) / (double)selected_samples.size();
    sum_c += weight * samples[index].p;
    sum_weights += weight;
  }

  return (sum_weights > 0 ? sum_c / sum_weights : Vector3::Zero());
}

#endif // CONNECTED_COMPONENTS

// Returns transform of samples1 to align them to samples2
AffineTransform3
initTransform(int32 selected_label,
              SampleArray const & samples1, NeighborSets const & nbrs1, Bvh & bvh1,
              SampleArray const & samples2, NeighborSets const & nbrs2, Bvh & bvh2)
{
  AffineTransform3 tr = AffineTransform3::identity();

#ifdef CONNECTED_COMPONENTS
  if (selected_label < 0)
#endif
  {
    AxisAlignedBox3 selected_bounds1 = computeBounds(samples1, selected_label);
    AxisAlignedBox3 selected_bounds2 = computeBounds(samples2, selected_label);

    tr = AffineTransform3::translation(selected_bounds2.getCenter());
    if (selected_bounds1.getExtent().squaredNorm() > 1e-20 && selected_bounds2.getExtent().squaredNorm() > 1e-20)
      tr = tr * AffineTransform3::scaling(selected_bounds2.getExtent().cwiseQuotient(selected_bounds1.getExtent()));
    tr = tr * AffineTransform3::translation(-selected_bounds1.getCenter());
  }
#ifdef CONNECTED_COMPONENTS
  else
  {
    Vector3 c1 = computeCentroid(samples1, selected_label, true);
    Vector3 c2 = computeCentroid(samples2, selected_label, true);

    tr = AffineTransform3::translation(c2 - c1);
  }
#endif

  return tr;
}

#ifdef CONNECTED_COMPONENTS

void
selectSamples(SampleArray const & all_samples, int selected_label, Array<size_t> & selected_samples)
{
  selected_samples.clear();

  for (size_t i = 0; i < all_samples.size(); ++i)
    if (selected_label < 0 || all_samples[i].label == selected_label)
      selected_samples.push_back(i);
}

void
findConnectedComponents(Array<size_t> const & selected_samples, NeighborSets const & nbrs,
                        UnionFind<size_t> & uf)
{
  for (size_t i = 0; i < selected_samples.size(); ++i)
  {
    size_t index = selected_samples[i];
    intx src_id = uf.getObjectId(index);
    for (size_t j = 0; j < nbrs[index].size(); ++j)
      uf.merge(src_id, uf.getObjectId(nbrs[index][j]));
  }
}

void
computeComponentProperties(SampleArray const & all_samples, Array<size_t> const & selected_samples,
                           UnionFind<size_t> const & uf, Array<intx> & cc_reps, Array<intx> & cc_counts,
                           Array<AxisAlignedBox3> & cc_bounds)
{
  size_t ncc = (size_t)uf.numSets();
  cc_reps.resize(ncc); fill(cc_reps.begin(), cc_reps.end(), -1);
  cc_counts.resize(ncc); fill(cc_counts.begin(), cc_counts.end(), 0);
  cc_bounds.resize(ncc); fill(cc_bounds.begin(), cc_bounds.end(), AxisAlignedBox3());

  typedef UnorderedMap<intx, size_t> RepSetMap;  // maps from ID of set's representative sample to set ID
  RepSetMap set_ids;

  for (size_t i = 0; i < selected_samples.size(); ++i)
  {
    intx rep_id = uf.find(uf.getObjectId(selected_samples[i]));
    RepSetMap::const_iterator existing = set_ids.find(rep_id);
    size_t set_id;
    if (existing == set_ids.end())
    {
      set_id = (size_t)set_ids.size();
      set_ids[rep_id] = set_id;
      cc_reps[set_id] = rep_id;
    }
    else
      set_id = existing->second;

    cc_counts[set_id]++;
    cc_bounds[set_id].merge(all_samples[selected_samples[i]].p);
  }
}

template <typename T>
void
sortIndices(Array<T> const & values, Array<size_t> & sorted_indices, bool descending = false)
{
  sorted_indices.resize(values.size());
  for (size_t i = 0; i < sorted_indices.size(); ++i)
    sorted_indices[i] = i;

  for (size_t i = 0; i < sorted_indices.size(); ++i)
    for (size_t j = i + 1; j < sorted_indices.size(); ++j)
      if ((values[sorted_indices[i]] < values[sorted_indices[j]]) == descending)
        swap(sorted_indices[i], sorted_indices[j]);
}

void
deactivateSmallComponents(Array<size_t> const & selected_samples, UnionFind<size_t> const & uf,
                          intx num_active_sets, Array<intx> const & cc_reps, Array<size_t> const & cc_sorted,
                          SampleArray & all_samples)
{
  for (size_t i = 0; i < selected_samples.size(); ++i)
  {
    intx rep = uf.find(uf.getObjectId(selected_samples[i]));
    for (size_t j = (size_t)num_active_sets; j < cc_sorted.size(); ++j)
      if (rep == cc_reps[cc_sorted[j]])
      {
        all_samples[selected_samples[i]].active = false;
        break;
      }
  }
}

#endif // CONNECTED_COMPONENTS

void
initOffsetsForLabel(int32 selected_label,
                    SampleArray & samples1, NeighborSets const & nbrs1, Bvh & bvh1,
                    SampleArray & samples2, NeighborSets const & nbrs2, Bvh & bvh2,
                    OffsetArray & offsets1, OffsetArray & offsets2)
{
  string label_name = "all";
  if (selected_label >= 0)
  {
    IndexLabelMap::const_iterator existing = LABELS.find(selected_label);
    label_name = (existing != LABELS.end() ? existing->second : "anonymous");
  }

#ifdef CONNECTED_COMPONENTS
  if (selected_label >= 0)
  {
    //=========================================================================================================================
    // Select samples to be processed
    //=========================================================================================================================

    Array<size_t> selected_samples1, selected_samples2;
    selectSamples(samples1, selected_label, selected_samples1);
    selectSamples(samples2, selected_label, selected_samples2);

    if (selected_samples1.empty() || selected_samples2.empty())
    {
      for (size_t i = 0; i < samples1.size(); ++i)
        if (selected_label < 0 || samples1[i].label == selected_label)
          offsets1[i].unset();

      for (size_t i = 0; i < samples2.size(); ++i)
        if (selected_label < 0 || samples2[i].label == selected_label)
          offsets2[i].unset();

      return;
    }

    //=========================================================================================================================
    // Find connected components
    //=========================================================================================================================

    UnionFind<size_t> uf1(selected_samples1.begin(), selected_samples1.end());
    findConnectedComponents(selected_samples1, nbrs1, uf1);

    UnionFind<size_t> uf2(selected_samples2.begin(), selected_samples2.end());
    findConnectedComponents(selected_samples2, nbrs2, uf2);

    THEA_CONSOLE << "Label '" << label_name << "' in shape 1 has " << uf1.numSets() << " connected component(s)";
    THEA_CONSOLE << "Label '" << label_name << "' in shape 2 has " << uf2.numSets() << " connected component(s)";

    //=========================================================================================================================
    // Map connected components to each other
    //=========================================================================================================================

    Array<intx> cc_reps1, cc_reps2;
    Array<intx> cc_counts1, cc_counts2;
    Array<AxisAlignedBox3> cc_bounds1, cc_bounds2;
    computeComponentProperties(samples1, selected_samples1, uf1, cc_reps1, cc_counts1, cc_bounds1);
    computeComponentProperties(samples2, selected_samples2, uf2, cc_reps2, cc_counts2, cc_bounds2);

    Array<size_t> cc_sorted1, cc_sorted2;
    sortIndices(cc_counts1, cc_sorted1, true);
    sortIndices(cc_counts2, cc_sorted2, true);

    if (cc_sorted1.size() > 2)
      deactivateSmallComponents(selected_samples1, uf1, 2, cc_reps1, cc_sorted1, samples1);

    if (cc_sorted2.size() > 2)
      deactivateSmallComponents(selected_samples2, uf2, 2, cc_reps2, cc_sorted2, samples2);
  }
#endif

  AffineTransform3 tr_1_to_2 = initTransform(selected_label, samples1, nbrs1, bvh1, samples2, nbrs2, bvh2);
  AffineTransform3 tr_2_to_1 = tr_1_to_2.inverse();

  // Compute forward offsets
  bvh2.setTransform(tr_2_to_1);  // align bounding boxes in first iteration
    updateOffsets(samples1, bvh2, offsets1, selected_label);
  bvh2.clearTransform();

  // Compute backward offsets
  bvh1.setTransform(tr_1_to_2);  // align bounding boxes in first iteration
    updateOffsets(samples2, bvh1, offsets2, selected_label);
  bvh1.clearTransform();

  if (selected_label >= 0)
    cout << "Initialized offsets for label: " << label_name << endl;
  else
    cout << "Initialized offsets for all samples" << endl;
}

void
initOffsets(SampleArray & samples1, NeighborSets const & nbrs1, Bvh & bvh1,
            SampleArray & samples2, NeighborSets const & nbrs2, Bvh & bvh2,
            OffsetArray & offsets1, OffsetArray & offsets2)
{
  // For initial correspondences, we'll match bounding boxes, per-label if labels are available, else globally for the shapes
  if (USE_LABELS)
  {
    for (IndexLabelMap::const_iterator li = LABELS.begin(); li != LABELS.end(); ++li)
      initOffsetsForLabel(li->first, samples1, nbrs1, bvh1, samples2, nbrs2, bvh2, offsets1, offsets2);

    // No label supplied for these samples
    initOffsetsForLabel(0, samples1, nbrs1, bvh1, samples2, nbrs2, bvh2, offsets1, offsets2);
  }
  else
  {
    initOffsetsForLabel(-1, samples1, nbrs2, bvh1, samples2, nbrs2, bvh2, offsets1, offsets2);
  }
}

void
enforceConstraints(SampleArray const & samples1, SampleArray const & samples2, Array<size_t> const & salient_indices1,
                   Array<size_t> const & salient_indices2, OffsetArray & offsets1)
{
  for (size_t i = 0; i < salient_indices1.size(); ++i)
    offsets1[salient_indices1[i]].set(samples2[salient_indices2[i]].p - samples1[salient_indices1[i]].p);
}

// Output offsets assume samples2 is shifted to have the same centroid as samples1
//
// TODO: Do salient points and point labels play nice with each other?
bool
alignNonRigid(SampleArray & samples1, SampleArray & samples2, Array<Vector3> const & salient1,
              Array<Vector3> const & salient2, OffsetArray & offsets1)
{
  // Initialize BVHs
  Bvh bvh1(SamplePtrIterator(samples1.begin()), SamplePtrIterator(samples1.end()));
  bvh1.enableNearestNeighborAcceleration();

  Bvh bvh2(SamplePtrIterator(samples2.begin()), SamplePtrIterator(samples2.end()));
  bvh2.enableNearestNeighborAcceleration();

  THEA_CONSOLE << "Initialized BVHs";

  // Compute sample neighbors
  NeighborSets nbrs1, nbrs2;
  findNeighbors(samples1, bvh1, nbrs1);
  findNeighbors(samples2, bvh2, nbrs2);

  // Map salient points to nearest samples
  if (salient1.size() != salient2.size())
  {
    THEA_ERROR << "Numbers of salient points in two sets don't match";
    return false;
  }

  Array<size_t> salient_indices1, salient_indices2;
  for (size_t i = 0; i < salient1.size(); ++i)
  {
    intx nn_index = bvh1.closestElement<MetricL2>(salient1[i]);
    if (nn_index < 0)
    {
      THEA_ERROR << "Could not map salient1[" << i << "] to a corresponding sample";
      return false;
    }
    salient_indices1.push_back((size_t)nn_index);
  }

  for (size_t i = 0; i < salient2.size(); ++i)
  {
    intx nn_index = bvh2.closestElement<MetricL2>(salient2[i]);
    if (nn_index < 0)
    {
      THEA_ERROR << "Could not map salient2[" << i << "] to a corresponding sample";
      return false;
    }
    salient_indices2.push_back((size_t)nn_index);
  }

  // Start with all offsets zero
  offsets1.clear();
  offsets1.resize(samples1.size());
  OffsetArray offsets2(samples2.size());

  // Initialize offsets per-label if available, else globally
  initOffsets(samples1, nbrs1, bvh1, samples2, nbrs2, bvh2, offsets1, offsets2);

  SampleArray offset_samples2;
  for (int round = 0; round < MAX_ROUNDS; ++round)
  {
    cout << "Round " << round << ':' << flush;

    // Compute forward offsets
    if (round > 0)
      updateOffsets(samples1, bvh2, offsets1);

    // Enforce hard constraints
    enforceConstraints(samples1, samples2, salient_indices1, salient_indices2, offsets1);

    if (round < MAX_ROUNDS - 1)  // don't blend backward offsets on the last round
    {
      // Compute backward offsets
      if (round > 0)
        updateOffsets(samples2, bvh1, offsets2);

      // Smooth backward offsets
      enforceConstraints(samples2, samples1, salient_indices2, salient_indices1, offsets2);
        smoothOffsets(samples2, nbrs2, offsets2, MAX_SMOOTHING_ROUNDS);
      enforceConstraints(samples2, samples1, salient_indices2, salient_indices1, offsets2);

      // Blend backward offsets with forward offsets
      offset_samples2 = samples2;
      for (size_t i = 0; i < samples2.size(); ++i)
        offset_samples2[i].p += offsets2[i].d();

      Bvh offset_bvh2(SamplePtrIterator(offset_samples2.begin()), SamplePtrIterator(offset_samples2.end()));
      offset_bvh2.enableNearestNeighborAcceleration();

      for (size_t i = 0; i < samples1.size(); ++i)
      {
        NNFilter filter(samples1[i].n, samples1[i].label);
        offset_bvh2.pushFilter(&filter);
          intx nn_index = offset_bvh2.closestElement<MetricL2>(samples1[i].p);
        offset_bvh2.popFilter();

        if (nn_index >= 0)
          offsets1[i].set((Real)0.5 * (offsets1[i].d() - offsets2[(size_t)nn_index].d()));
      }
    }

    // Smooth final forward offsets
    if (round < MAX_ROUNDS - 1)
    {
      smoothOffsets(samples1, nbrs1, offsets1, MAX_SMOOTHING_ROUNDS);
      enforceConstraints(samples1, samples2, salient_indices1, salient_indices2, offsets1);
    }
    else
      smoothOffsets(samples1, nbrs1, offsets1, LAST_SMOOTHING_ROUNDS);

    cout << endl;
  }

  return true;
}

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [OPTIONS] <ptfile1> <ptfile2> <offset-outfile1>";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  [--rounds <num-rounds>]";
  THEA_CONSOLE << "  [--smooth <num-rounds>]";
  THEA_CONSOLE << "  [--smooth-last <num-rounds>]";
  THEA_CONSOLE << "  [--no-normals]";
  THEA_CONSOLE << "  [--labels]";
  THEA_CONSOLE << "  [--salient <file1> <file2>]";
  THEA_CONSOLE << "  [--max-salient <num-points>]";
  THEA_CONSOLE << "  [--salient-exclude-prefix <prefix>]+";
  THEA_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 4)
    return usage(argc, argv);

  string samples_path1;
  string samples_path2;
  string offsets_path1;
  string salient_path1, salient_path2;
  intx max_salient = -1;
  Array<string> salient_exclude_prefixes;  // all salient points with tags with these prefixes will be ignored

  int positional = 0;

  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];

    if (arg == "--rounds")
    {
      if (i >= argc - 1)
      {
        THEA_ERROR << "--rounds requires 1 argument";
        return -1;
      }

      MAX_ROUNDS = atoi(argv[++i]);
      if (MAX_ROUNDS < 1)
      {
        THEA_ERROR << "Invalid number of iterations: " << MAX_ROUNDS;
        return -1;
      }
    }
    else if (arg == "--smooth")
    {
      if (i >= argc - 1)
      {
        THEA_ERROR << "--smooth requires 1 argument";
        return -1;
      }

      MAX_SMOOTHING_ROUNDS = atoi(argv[++i]);
      if (MAX_SMOOTHING_ROUNDS < 0)
      {
        THEA_ERROR << "Invalid number of smoothing iterations: " << MAX_SMOOTHING_ROUNDS;
        return -1;
      }
    }
    else if (arg == "--smooth-last")
    {
      if (i >= argc - 1)
      {
        THEA_ERROR << "--smooth-last requires 1 argument";
        return -1;
      }

      LAST_SMOOTHING_ROUNDS = atoi(argv[++i]);
      if (LAST_SMOOTHING_ROUNDS < 0)
      {
        THEA_ERROR << "Invalid number of last-round smoothing iterations: " << LAST_SMOOTHING_ROUNDS;
        return -1;
      }
    }
    else if (arg == "--labels")
    {
      USE_LABELS = true;
    }
    else if (arg == "--no-normals")
    {
      USE_NORMALS = false;
    }
    else if (arg == "--salient")
    {
      if (i >= argc - 2)
      {
        THEA_ERROR << "--salient requires 2 arguments";
        return -1;
      }

      salient_path1 = argv[++i];
      salient_path2 = argv[++i];
    }
    else if (arg == "--max-salient")
    {
      if (i >= argc - 1)
      {
        THEA_ERROR << "--max-salient requires 1 argument";
        return -1;
      }

      max_salient = atoi(argv[++i]);
    }
    else if (arg == "--salient-exclude-prefix")
    {
      if (i >= argc - 1)
      {
        THEA_ERROR << "--salient-exclude-prefix requires 1 argument";
        return -1;
      }

      salient_exclude_prefixes.push_back(argv[++i]);
    }
    else
    {
      switch (positional)
      {
        case 0: samples_path1 = arg; break;
        case 1: samples_path2 = arg; break;
        case 2: offsets_path1 = arg; break;
        default: break;
      }

      positional++;
    }
  }

  if (positional < 3)
    return usage(argc, argv);

  THEA_CONSOLE << "Max iterations: " << MAX_ROUNDS;
  THEA_CONSOLE << "Max smoothing iterations: " << MAX_SMOOTHING_ROUNDS;
  THEA_CONSOLE << "Smoothing iterations in last iteration: " << LAST_SMOOTHING_ROUNDS;
  THEA_CONSOLE << "Match normal orientations: " << USE_NORMALS;
  THEA_CONSOLE << "Match labels: " << USE_LABELS;
  THEA_CONSOLE << "Max salient points: " << max_salient;

  SampleArray samples1, samples2;
  if (!loadSamples(samples_path1, samples1, USE_LABELS) || !loadSamples(samples_path2, samples2, USE_LABELS))
    return -1;

  bool has_salient = (!salient_path1.empty() && !salient_path2.empty());
  Array<Vector3> salient1, salient2;
  if (has_salient)
  {
    ifstream in1(salient_path1.c_str());
    if (!in1)
    {
      THEA_ERROR << "Could not open file '" << salient_path1 << "' containing first set of salient points";
      return -1;
    }

    ifstream in2(salient_path2.c_str());
    if (!in2)
    {
      THEA_ERROR << "Could not open file '" << salient_path2 << "' containing second set of salient points";
      return -1;
    }

    Array<string> salient_lines1, salient_lines2;
    string line;
    while (getline(in1, line)) salient_lines1.push_back(line);
    while (getline(in2, line)) salient_lines2.push_back(line);

    for (size_t i = 0; i < salient_lines1.size(); ++i)
    {
      if (max_salient >= 0 && (intx)salient1.size() >= max_salient) break;

      if (salient_lines1[i].empty())
        continue;

      istringstream line_in1(salient_lines1[i]);
      string tag1; line_in1 >> tag1;
      if (tag1.empty()) continue;

      bool exclude = false;
      for (size_t j = 0; j < salient_exclude_prefixes.size(); ++j)
        if (beginsWith(tag1, salient_exclude_prefixes[j]))
        {
          exclude = true;
          break;
        }

      if (exclude)
        continue;

      for (size_t j = 0; j < salient_lines2.size(); ++j)
      {
        istringstream line_in2(salient_lines2[j]);
        string tag2; line_in2 >> tag2;
        if (tag1 == tag2)
        {
          Vector3 s1, s2;
          if ((line_in1 >> s1[0] >> s1[1] >> s1[2])
           && (line_in2 >> s2[0] >> s2[1] >> s2[2]))
          {
            salient1.push_back(s1);
            salient2.push_back(s2);
            break;
          }
        }
      }
    }

    THEA_CONSOLE << "Read " << salient1.size() << " salient point(s)";
  }

  OffsetArray offsets1;
  if (!alignNonRigid(samples1, samples2, salient1, salient2, offsets1))
    return -1;

  if (offsets1.size() != samples1.size())
  {
    THEA_ERROR << "Incorrect number of offsets produced (expected " << samples1.size() << ", computed " << offsets1.size()
               << ')';
    return -1;
  }

  //===========================================================================================================================
  // Save offsets, 3 numbers per line
  //===========================================================================================================================
  {
    ofstream out(offsets_path1.c_str(), ios::binary);
    if (!out)
    {
      THEA_ERROR << "Could not open output file '" << offsets_path1 << "' for writing";
      return -1;
    }

    for (size_t i = 0; i < offsets1.size(); ++i)
    {
      Vector3 const & d = offsets1[i].d();
      out << d[0] << ' ' << d[1] << ' ' << d[2] << endl;
    }

    THEA_CONSOLE << "Wrote offsets to " << offsets_path1;
  }

  //===========================================================================================================================
  // Save source points plus offsets, 6 numbers per line
  //===========================================================================================================================
  {
    string pts_with_offsets_path1 = FilePath::concat(FilePath::parent(offsets_path1),
                                                     FilePath::baseName(samples_path1) + "_with_offsets.pts");
    ofstream out_pts(pts_with_offsets_path1.c_str(), ios::binary);
    for (size_t i = 0; i < samples1.size(); ++i)
    {
      Vector3 const & p = samples1[i].p;
      Vector3 const & n = offsets1[i].d();
      out_pts << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
              << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
    }

    THEA_CONSOLE << "Wrote points with offsets to " << pts_with_offsets_path1;
  }

  //===========================================================================================================================
  // Save deformed source points plus labels, and the labels individually
  //===========================================================================================================================
  {
    string offset_pts_path1 = FilePath::concat(FilePath::parent(offsets_path1),
                                               FilePath::baseName(samples_path1) + "_deformed.pts");
    ofstream out_def_pts(offset_pts_path1.c_str(), ios::binary);

    ofstream out_def_labels;
    if (USE_LABELS)
    {
      string offset_labels_path1 = FilePath::changeExtension(offset_pts_path1, "labels");
      out_def_labels.open(offset_labels_path1.c_str(), ios::binary);
    }

    for (size_t i = 0; i < samples1.size(); ++i)
    {
      Vector3 p = samples1[i].p + offsets1[i].d();
      Vector3 const & n = samples1[i].n;
      out_def_pts << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
                  << n[0] << ' ' << n[1] << ' ' << n[2];

      if (USE_LABELS)
      {
        int32 label_hash = samples1[i].label;
        if (label_hash >= 0)
        {
          IndexLabelMap::const_iterator existing = LABELS.find(label_hash);
          if (existing != LABELS.end())
          {
            out_def_pts << " \"" << existing->second << "\"\n";
            out_def_labels << existing->second << '\n';
          }
          else
          {
            out_def_pts << " \"\"\n";
            out_def_labels << '\n';
          }
        }
        else
        {
          out_def_pts << " \"\"\n";
          out_def_labels << '\n';
        }
      }
      else
      {
        out_def_pts << '\n';
        out_def_labels << '\n';
      }
    }

    THEA_CONSOLE << "Wrote offset (deformed) points to " << offset_pts_path1;
  }

  //===========================================================================================================================
  // Save source and target points, colored by target point position
  //===========================================================================================================================
  {
    string colored_pts_path1 = FilePath::concat(FilePath::parent(offsets_path1),
                                                FilePath::baseName(samples_path1) + "_colored.pts");
    ofstream out_colored_pts1(colored_pts_path1.c_str(), ios::binary);

    AxisAlignedBox3 bounds2 = computeBounds(samples2);
    Vector3 center2 = bounds2.getCenter(), half_ext2 = 0.5f * bounds2.getExtent();

    for (size_t i = 0; i < samples1.size(); ++i)
    {
      Vector3 p = samples1[i].p;
      Vector3 p_def = p + offsets1[i].d();
      Vector3 n = (p_def - center2).cwiseQuotient(half_ext2);
      out_colored_pts1 << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
                       << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
    }

    out_colored_pts1.close();

    string colored_pts_path2 = FilePath::concat(FilePath::parent(offsets_path1),
                                                FilePath::baseName(samples_path2) + "_colored.pts");
    ofstream out_colored_pts2(colored_pts_path2.c_str(), ios::binary);

    for (size_t i = 0; i < samples2.size(); ++i)
    {
      Vector3 const & p = samples2[i].p;
      Vector3 n = (p - center2).cwiseQuotient(half_ext2);
      out_colored_pts2 << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
                       << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
    }

    THEA_CONSOLE << "Wrote points colored by correspondences to " << colored_pts_path1 << " and " << colored_pts_path2;
  }


  //===========================================================================================================================
  // Write pairs of corresponding points
  //===========================================================================================================================
  {
    string corr_path = FilePath::concat(FilePath::parent(offsets_path1),
                           FilePath::baseName(samples_path1) + "_corr_" + FilePath::baseName(samples_path2) + ".txt");
    ofstream out_corr(corr_path.c_str(), ios::binary);

    Bvh bvh2(SamplePtrIterator(samples2.begin()), SamplePtrIterator(samples2.end()));
    bvh2.enableNearestNeighborAcceleration();

    intx num_matched = 0;
    double sum_sqdist = 0;
    for (size_t i = 0; i < samples1.size(); ++i)
    {
      Vector3 p1 = samples1[i].p;

      Vector3 deformed_p1 = p1 + offsets1[i].d();
      NNFilterLabelOnly filter(samples1[i].label);
      bvh2.pushFilter(&filter);
        intx nn_index = bvh2.closestElement<MetricL2>(deformed_p1);
      bvh2.popFilter();

      if (nn_index < 0)
      {
        // THEA_WARNING << "No corresponding target point found for source point " << toString(p1);
        continue;
      }

      out_corr << i << ' ' << nn_index << '\n';

      num_matched++;
      sum_sqdist += (deformed_p1 - samples2[(size_t)nn_index].p).squaredNorm();
    }

    THEA_CONSOLE << "Wrote correspondences for " << num_matched << '/' << samples1.size() << " point(s) to " << corr_path;
    THEA_CONSOLE << "Root mean squared separation of matched pairs = "
                 << (num_matched > 0 ? sqrt(sum_sqdist / num_matched) : 0);
  }

  return 0;
}
