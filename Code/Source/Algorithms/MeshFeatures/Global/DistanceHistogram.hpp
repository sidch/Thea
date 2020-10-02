//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2015
//
//============================================================================

#ifndef __Thea_Algorithms_MeshFeatures_Global_DistanceHistogram_hpp__
#define __Thea_Algorithms_MeshFeatures_Global_DistanceHistogram_hpp__

#include "../Local/LocalDistanceHistogram.hpp"
#include "../../Histogram.hpp"
#include "../../../Math.hpp"
#include "../../../Random.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {

/** Namespace for classes that compute global features on a mesh. */
namespace Global {

/** Compute the histogram of distances between pairs of points on the shape. */
template < typename ExternalSampleKdTreeT = KdTreeN<Vector3, 3> >
class DistanceHistogram
{
  public:
    typedef ExternalSampleKdTreeT ExternalSampleKdTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    static intx const DEFAULT_NUM_SAMPLES = 5000;  ///< Default number of points to initially sample from the shape.
    static intx const DEFAULT_APPROX_NUM_PAIRS = 500 * 500;  /**< Default (approximate) number of point pairs to use for the
                                                                  histogram. */

  public:
    /**
     * Constructs the object to compute the histogram of distances between sample points on a given mesh. The mesh need not
     * persist beyond the execution of the constructor.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    DistanceHistogram(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : ldh(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances between sample points on a given mesh group. The mesh group
     * need not persist beyond the execution of the constructor.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    DistanceHistogram(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : ldh(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances between sample points of a shape with a precomputed kd-tree
     * on these points. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree The kd-tree representing the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    DistanceHistogram(ExternalSampleKdTree const * sample_kdtree, Real normalization_scale = -1)
    : ldh(sample_kdtree, normalization_scale)
    {}

    /** Get the number of surface samples. */
    intx numSamples() const { return ldh.numSamples(); }

    /** Get the position of the surface sample with index \a index. */
    Vector3 getSamplePosition(intx index) const { return ldh.getSamplePosition(index); }

    /**
     * Compute the histogram of distances between sample points on the shape. The histogram bins uniformly subdivide the range
     * of distances from zero to \a max_distance. If \a max_distance is negative, the shape scale specified in the constructor
     * will be used.
     *
     * @param histogram The histogram to be computed.
     * @param dist_type The type of distance metric.
     * @param max_distance The maximum separation between points to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor). The histogram range is set appropriately.
     * @param pair_reduction_ratio The fraction of the available set of point pairs -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. Note that the final set is a subsampling of
     *   the set of all possible pairs, not the set of all possible pairs of a subsampled set of points. This may be useful for
     *   getting a more evenly sampled set of pairwise distances, with an extra-large initial set of points. A value of 1
     *   indicates all ordered pairs will be used, but this counts every distance twice, so a value of 0.5 or so may be more
     *   appropriate. A negative value picks a default of 0.5, or the ratio that gives a maximum of ~1M ordered pairs, whichever
     *   is smaller.
     */
    void compute(Histogram & histogram, DistanceType dist_type, Real max_distance = -1, Real pair_reduction_ratio = -1) const
    {
      intx num_samples = ldh.numSamples();
      intx num_distinct_unordered = (num_samples - 1) * num_samples;

      if (pair_reduction_ratio < 0)
        pair_reduction_ratio = (Real)std::min(0.5, 1000000.0 / num_distinct_unordered);  // don't count (x, x)

      histogram.setZero();  // don't bother setting the range

      Real local_reduction_ratio = std::sqrt(pair_reduction_ratio);
      intx num_queries = Math::clamp((intx)std::ceil(local_reduction_ratio * num_samples), 0, num_samples - 1);
      if (num_queries <= 0)
        return;

      Array<int32> query_indices((size_t)num_queries);
      Random::common().sortedIntegers(0, (int32)num_samples - 1, (int32)num_queries, &query_indices[0]);

      Histogram local_histogram(histogram.numBins());
      for (size_t i = 0; i < query_indices.size(); ++i)
      {
        Vector3 p = ldh.getSamplePosition(query_indices[i]);
        ldh.compute(p, local_histogram, dist_type, max_distance, local_reduction_ratio);

        // Remove the zero distance from the query point to itself
        local_histogram.remove(0.0);

        if (i == 0)
          histogram.setRange(local_histogram.minValue(), local_histogram.maxValue());

        histogram.insert(local_histogram);
      }
    }

  private:
    Local::LocalDistanceHistogram<ExternalSampleKdTree> ldh;  ///< Used to compute histograms from a single query point.

}; // class DistanceHistogram

} // namespace Global
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
