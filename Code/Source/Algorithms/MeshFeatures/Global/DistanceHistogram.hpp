//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2015, Siddhartha Chaudhuri
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
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
namespace Global {

/** Compute the histogram of distances between pairs of points on the shape. */
template < typename ExternalSampleKDTreeT = KDTreeN<Vector3, 3> >
class DistanceHistogram
{
  public:
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    static long const DEFAULT_NUM_SAMPLES = 5000;  ///< Default number of points to initially sample from the shape.
    static long const DEFAULT_APPROX_NUM_PAIRS = 500 * 500;  /**< Default (approximate) number of point pairs to use for the
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
    DistanceHistogram(MeshT const & mesh, long num_samples = -1, Real normalization_scale = -1)
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
    DistanceHistogram(Graphics::MeshGroup<MeshT> const & mesh_group, long num_samples = -1, Real normalization_scale = -1)
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
    DistanceHistogram(ExternalSampleKDTree const * sample_kdtree, Real normalization_scale = -1)
    : ldh(sample_kdtree, normalization_scale)
    {}

    /** Get the number of surface samples. */
    long numSamples() const { return ldh.numSamples(); }

    /** Get the position of the surface sample with index \a index. */
    Vector3 getSamplePosition(long index) const { return ldh.getSamplePosition(index); }

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
      long num_samples = ldh.numSamples();
      long num_distinct_unordered = (num_samples - 1) * num_samples;

      if (pair_reduction_ratio < 0)
        pair_reduction_ratio = (Real)std::min(0.5, 1000000.0 / num_distinct_unordered);  // don't count (x, x)

      histogram.setZero();  // don't bother setting the range

      Real local_reduction_ratio = std::sqrt(pair_reduction_ratio);
      long num_queries = Math::clamp((long)std::ceil(local_reduction_ratio * num_samples), 0, num_samples - 1);
      if (num_queries <= 0)
        return;

      TheaArray<int32> query_indices((size_t)num_queries);
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
    Local::LocalDistanceHistogram<ExternalSampleKDTree> ldh;  ///< Used to compute histograms from a single query point.

}; // class DistanceHistogram

} // namespace Global
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
