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

#ifndef __Thea_Algorithms_SurfaceFeatures_Global_DistanceHistogram_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Global_DistanceHistogram_hpp__

#include "../Local/LocalDistanceHistogram.hpp"
#include "../../Histogram.hpp"

namespace Thea {
namespace Algorithms {

/** Namespace for classes that compute features of a surface. */
namespace SurfaceFeatures {

/** Namespace for classes that compute global features of a surface. */
namespace Global {

/** Compute the histogram of distances between pairs of points on a shape. */
class DistanceHistogram
{
  public:
    /**
     * Constructs the object to compute the histogram of distances between sample points on a given surface. The sampled surface
     * must persist as long as this object does.
     */
    DistanceHistogram(PointSet3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointSet3 const * getSurface() const { return ldh.getSurface(); }

    /**
     * Compute the histogram of distances between sample points on the shape. The histogram bins uniformly subdivide the range
     * of distances from zero to \a max_distance. If \a max_distance is negative, the shape's native scale is used.
     *
     * @param histogram The histogram to be computed.
     * @param dist_type The type of distance metric.
     * @param max_distance The maximum separation between points to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale). The histogram range is set
     *   appropriately.
     * @param pair_reduction_ratio The fraction of the available set of point pairs -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. Note that the final set is a subsampling of
     *   the set of all possible pairs, not the set of all possible pairs of a subsampled set of points. This may be useful for
     *   getting a more evenly sampled set of pairwise distances, with an extra-large initial set of points. A value of 1
     *   indicates all ordered pairs will be used, but this counts every distance twice, so a value of 0.5 or so may be more
     *   appropriate. A negative value picks a default of 0.5, or the ratio that gives a maximum of ~1M ordered pairs, whichever
     *   is smaller.
     */
    void compute(Histogram & histogram, DistanceType dist_type, Real max_distance = -1, Real pair_reduction_ratio = -1) const;

  private:
    Local::LocalDistanceHistogram ldh;  ///< Used to compute histograms from a single query point.

}; // class DistanceHistogram

} // namespace Global
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
