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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_LocalDistanceHistogram_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_LocalDistanceHistogram_hpp__

#include "../../PointCloud3.hpp"
#include "../../Histogram.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/** Compute the histogram of distances from a query point to other points on a shape. */
class LocalDistanceHistogram
{
  public:
    /**
     * Constructs the object to compute the histogram of distances to sample points on a given surface. The sampled surface must
     * persist as long as this object does.
     */
    LocalDistanceHistogram(PointCloud3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointCloud3 const * getSurface() const { return surf; }

    /**
     * Compute the histogram of distances from a query point to sample points on the shape. The histogram bins uniformly
     * subdivide the range of distances from zero to \a max_distance. If \a max_distance is negative, the shape's native scale
     * will be used.
     *
     * @param position The position of the query point.
     * @param histogram The histogram to be computed.
     * @param dist_type The type of distance metric.
     * @param max_distance The distance to the furthest point to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale. The histogram range is set
     *   appropriately.
     * @param sample_reduction_ratio The fraction of the available set of samples -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. This may be useful for getting a more evenly
     *   sampled set of pairwise distances when calling this function with multiple query points (and an extra-large initial set
     *   of points). A negative value, or a value of 1, indicates all sample points will be used.
     */
    void compute(Vector3 const & position, Histogram & histogram, DistanceType dist_type = DistanceType::EUCLIDEAN,
                 Real max_distance = -1, Real sample_reduction_ratio = -1) const;

  private:
    PointCloud3 const * surf;  ///< The point-sampled surface.

}; // class LocalDistanceHistogram

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
