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
// First version: 2016
//
//============================================================================

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_AverageDistance_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_AverageDistance_hpp__

#include "../../PointCloud3.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/** Compute the average distance from a query point to other points on a shape. */
class AverageDistance
{
  public:
    /**
     * Constructs the object to compute the histogram of distances to sample points on a given surface. The sampled surface must
     * persist as long as this object does.
     */
    AverageDistance(PointCloud3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointCloud3 const * getSurface() const { return surf; }

    /**
     * Compute the average distance from a query point to sample points on the shape. The returned distance is normalized by
     * dividing by the normalization scale, which is \a max_distance (if non-negative), else the shape's native scale.
     *
     * @param position The position of the query point.
     * @param dist_type The type of distance metric.
     * @param max_distance The distance to the furthest point to consider. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor).
     */
    double compute(Vector3 const & position, DistanceType dist_type = DistanceType::EUCLIDEAN, Real max_distance = -1) const;

  private:
    PointCloud3 const * surf;  ///< The point-sampled surface.

}; // class AverageDistance

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
