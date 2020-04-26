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
// First version: 2018
//
//============================================================================

#ifndef __Thea_Algorithms_FurthestPointSampling_hpp__
#define __Thea_Algorithms_FurthestPointSampling_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Subsample a set of 3D points to pick a subset of points distributed evenly, that is, neighboring points are (roughly) equally
 * spaced. The algorithm repeatedly picks the left-over point that is furthest from any of the previously selected ones. The
 * returned list of points has the property that any prefix of the list is also an evenly spaced set, making further subsampling
 * trivial.
 */
class FurthestPointSampling
{
  public:
    /**
     * Samples approximately uniformly separated points from a larger set. The algorithm repeatedly picks the left-over point
     * that is furthest from any of the previously selected ones. The function returns an array of points ordered so that for
     * any K, the first K points form an approximately uniformly separated subsampling.
     *
     * @param num_orig_points The number of input points.
     * @param orig_points The set of input points.
     * @param num_desired_points The number of points to be subsampled.
     * @param selected_indices The indices of the subsampled points. This array must be preallocated to (at least)
     *   \a num_desired_points elements.
     * @param dist_type The distance metric to be used. Currently only DistanceType::GEODESIC is supported.
     * @param verbose If true, prints progress messages.
     *
     * @return The number of subsampled points. A negative value, or a value less than \a num_desired_points in general,
     *   indicates an error occurred.
     */
    static intx subsample(intx num_orig_points, Vector3 const * orig_points, intx num_desired_points, intx * selected_indices,
                          DistanceType dist_type = DistanceType::GEODESIC, bool verbose = false);

}; // class FurthestPointSampling

} // namespace Algorithms
} // namespace Thea

#endif
