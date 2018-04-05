//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2018, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_FurthestPointSampling_hpp__
#define __Thea_Algorithms_FurthestPointSampling_hpp__

#include "../Common.hpp"
#include "../Vector3.hpp"

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
    static long subsample(long num_orig_points, Vector3 const * orig_points, long num_desired_points, long * selected_indices,
                          DistanceType dist_type = DistanceType::GEODESIC, bool verbose = false);

}; // class FurthestPointSampling

} // namespace Algorithms
} // namespace Thea

#endif
