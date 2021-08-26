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

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_RandomWalks_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_RandomWalks_hpp__

#include "../../PointCloud3.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/** Compute the average offset, from the query position, after each step of an n-step random walk on a surface. */
class RandomWalks
{
  public:
    /**
     * Constructs the object to compute random walk patterns on a given surface. The sampled surface must persist as long as
     * this object does.
     */
    RandomWalks(PointCloud3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointCloud3 const * getSurface() const { return surf; }

    /**
     * Compute the average offset, from the query position, after each step of an n-step random walk on the shape's sample
     * graph.
     *
     * @param position The position of the query point.
     * @param num_steps The number of random steps to take from the query point.
     * @param features Used to return the point features. Should be preallocated to \a num_steps * 3 entries.
     * @param num_walks The number of random walks over which to take averages.
     */
    void compute(Vector3 const & position, intx num_steps, double * features, intx num_walks = -1) const;

  private:
    PointCloud3 const * surf;  ///< The point-sampled surface.

}; // class RandomWalks

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
