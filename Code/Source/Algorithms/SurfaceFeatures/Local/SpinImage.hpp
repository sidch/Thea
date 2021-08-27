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

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_SpinImage_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_SpinImage_hpp__

#include "../../PointSet3.hpp"
#include "../../../MatVec.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/**
 * Compute the spin image at a point on a shape. The spin image at a point is a histogram of the rest of the surface, where each
 * bin corresponds to an annular ring with the normal at the central point as its axis. The bins are parametrized by radial
 * distance from this axis, and height above/below the tangent plane.
 *
 * A. Johnson and M. Hebert, "Using Spin Images for Efficient Object Recognition in Cluttered 3D Scenes",
 * IEEE Trans. PAMI 21(5), 433-449, 1999.
 */
class SpinImage
{
  public:
    /**
     * Constructs the object to compute spin images at points on a given surface. The sampled surface must persist as long as
     * this object does.
     */
    SpinImage(PointSet3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointSet3 const * getSurface() const { return surf; }

    /**
     * Compute the spin image at a query point on the mesh.
     *
     * This version of the function explicitly computes the normal at the query point -- the other version of the function
     * should be used if the normal is known in advance. The normal is computed from samples, so <b>may be quite inaccurate</b>
     * especially in thin areas.
     *
     * @param position The point at which the spin image is to be calculated.
     * @param num_radial_bins The number of bins in the radial direction (distance from axis).
     * @param num_height_bins The number of bins in the height direction (distance from tangent plane).
     * @param spin_image Used to return the computed spin image, with rows corresponding to radial divisions and columns to
     *   height divisions.
     */
    void compute(Vector3 const & position, int num_radial_bins, int num_height_bins, MatrixX<double> & spin_image) const;

    /**
     * Compute the spin image at a query point on the mesh with a known normal.
     *
     * @param position The point at which the spin image is to be calculated.
     * @param normal The normal at this point.
     * @param num_radial_bins The number of bins in the radial direction (distance from axis).
     * @param num_height_bins The number of bins in the height direction (distance from tangent plane).
     * @param spin_image Used to return the computed spin image, with rows corresponding to radial divisions and columns to
     *   height divisions.
     */
    void compute(Vector3 const & position, Vector3 const & normal, int num_radial_bins, int num_height_bins,
                 MatrixX<double> & spin_image) const;

  private:
    PointSet3 const * surf;  ///< The point-sampled surface.

}; // class SpinImage

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
