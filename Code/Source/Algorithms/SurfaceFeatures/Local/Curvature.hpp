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
// First version: 2013
//
//============================================================================

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_Curvature_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_Curvature_hpp__

#include "../../PointCloud3.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/** Compute curvature measures at a point on a shape. */
class Curvature
{
  public:
    /**
     * Constructs the object to compute curvature measures at points on a given surface. The sampled surface must persist as
     * long as this object does.
     */
    Curvature(PointCloud3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointCloud3 const * getSurface() const { return surf; }

    /**
     * Compute the <em>projected</em> curvature at a query point on the mesh. The projected curvature is an approximation to the
     * actual curvature, obtained by projecting sample points in the neighborhood of the query point onto the normal.
     *
     * This version of the function explicitly computes the normal at the query point -- the other version of the function
     * should be used if the normal is known in advance. The normal is computed from the nearest sample, so <b>may be quite
     * inaccurate</b> especially in thin areas.
     *
     * @param position The position at which to compute curvature.
     * @param nbd_radius The size of the neighborhood over which to compute curvature, specified as a multiple of the shape
     *   scale. A negative argument selects a default size.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal, distinguishing
     *   convex from concave. This is completely different from the sign of the Gaussian curvature.
     */
    double computeProjectedCurvature(Vector3 const & position, Real nbd_radius = -1) const;

    /**
     * Compute the <em>projected</em> curvature at a query point with a known normal on the mesh. The projected curvature is an
     * approximation to the actual curvature, obtained by projecting sample points in the neighborhood of the query point onto
     * the normal.
     *
     * @param position The position at which to compute curvature.
     * @param normal The normal at this position.
     * @param nbd_radius The size of the neighborhood over which to compute curvature, specified as a multiple of the shape
     *   scale. A negative argument selects a default size.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal, distinguishing
     *   convex from concave. This is completely different from the sign of the Gaussian curvature.
     */
    double computeProjectedCurvature(Vector3 const & position, Vector3 const & normal, Real nbd_radius = -1) const;

  private:
    PointCloud3 const * surf;  ///< The point-sampled surface.

}; // class Curvature

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
