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
// First version: 2009
//
//============================================================================

#ifndef __Thea_Algorithms_BestFitEllipsoid3_hpp__
#define __Thea_Algorithms_BestFitEllipsoid3_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Box3.hpp"

namespace Thea {
namespace Algorithms {

/** Approximate best-fit ellipsoid. */
class THEA_API BestFitEllipsoid3
{
  public:
    THEA_DECL_SMART_POINTERS(BestFitEllipsoid3)

    /**
     * Constructor. The optional argument (which may be ignored by the implementation) controls the approximation ratio
     * (1 + eps_) of the computed ellipsoid. A smaller value gives a better fit but may take longer to compute.
     */
    BestFitEllipsoid3(Real eps_ = 0.05);

    /** Add a point to the set. */
    void addPoint(Vector3 const & point);

    /** Remove all data and (lazily) set the ellipsoid to null. */
    void clear();

    /** Remove all cached data to free memory, but do <b>not</b> mark the ellipsoid for recomputation. */
    void releaseMemoryWithoutUpdate();

    /**
     * Get the axes of the ellipsoid. The lengths of the axes represent the shape of the ellipsoid. The axes will always be
     * returned in decreasing order of length, i.e. |axis0| >= |axis1| >= |axis2|. Also,
     * (axis0.cross(axis1)).normalized() == axis2.normalized().
     */
    void getAxes(Vector3 & axis0, Vector3 & axis1, Vector3 & axis2) const;

    /** Get the center of the ellipsoid. */
    Vector3 const & getCenter() const;

    /**
     * Get a tight bounding box for the input data (<b>not</b> for the ellipsoid) oriented along the principal axes of the
     * ellipsoid.
     */
    Box3 const & getOrientedBoundingBox() const;

  private:
    /** Recompute the best-fit ellipsoid. */
    void update() const;

    Array<Vector3> points;
    Real eps;
    mutable Vector3 center, axis[3];
    mutable Box3 obb;
    mutable bool updated;

}; // class BestFitEllipsoid3

} // namespace Algorithms
} // namespace Thea

#endif
