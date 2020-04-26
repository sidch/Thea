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

#ifndef __Thea_Algorithms_BestFitBox3_hpp__
#define __Thea_Algorithms_BestFitBox3_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Box3.hpp"

namespace Thea {
namespace Algorithms {

/** Approximate best-fit oriented bounding box. */
class THEA_API BestFitBox3
{
  public:
    THEA_DECL_SMART_POINTERS(BestFitBox3)

    /** Constructor. */
    BestFitBox3();

    /** Add a point to the set. */
    void addPoint(Vector3 const & point);

    /** Remove all data and (lazily) set the box to null. */
    void clear();

    /** Set the up vector. The computed box will only consider orientations with the given up vector. */
    void setUpVector(Vector3 const & up_) { up = up_.normalized(); has_up = true; }

    /** Check if the up vector has been set. */
    bool hasUpVector() const { return has_up; }

    /**
     * Get the up vector, if it has been set.
     *
     * @see hasUpVector();
     */
    Vector3 const & getUpVector() const { return up; }

    /** Clear the up vector. Subsequent alignments will be unconstrained. */
    void clearUpVector() { has_up = false; }

    /** Remove all cached data to free memory, but do <b>not</b> mark the box for recomputation. */
    void releaseMemoryWithoutUpdate();

    /** Get the box object. Will only force a recomputation if new data have been added since the last call to this function. */
    Box3 const & getBox() const;

  private:
    /** Recompute the best-fit oriented bounding box. */
    void update() const;

    Array<Vector3> points;
    bool has_up;
    Vector3 up;

    mutable Box3 box;
    mutable bool updated;

}; // class BestFitBox3

} // namespace Algorithms
} // namespace Thea

#endif
