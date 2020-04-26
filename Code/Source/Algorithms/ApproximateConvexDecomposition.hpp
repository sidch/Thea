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

#ifndef __Thea_Algorithms_ApproximateConvexDecomposition_hpp__
#define __Thea_Algorithms_ApproximateConvexDecomposition_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

/** Approximate convex decomposition of point clouds. */
class THEA_API ApproximateConvexPointDecomposition
{
  public:
    THEA_DECL_SMART_POINTERS(ApproximateConvexPointDecomposition)

    /** Add a point to the cloud to be decomposed. */
    void addPoint(Vector3 const & p) { points.push_back(p); }

    /** Add a set of points to the cloud to be decomposed. */
    template <typename PointInputIterator> void addPoint(PointInputIterator begin, PointInputIterator end)
    { points.insert(points.end(), begin, end); }

    /** Compute the approximate convex decomposition. */
    void decompose(double convexity_threshold, Array< Array<int> > & results);

  private:
    Array<Vector3> points;

}; // class ApproximateConvexPointDecomposition

} // namespace Algorithms
} // namespace Thea

#endif
