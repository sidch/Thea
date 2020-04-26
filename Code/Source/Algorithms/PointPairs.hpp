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

#ifndef __Thea_Algorithms_PointPairs_hpp__
#define __Thea_Algorithms_PointPairs_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"
#include <utility>

namespace Thea {
namespace Algorithms {

/** A pair of points in 1-space (single-precision) */
typedef std::pair<float, float> FloatPair;

/** A pair of points in 1-space (double-precision) */
typedef std::pair<double, double> DoublePair;

/** A pair of points in 1-space (default precision) */
typedef std::pair<Real, Real> RealPair;

/** A pair of points in 2-space. */
typedef std::pair<Vector2, Vector2> PointPair2;

/** A pair of points in 3-space. */
typedef std::pair<Vector3, Vector3> PointPair3;

/** A pair of points in 4-space. */
typedef std::pair<Vector4, Vector4> PointPair4;

/** A pair of points in n-dimensional space. */
template <size_t N, typename T>
class /* THEA_API */ PointPairN : public std::pair< Vector<N, T>, Vector<N, T> >
{
  private:
    typedef std::pair< Vector<N, T>, Vector<N, T> > BaseType;

  public:
    /** Default constructor. */
    PointPairN() : BaseType() {}

    /** Initializing constructor. */
    PointPairN(Vector<N, T> const & a, Vector<N, T> const & b) : BaseType(a, b) {}

}; // class PointPairN

} // namespace Algorithms
} // namespace Thea

#endif
