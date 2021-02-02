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
// First version: 2011
//
//============================================================================

#ifndef __Thea_LineN_hpp__
#define __Thea_LineN_hpp__

#include "Common.hpp"
#include "Math.hpp"
#include "MatVec.hpp"

namespace Thea {

// Forward declarations
template <int N, typename T> class LineN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for straight lines in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and
 * T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ LineNBase
{
  public:
    typedef LineN<N, T>   LineT;    ///< N-dimensional straight line.
    typedef Vector<N, T>  VectorT;  ///< N-dimensional vector.

    THEA_DECL_SMART_POINTERS(LineT)

    /**
     * Construct a line from a point on it, and the direction vector of the line (need not be a unit vector). The \a normalize
     * argument suppresses rescaling of the direction vector to unit length if set to false.
     */
    static LineT fromPointAndDirection(VectorT const & point_, VectorT const & direction_, bool normalize = true)
    {
      if (Math::fuzzyEq(direction_.squaredNorm(), static_cast<T>(0)))
        throw Error("LineN: Direction vector has zero (or nearly zero) length");

      LineT line;
      line.point = point_;
      line.direction = direction_; if (normalize) line.direction.normalize();
      return line;
    }

    /** Construct a line from two points on it. */
    static LineT fromTwoPoints(VectorT const & point1, VectorT const & point2)
    {
      return fromPointAndDirection(point1, point2 - point1);
    }

    /** Cast the line to a different scalar type. */
    template <typename U> LineN<N, U> cast() const
    {
      return LineN<N, U>::fromPointAndDirection(point.template cast<U>(), direction.template cast<U>(),
                                                /* normalize = */ false);
    }

    /** Get a point on the line. */
    VectorT const & getPoint() const { return point; }

    /** Get the unit direction vector of the line. */
    VectorT const & getDirection() const { return direction; }

    /** Get the distance of the line from a given point. */
    T distance(VectorT const & p) const
    {
      return std::sqrt(squaredDistance(p));
    }

    /** Get the square of the distance of the line from a given point. */
    T squaredDistance(VectorT const & p) const
    {
      return (p - closestPoint(p)).squaredNorm();
    }

    /** Get the point on the line closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      T t = direction.dot(p - point);
      return point + t * direction;
    }

    /** Get the distance of this line from another line. */
    T distance(LineT const & other) const
    {
      return std::sqrt(squaredDistance(other));
    }

    /**
     * Get the point on this line and the point on another line closest to each other, and return the squared distance between
     * them.
     */
    T squaredDistance(LineT const & other, VectorT * this_pt = nullptr, VectorT * other_pt = nullptr) const
    {
      // Adapted from Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.

      VectorT r = point - other.point;
      T b = direction.dot(other.direction);
      T f = other.direction.dot(r);
      T denom = 1 - b * b; // always nonnegative

      // If segments not parallel, compute closest point on L1 to L2. Else pick arbitrary s (here 0).
      T s = 0;
      if (Math::fuzzyGt(denom, static_cast<T>(0)))
      {
        T c = direction.dot(r);
        s = (b * f - c) / denom;
      }

      VectorT c1 = point + s * direction;
      if (this_pt)
        *this_pt = c1;

      // Compute point on L2 closest to S1(s) using t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f)
      T t = b * s + f;
      VectorT c2 = other.point + t * direction;
      if (other_pt)
        *other_pt = c2;

      return (c1 - c2).squaredNorm();
    }

    /** Get a textual description of the line. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[P: " << Thea::toString(point) << ", U: " << Thea::toString(direction) << ']';
      return oss.str();
    }

  private:
    VectorT point;      ///< A point on the line.
    VectorT direction;  ///< A unit vector along the direction of the line.

}; // class LineNBase

} // namespace Internal

/** A straight line in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <int N, typename T = Real>
class /* THEA_API */ LineN : public Internal::LineNBase<N, T>
{
}; // class LineN

} // namespace Thea

#include "Line2.hpp"
#include "Line3.hpp"

#endif
