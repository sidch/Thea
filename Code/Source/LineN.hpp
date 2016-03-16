//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_LineN_hpp__
#define __Thea_LineN_hpp__

#include "Common.hpp"
#include "Math.hpp"
#include "VectorN.hpp"

namespace Thea {

// Forward declarations
template <long N, typename T> class LineN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for straight lines in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and
 * T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ LineNBase
{
  public:
    typedef LineN<N, T>    LineT;    ///< N-dimensional straight line.
    typedef VectorN<N, T>  VectorT;  ///< N-dimensional vector.

    THEA_DEF_POINTER_TYPES(LineT, shared_ptr, weak_ptr)

    /** Construct a line from a point on it, and the direction vector of the line (need not be a unit vector). */
    static LineT fromPointAndDirection(VectorT const & point_, VectorT const & direction_)
    {
      if (Math::fuzzyEq(direction_.squaredLength(), static_cast<T>(0)))
        throw Error("LineN: Direction vector has zero (or nearly zero) length");

      LineT line;
      line.point = point_;
      line.direction = direction_.unit();
      return line;
    }

    /** Construct a line from two points on it. */
    static LineT fromTwoPoints(VectorT const & point1, VectorT const & point2)
    {
      return fromPointAndDirection(point1, point2 - point1);
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
      return (p - closestPoint(p)).squaredLength();
    }

    /** Get the point on the line closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      T t = direction.dot(p - point);
      return point + t * direction;
    }

    /**
     * Get the point on this line closest to a given line. Optionally also returns the point on the other line closest to this
     * line.
     */
    VectorT closestPoint(LineT const & other, VectorT * other_closest_point = NULL) const
    {
      // Adapted from Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.

      VectorT r = point - other.point;
      T b = direction.dot(other.direction);
      T f = other.direction.dot(r);
      T denom = 1 - b * b; // always nonnegative

      // If segments not parallel, compute closest point on L1 to L2. Else pick arbitrary s (here 0)
      T s = 0;
      if (Math::fuzzyGt(denom, static_cast<T>(0)))
      {
        T c = direction.dot(r);
        s = (b * f - c) / denom;
      }

      if (other_closest_point)
      {
        // Compute point on L2 closest to S1(s) using t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f)
        T t = b * s + f;
        *other_closest_point = other.point + t * direction;
      }

      return point + s * direction;
    }

    /** Get a textual description of the line. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[P: " << point.toString() << ", U: " << direction << ']';
      return oss.str();
    }

  private:
    VectorT point;      ///< A point on the line.
    VectorT direction;  ///< A unit vector along the direction of the line.

}; // class LineNBase

} // namespace Internal

/** A straight line in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <long N, typename T>
class /* THEA_API */ LineN : public Internal::LineNBase<N, T>
{
}; // class LineN

} // namespace Thea

#include "Line2.hpp"
#include "Line3.hpp"

#endif
