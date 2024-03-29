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

#ifndef __Thea_LineSegmentN_hpp__
#define __Thea_LineSegmentN_hpp__

#include "Common.hpp"
#include "LineN.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayN.hpp"

namespace Thea {

// Forward declarations
template <int N, typename T> class LineSegmentN;
template <int N, typename T> class AxisAlignedBoxN;

namespace Internal {

// Get the closest pair of points between two line segments, and the square of the distance between them. Adapted from Christer
// Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.
template <int N, typename T> T closestPtSegmentSegment(Vector<N, T> const & p1, Vector<N, T> const & q1, bool is_line1,
                                                       Vector<N, T> const & p2, Vector<N, T> const & q2, bool is_line2,
                                                       T & s, T & t, Vector<N, T> & c1, Vector<N, T> & c2);

/**
 * <b>[Internal]</b> Base class for straight line segments in N-dimensional space, where N is any <b>positive</b> (non-zero)
 * integer and T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ LineSegmentNBase
{
  public:
    typedef LineSegmentN<N, T>  LineSegmentT;  ///< N-dimensional straight line.
    typedef Vector<N, T>        VectorT;       ///< N-dimensional vector.

    THEA_DECL_SMART_POINTERS(LineSegmentT)

    /** Default constructor, does not initialize the segment. */
    LineSegmentNBase() {}

    /** Construct the line segment from its endpoints. */
    LineSegmentNBase(VectorT const & point1, VectorT const & point2)
    : point(point1), direction(point2 - point1)
    {}

    /**
     * Construct the line segment from its first endpoint and its direction vector, with length equal to that of the segment.
     */
    static LineSegmentT fromPointAndDirection(VectorT const & point_, VectorT const & direction_)
    {
      LineSegmentT seg;
      seg.point = point_;
      seg.direction = direction_;
      return seg;
    }

    /** Cast the line segment to a different scalar type. */
    template <typename U> LineSegmentN<N, U> cast() const
    {
      return LineSegmentN<N, U>::fromPointAndDirection(point.template cast<U>(), direction.template cast<U>());
    }

    /** Get an endpoint of the line segment: 0 returns the first endpoint and 1 returns the second. */
    VectorT getEndpoint(int i) const { return i == 0 ? point : point + direction; }

    /** Get the unnormalized direction vector of the segment from the first endpoint to the second. */
    VectorT const & getDirection() const { return direction; }

    /** Get a point on the line segment: \a t = 0 maps to the first endpoint and \a t = 1 maps to the second. */
    VectorT getPoint(Real t) const { return point + t * direction; }

    /** Get the length of the line segment. */
    T length() const { return direction.norm(); }

    /** Get the square of the length of the line segment. */
    T squaredLength() const { return direction.squaredNorm(); }

    /** Get the distance of the line segment from a given point. */
    T distance(VectorT const & p) const
    {
      return std::sqrt(squaredDistance(p));
    }

    /** Get the square of the distance of the line segment from a given point. */
    T squaredDistance(VectorT const & p) const
    {
      return (p - closestPoint(p)).squaredNorm();
    }

    /** Get the point on the line segment closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      // Taken from G3D::LineSegment

      T d2 = direction.squaredNorm();
      if (Math::fuzzyEq(d2, static_cast<T>(0)))
        return point;

      // The vector from the end of the segment to the point in question.
      VectorT v(p - point);

      // Projection of v onto the line segment scaled by the length of the segment.
      T t = direction.dot(v);

      // Avoid some square roots. Derivation:
      //    t / direction.norm() <= direction.norm()
      //    t <= direction.squaredNorm()

      if (t >= 0 && t <= d2)
      {
        // The point falls within the segment. Normalize direction, divide t by the length of direction.
        return point + (t / d2) * direction;
      }
      else
      {
        // The point does not fall within the segment; see which end is closer.

        // Distance from 0, squared
        T d0_squared = v.squaredNorm();

        // Distance from 1, squared
        T d1_squared = (v - direction).squaredNorm();

        if (d0_squared < d1_squared)  // point 0 is closer
          return point;
        else  // point 1 is closer
          return point + direction;
      }
    }

    /** Get the distance of this segment from another segment. */
    T distance(LineSegmentT const & other) const
    {
      return std::sqrt(squaredDistance(other));
    }

    /** Get the squared distance between this segment and another segment, and optionally return the closest pair of points. */
    T squaredDistance(LineSegmentT const & other, VectorT * this_pt = nullptr, VectorT * other_pt = nullptr) const
    {
      VectorT c1, c2;
      T s, t;
      auto sqdist = Internal::closestPtSegmentSegment<N, T>(point, VectorT(point + direction), false, other.point,
                                                            VectorT(other.point + other.direction), false, s, t, c1, c2);

      if (this_pt)  *this_pt  = c1;
      if (other_pt) *other_pt = c2;

      return sqdist;
    }

    /** Get the distance of this segment from an infinite line. */
    T distance(LineN<N, T> const & line) const
    {
      return std::sqrt(squaredDistance(line));
    }

    /** Get the squared distance between this segment and an infinite line, and optionally return the closest pair of points. */
    T squaredDistance(LineN<N, T> const & line, VectorT * this_pt = nullptr, VectorT * line_pt = nullptr) const
    {
      VectorT c1, c2;
      T s, t;
      auto sqdist = Internal::closestPtSegmentSegment<N, T>(point, VectorT(point + direction), false, line.getPoint(),
                                                            VectorT(line.getPoint() + line.getDirection()), true, s, t, c1, c2);

      if (this_pt) *this_pt = c1;
      if (line_pt) *line_pt = c2;

      return sqdist;
    }

    /** Get the distance of this segment from a ray. */
    T distance(RayN<N, T> const & ray) const
    {
      return std::sqrt(squaredDistance(ray));
    }

    /** Get the squared distance between this segment and a ray, and optionally return the closest pair of points. */
    T squaredDistance(RayN<N, T> const & ray, VectorT * this_pt = nullptr, VectorT * ray_pt = nullptr) const
    {
      VectorT c1, c2;
      T s, t;
      auto sqdist = Internal::closestPtSegmentSegment<N, T>(point, VectorT(point + direction), false, ray.getOrigin(),
                                                            VectorT(ray.getOrigin() + ray.getDirection()), true, s, t, c1, c2);
      if (t < 0)
        c2 = ray.getOrigin();

      if (this_pt) *this_pt = c1;
      if (ray_pt)  *ray_pt  = c2;

      return sqdist;
    }

    /**
     * Get the parameter <tt>t</tt> of a point \a p, assumed to be on the line segment, such that <tt>p = getPoint(t)</tt>.
     *
     * @see getPoint()
     */
    T parametrize(VectorT const & p) const
    {
      return (p - point).dot(direction) / direction.squaredNorm();
    }

    /** Get a bounding box for the line segment. */
    AxisAlignedBoxN<N, T> getBounds() const
    {
      AxisAlignedBoxN<N, T> box(point);
      box.merge(point + direction);
      return box;
    }

    /** Get a textual representation of the line segment. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << '[' << Thea::toString(getEndpoint(0)) << ", " << Thea::toString(getEndpoint(1)) << ']';
      return oss.str();
    }

  private:
    VectorT point;      ///< The first endpoint of the segment.
    VectorT direction;  ///< The unnormalized vector from the first endpoint to the second.

}; // class LineSegmentNBase

} // namespace Internal

/** A straight line in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <int N, typename T = Real>
class /* THEA_API */ LineSegmentN : public Internal::LineSegmentNBase<N, T>
{
  private:
    typedef Internal::LineSegmentNBase<N, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    /** Default constructor, does not initialize the segment. */
    LineSegmentN() {}

    /** Construct the line segment from its endpoints. */
    LineSegmentN(VectorT const & point1, VectorT const & point2) : BaseT(point1, point2) {}

}; // class LineSegmentN

namespace Internal {

template <int N, typename T>
T
closestPtSegmentSegment(Vector<N, T> const & p1, Vector<N, T> const & q1, bool is_line1,
                        Vector<N, T> const & p2, Vector<N, T> const & q2, bool is_line2,
                        T & s, T & t, Vector<N, T> & c1, Vector<N, T> & c2)
{
  typedef Vector<N, T> VectorT;

  VectorT d1 = q1 - p1;  // Direction vector of segment S1
  VectorT d2 = q2 - p2;  // Direction vector of segment S2
  VectorT r = p1 - p2;
  T a = d1.squaredNorm();  // Squared length of segment S1, always nonnegative
  T e = d2.squaredNorm();  // Squared length of segment S2, always nonnegative
  T f = d2.dot(r);

  // Check if either or both segments degenerate into points
  if (Math::fuzzyEq(a, static_cast<T>(0)))
  {
    if (Math::fuzzyEq(e, static_cast<T>(0)))
    {
      // Both segments degenerate into points
      s = t = 0;
      c1 = p1;
      c2 = p2;
      return (c1 - c2).dot(c1 - c2);
    }
    else
    {
      // First segment degenerates into a point
      s = 0;
      t = f / e;  // s = 0 => t = (b*s + f) / e = f / e

      if (!is_line2)
        t = Math::clamp(t, static_cast<T>(0), static_cast<T>(1));
    }
  }
  else
  {
    T c = d1.dot(r);
    if (Math::fuzzyEq(e, static_cast<T>(0)))
    {
      // Second segment degenerates into a point
      t = 0;
      s = -c / a;

      if (!is_line1)
        s = Math::clamp(s, static_cast<T>(0), static_cast<T>(1));  // t = 0 => s = (b*t - c) / a = -c / a
    }
    else
    {
      // The general nondegenerate case starts here
      T b = d1.dot(d2);
      T denom = a * e - b * b; // Always nonnegative

      // If segments not parallel, compute closest point on L1 to L2, and clamp to segment S1. Else pick arbitrary s (here 0)
      if (denom != 0)
      {
        s = (b * f - c * e) / denom;

        if (!is_line1)
          s = Math::clamp(s, static_cast<T>(0), static_cast<T>(1));
      }
      else
        s = 0;

      // Compute point on L2 closest to S1(s) using t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
      t = (b * s + f) / e;

      if (!is_line2)
      {
        // If t in [0,1] done. Else clamp t, recompute s for the new value of t using
        // s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a and clamp s to [0, 1]
        if (t < 0)
        {
          t = 0;
          s = -c / a;

          if (!is_line1)
            s = Math::clamp(s, static_cast<T>(0), static_cast<T>(1));
        }
        else if (t > 1)
        {
          t = 1;
          s = (b - c) / a;

          if (!is_line1)
            s = Math::clamp(s, static_cast<T>(0), static_cast<T>(1));
        }
      }
    }
  }

  c1 = p1 + s * d1;
  c2 = p2 + t * d2;
  return (c1 - c2).squaredNorm();
}

} // namespace Internal

} // namespace Thea

#include "LineSegment2.hpp"
#include "LineSegment3.hpp"

#endif
