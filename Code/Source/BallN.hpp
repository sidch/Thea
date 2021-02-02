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

#ifndef __Thea_BallN_hpp__
#define __Thea_BallN_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"
#include <sstream>

namespace Thea {

/**
 * A ball in N-space.
 *
 * @note A (N - 1)-sphere is the (N - 1)-dimensional surface of a N-dimensional ball.
 */
template <int N, typename T = Real>
class /* THEA_API */ BallN : public RayIntersectableN<N, T>
{
  public:
    typedef Vector<N, T> VectorT;

    THEA_DECL_SMART_POINTERS(BallN)

    /** Default constructor. Does not initialize anything. */
    BallN() {}

    /** Initialize with a center and a radius. */
    BallN(VectorT const & center_, T const & radius_) : center(center_), radius(radius_) {}

    /** Copy constructor. */
    BallN(BallN const & src) : center(src.center), radius(src.radius) {}

    /** Cast the ball to a different scalar type. */
    template <typename U> BallN<N, U> cast() const { return BallN<N, U>(center.template cast<U>(), static_cast<U>(radius)); }

    /** Get the center of the ball. */
    VectorT const & getCenter() const { return center; }

    /** Set the center of the ball. */
    void setCenter(VectorT const & center_) { center = center_; }

    /** Get the radius of the ball. */
    T getRadius() const { return radius; }

    /** Get the diameter of the ball. */
    T getDiameter() const { return 2 * radius; }

    /** Set the radius of the ball. */
    void setRadius(T radius_) { radius = radius_; }

    /** Test if this ball intersects (contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the ball intersects another ball. */
    bool intersects(BallN const & other) const
    {
      return (center - other.center).norm() < radius + other.radius;
    }

    /** Check if the ball intersects an axis-aligned box. */
    bool intersects(AxisAlignedBoxN<N, T> const & aab) const
    {
      return aab.squaredDistance(center) <= radius * radius;
    }

    /** Check if the ball contains a point. */
    bool contains(VectorT const & p) const { return (p - center).squaredNorm() <= radius * radius; }

    /** Check if the ball contains another ball. */
    bool contains(BallN const & other) const
    {
      return radius >= other.radius && (center - other.center).norm() < radius - other.radius;
    }

    /** Check if the ball contains an axis-aligned box. */
    bool contains(AxisAlignedBoxN<N, T> const & aab) const
    {
      // FIXME: Currently works only for N < sizeof(uintx)
      for (uintx i = 0; i < (1 << N); ++i)
        if (!contains(aab.getCorner(i)))
          return false;

      return true;
    }

    /** Get the distance of the ball from a point. */
    T distance(VectorT const & p) const
    {
      return std::max((p - center).norm() - radius, static_cast<T>(0));
    }

    /** Get the distance of the ball from a point. */
    T distance(BallN const & other) const
    {
      return std::max((other.center - center).norm() - radius - other.radius, static_cast<T>(0));
    }

    /** Get the distance of the ball from an axis-aligned box. */
    T distance(AxisAlignedBoxN<N, T> const & aab) const
    {
      return std::max(aab.distance(center) - radius, static_cast<T>(0));
    }

    /** Get the squared distance of the ball from a point, ball or axis-aligned box. */
    template <typename OtherType> T squaredDistance(OtherType const & other) const
    { T x = distance(other); return x * x; }

    /** Get the squared distance of the ball from an axis-aligned box and optionally return the closest pair of points. */
    T squaredDistance(AxisAlignedBoxN<N, T> const & aab, VectorT * this_pt = nullptr, VectorT * other_pt = nullptr) const
    {
      if (this_pt || other_pt)
        throw FatalError("BallN: Finding closest points between ball and box not yet implemented");

      T x = distance(aab); return x * x;
    }

    /** Get a bounding box for the ball. */
    AxisAlignedBoxN<N, T> getBounds() const
    {
      VectorT half_ext(radius);
      return AxisAlignedBoxN<N, T>(center - half_ext, center + half_ext);
    }

    /** Get a textual representation of the ball. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[center: " << Thea::toString(center) << ", radius: " << radius << ']';
      return oss.str();
    }

    bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    {
      if (max_time >= 0)
        return rayIntersectionTime(ray, max_time) >= 0;

      VectorT co = ray.getOrigin() - center;

      T c = co.squaredNorm() - radius * radius;
      if (c <= 0)  // origin is inside ball
        return true;

      T a = ray.getDirection().squaredNorm();
      T b = 2 * co.dot(ray.getDirection());

      // Solve quadratic a * t^2 + b * t + c = 0
      T b2 = b * b;
      T det = b2 - 4 * a * c;
      if (det < 0) return false;

      if (a > 0)
        return b <= 0 || det >= b2;
      else if (a < 0)
        return b >= 0 || det >= b2;
      else
        return false;
    }

    T rayIntersectionTime(RayN<N, T> const & ray, T max_time = -1) const
    {
      VectorT co = ray.getOrigin() - center;

      T c = co.squaredNorm() - radius * radius;
      if (c <= 0)  // origin is inside ball
        return 0;

      // We could do an early test to see if the distance from the ray origin to the ball is less than
      // max_time * ray.getDirection().norm(), but it would involve a square root so might as well solve the quadratic.

      T a = ray.getDirection().squaredNorm();
      T b = 2 * co.dot(ray.getDirection());

      // Solve quadratic a * t^2 + b * t + c = 0
      T det = b * b - 4 * a * c;
      if (det < 0) return -1;

      T d = std::sqrt(det);
      T t = -1;
      if (a > 0)
      {
        T s0 = -b - d;
        if (s0 >= 0)
          t = s0 / (2 * a);
        else
        {
          T s1 = -b + d;
          if (s1 >= 0)
            t = s1 / (2 * a);
        }
      }
      else if (a < 0)
      {
        T s0 = -b + d;
        if (s0 <= 0)
          t = s0 / (2 * a);
        else
        {
          T s1 = -b - d;
          if (s1 <= 0)
            t = s1 / (2 * a);
        }
      }

      if (max_time >= 0 && t > max_time)
        return -1;
      else
        return t;
    }

    RayIntersectionN<N, T> rayIntersection(RayN<N, T> const & ray, T max_time = -1) const
    {
      T t = rayIntersectionTime(ray, max_time);
      if (t >= 0)
      {
        VectorT p = ray.getPoint(t);
        VectorT n = p - center;
        return RayIntersectionN<N, T>(t, &n);
      }
      else
        return RayIntersectionN<N, T>(-1);
    }

  private:
    VectorT  center;  ///< Center of the ball.
    T        radius;  ///< Radius of the ball.

}; // class BallN

} // namespace Thea

#include "Ball3.hpp"

#endif
