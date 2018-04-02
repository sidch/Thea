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

#ifndef __Thea_BallN_hpp__
#define __Thea_BallN_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include "Math.hpp"
#include "RayIntersectableN.hpp"
#include "VectorN.hpp"
#include <sstream>

namespace Thea {

/**
 * A ball in N-space.
 *
 * @note A (N - 1)-sphere is the (N - 1)-dimensional surface of a N-dimensional ball.
 */
template <long N, typename T = Real>
class /* THEA_API */ BallN : public RayIntersectableN<N, T>
{
  public:
    typedef VectorN<N, T> VectorT;

    THEA_DEF_POINTER_TYPES(BallN, shared_ptr, weak_ptr)

    /** Default constructor. Does not initialize anything. */
    BallN() {}

    /** Initialize with a center and a radius. */
    BallN(VectorT const & center_, T const & radius_) : center(center_), radius(radius_) {}

    /** Copy constructor. */
    BallN(BallN const & src) : center(src.center), radius(src.radius) {}

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
      return (center - other.center).length() < radius + other.radius;
    }

    /** Check if the ball intersects an axis-aligned box. */
    bool intersects(AxisAlignedBoxN<N, T> const & aab) const
    {
      return aab.squaredDistance(center) <= radius * radius;
    }

    /** Check if the ball contains a point. */
    bool contains(VectorT const & p) const { return (p - center).squaredLength() <= radius * radius; }

    /** Check if the ball contains another ball. */
    bool contains(BallN const & other) const
    {
      return radius >= other.radius && (center - other.center).length() < radius - other.radius;
    }

    /** Check if the ball contains an axis-aligned box. */
    bool contains(AxisAlignedBoxN<N, T> const & aab) const
    {
      // FIXME: Currently works only for N < sizeof(unsigned long)
      for (unsigned long i = 0; i < (1 << N); ++i)
        if (!contains(aab.getCorner(i)))
          return false;

      return true;
    }

    /** Get the distance of the ball from a point. */
    T distance(VectorT const & p) const
    {
      return std::max((p - center).length() - radius, static_cast<T>(0));
    }

    /** Get the distance of the ball from a point. */
    T distance(BallN const & other) const
    {
      return std::max((other.center - center).length() - radius - other.radius, static_cast<T>(0));
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
    T squaredDistance(AxisAlignedBoxN<N, T> const & aab, VectorT * this_pt = NULL, VectorT * other_pt = NULL) const
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
      oss << "[center: " << center.toString() << ", radius: " << radius << ']';
      return oss.str();
    }

    bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    {
      if (max_time >= 0)
        return rayIntersectionTime(ray, max_time) >= 0;

      VectorT co = ray.getOrigin() - center;

      T c = co.squaredLength() - radius * radius;
      if (c <= 0)  // origin is inside ball
        return true;

      T a = ray.getDirection().squaredLength();
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

      T c = co.squaredLength() - radius * radius;
      if (c <= 0)  // origin is inside ball
        return 0;

      // We could do an early test to see if the distance from the ray origin to the ball is less than
      // max_time * ray.getDirection().length(), but it would involve a square root so might as well solve the quadratic.

      T a = ray.getDirection().squaredLength();
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
