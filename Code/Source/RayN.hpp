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

#ifndef __Thea_RayN_hpp__
#define __Thea_RayN_hpp__

#include "Common.hpp"
#include "CoordinateFrameN.hpp"
#include "MatVec.hpp"

namespace Thea {

/** A ray in N-dimensional space, having an originating point and a direction vector (not necessarily unit length). */
template <int N, typename T = Real>
class /* THEA_API */ RayN
{
  public:
    THEA_DECL_SMART_POINTERS(RayN)

    typedef Vector<N, T> VectorT;  ///< N-dimensional vector.

    /** Default constructor. Does not initialize anything. */
    RayN() {}

    /** Initialize with an origin and a direction. The direction will <b>NOT</b> be normalized. */
    RayN(VectorT const & origin_, VectorT direction_) : origin(origin_), direction(direction_) {}

    /** Copy constructor. */
    RayN(RayN const & src) : origin(src.origin), direction(src.direction) {}

    /** Cast the ray to a different scalar type. */
    template <typename U> RayN<N, U> cast() const
    {
      return RayN<N, U>(origin.template cast<U>(), direction.template cast<U>());
    }

    /** Get the origin of the ray. */
    VectorT const & getOrigin() const { return origin; }

    /** Set the origin of the ray. */
    void setOrigin(VectorT const & origin_) { origin = origin_; }

    /** Get the direction of the ray. */
    VectorT const & getDirection() const { return direction; }

    /** Set the direction of the ray. The direction will <b>NOT</b> be normalized. */
    void setDirection(VectorT const & direction_) { direction = direction_; }

    /** Make the direction vector unit length. */
    void normalizeDirection() { direction.normalize(); }

    /** Transform the ray out of a local frame to world space. */
    RayN toWorldSpace(CoordinateFrameN<N, T> const & frame) const
    {
      return RayN(frame.pointToWorldSpace(origin), frame.vectorToWorldSpace(direction));
    }

    /** Transform the ray from world space into a local frame. */
    RayN toObjectSpace(CoordinateFrameN<N, T> const & frame) const
    {
      return RayN(frame.pointToObjectSpace(origin), frame.vectorToObjectSpace(direction));
    }

    /**
     * Get a parametrized point on the ray.
     *
     * @return getOrigin() + \a t * getDirection()
     */
    VectorT getPoint(T const & t) const { return origin + t * direction; }

    /** Get the distance of a point from the ray. */
    T distance(VectorT const & p) const { return (closestPoint(p) - p).norm(); }

    /** Get the square of the distance of a point from the ray. */
    T squaredDistance(VectorT const & p) const { return (closestPoint(p) - p).squaredNorm(); }

    /** Get the the point on the ray closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      T t = (p - origin).dot(direction);  // we'll normalize direction later
      if (t < 0)
        return origin;
      else
      {
        T dir_sqlen = direction.squaredNorm();
        if (Math::fuzzyEq(dir_sqlen, static_cast<T>(0)))
          return origin;

        // Direction has to be normalized twice -- once during the dot product above and once in the scaling below. We can
        // combine the two divisions and avoid the square roots.
        return origin + (t / dir_sqlen) * direction;
      }
    }

    /** Get a textual representation of the ray. */
    std::string toString() const
    {
      return "[origin: " + Thea::toString(origin) + ", direction: " + Thea::toString(direction) + ']';
    }

  private:
    VectorT origin;     ///< Origin of the ray.
    VectorT direction;  ///< Direction of the ray.

}; // class RayN

} // namespace Thea

#include "Ray3.hpp"

#endif
