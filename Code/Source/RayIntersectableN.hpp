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

#ifndef __Thea_RayIntersectableN_hpp__
#define __Thea_RayIntersectableN_hpp__

#include "Common.hpp"
#include "RayN.hpp"

namespace Thea {

/**
 * A description of the intersection point of a ray with an object. Specifies the hit time and the normal at the intersection
 * point.
 */
template <int N, typename T = Real>
class /* THEA_API */ RayIntersectionN
{
  public:
    typedef Vector<N, T> VectorT;  ///< N-dimensional vector.

  private:
    T time;
    bool has_normal;
    VectorT normal;

  public:
    /** Constructor. */
    RayIntersectionN(T time_ = -1, VectorT const * normal_ = nullptr)
    : time(time_), has_normal(normal_ != nullptr), normal(normal_ ? *normal_ : VectorT::Zero())
    {}

    /** Check if the intersection is valid. */
    bool isValid() const { return time >= 0; }

    /** Get the intersection time. */
    T getTime() const { return time; }

    /** Set the intersection time. */
    void setTime(T time_) { time = time_; }

    /** Check if the normal at the intersection point is known. */
    bool hasNormal() const { return has_normal; }

    /** Get the normal at the intersection point. The return value is undefined if hasNormal() returns false. */
    VectorT const & getNormal() const { return normal; }

    /** Set the normal at the intersection point. hasNormal() will subsequently return true. */
    void setNormal(VectorT const & normal_) { normal = normal_; has_normal = true; }

}; // class RayIntersectionN

/** Abstract base class for an object that supports ray intersection queries in N-space. */
template <int N, typename T = Real>
class /* THEA_API */ RayIntersectableN
{
  public:
    THEA_DECL_SMART_POINTERS(RayIntersectableN)

    /** Destructor. */
    virtual ~RayIntersectableN() {}

    /**
     * Check if a ray intersects the object in the forward direction.
     *
     * @param ray The ray to test for intersection.
     * @param max_time Maximum allowable hit time, ignored if negative.
     */
    virtual bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    { return rayIntersectionTime(ray, max_time) >= 0; }

    /**
     * Get the time taken for a ray to intersect the object, or a negative value if there was no intersection in the forward
     * direction. All subclasses must reimplement this method. If the return value is negative, it should be at least
     *
     * @param ray The ray to test for intersection.
     * @param max_time Maximum allowable hit time, ignored if negative.
     */
    virtual T rayIntersectionTime(RayN<N, T> const & ray, T max_time = -1) const = 0;

    /**
     * Get the intersection of a ray with the object, including the hit time and the normal at the intersection point. A
     * negative time is returned if there was no intersection in the forward direction. If the normal cannot be computed, the
     * zero vector is returned.
     *
     * @param ray The ray to test for intersection.
     * @param max_time Maximum allowable hit time, ignored if negative.
     *
     * @note The returned normal need not have unit length.
     */
    virtual RayIntersectionN<N, T> rayIntersection(RayN<N, T> const & ray, T max_time = -1) const
    { return RayIntersectionN<N, T>(rayIntersectionTime(ray, max_time)); }

}; // class RayIntersectableN

} // namespace Thea

#include "RayIntersectable3.hpp"

#endif
