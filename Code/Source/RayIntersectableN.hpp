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

#ifndef __Thea_RayIntersectableN_hpp__
#define __Thea_RayIntersectableN_hpp__

#include "Common.hpp"
#include "RayN.hpp"

namespace Thea {

/**
 * A description of the intersection point of a ray with an object. Specifies the hit time and the normal at the intersection
 * point.
 */
template <long N, typename T = Real>
class /* THEA_API */ RayIntersectionN
{
  public:
    typedef VectorN<N, T> VectorT;  ///< N-dimensional vector.

  private:
    T time;
    bool has_normal;
    VectorT normal;

  public:
    /** Constructor. */
    RayIntersectionN(T time_ = -1, VectorT const * normal_ = NULL)
    : time(time_), has_normal(normal_ != NULL), normal(normal_ ? *normal_ : VectorT::zero())
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

/** Interface for an object that supports ray intersection queries in N-space. */
template <long N, typename T = Real>
class /* THEA_API */ RayIntersectableN
{
  public:
    THEA_DEF_POINTER_TYPES(RayIntersectableN, shared_ptr, weak_ptr)

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
