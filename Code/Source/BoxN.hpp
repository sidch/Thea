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

#ifndef __Thea_BoxN_hpp__
#define __Thea_BoxN_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include "RayIntersectableN.hpp"
#include "VectorN.hpp"

namespace Thea {

/** An arbitrarily oriented box, implemented as an axis-aligned box in a coordinate frame. */
template <long N, typename T = Real>
class /* THEA_API */ BoxN : public RayIntersectableN<N, T>
{
  public:
    typedef VectorN<N, T>           VectorT;           ///< N-dimensional vector.
    typedef AxisAlignedBoxN<N, T>   AxisAlignedBoxT;   ///< N-dimensional axis-aligned box.
    typedef CoordinateFrameN<N, T>  CoordinateFrameT;  ///< N-dimensional coordinate frame.

    THEA_DEF_POINTER_TYPES(BoxN, shared_ptr, weak_ptr)

    /** Default constructor. Creates a null box. */
    BoxN() {}

    /** Initialize to an axis-aligned box in a coordinate frame. */
    BoxN(AxisAlignedBoxT const & aab_, CoordinateFrameT const & frame_ = CoordinateFrameT::identity())
    : aab(aab_), frame(frame_) {}

    /** Get the axis-aligned box in the local frame for the box. */
    AxisAlignedBoxT const & getLocalAAB() const { return aab; }

    /** Get the axis-aligned box in the local frame for the box. */
    AxisAlignedBoxT & getLocalAAB() { return aab; }

    /** Set the axis-aligned box in the local frame for the box. */
    void setLocalAAB(AxisAlignedBoxT const & aab_) { aab = aab_; }

    /** Get the local frame for the box. */
    CoordinateFrameT const & getLocalFrame() const { return frame; }

    /** Get the local frame for the box. */
    CoordinateFrameT & getLocalFrame() { return frame; }

    /** Set the local frame for the box. */
    void setLocalFrame(CoordinateFrameT const & frame_) { frame = frame_; }

    /** Get the center of the box. */
    VectorT getCenter() const { return frame.pointToWorldSpace(aab.getCenter()); }

    /** Get the extent of the box in the local frame. */
    VectorT getExtent() const { return aab.getExtent(); }

    /**
     * Get the i'th corner of the box, for \a i in the range [0, 2^N - 1]. This function works as expected only if
     * N <= sizeof(unsigned long).
     */
    VectorT getCorner(unsigned long i) const
    {
      return frame.pointToWorldSpace(aab.getCorner(i));
    }

    /** Check if the box intersects (i.e. contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the box intersects an axis-aligned box. */
    bool intersects(AxisAlignedBoxT const & aab) const
    {
      throw Error("BoxN: Intersection with axis-aligned box not implemented");
    }

    /** Check if the box intersects another. */
    bool intersects(BoxN const & other) const
    {
      throw Error("BoxN: Box-box intersection not implemented");
    }

    /** Check if the box contains a point. */
    bool contains(VectorT const & p) const
    {
      return aab.contains(frame.pointToObjectSpace(p));
    }

    /** Check if the box contains an axis-aligned box. */
    bool contains(AxisAlignedBoxT const & aab_) const
    {
      // FIXME: Currently works only for N < sizeof(unsigned long)
      for (unsigned long i = 0; i < (1 << N); ++i)
        if (!contains(aab_.getCorner(i)))
          return false;

      return true;
    }

    /** Check if the box contains another. */
    bool contains(BoxN const & other) const
    {
      CoordinateFrameT combined_transform = frame.inverse() * other.frame;

      // FIXME: Currently works only for N < sizeof(unsigned long)
      for (unsigned long i = 0; i < (1 << N); ++i)
        if (!aab.contains(combined_transform * other.aab.getCorner(i)))
          return false;

      return true;
    }

    /** Get an axis-aligned bounding box for the box. */
    AxisAlignedBoxT getBounds() const
    {
      return aab.transformAndBound(frame);
    }

    bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    {
      return aab.rayIntersects(ray.toObjectSpace(frame), max_time);
    }

    T rayIntersectionTime(RayN<N, T> const & ray, T max_time = -1) const
    {
      return aab.rayIntersectionTime(ray.toObjectSpace(frame), max_time);
    }

    RayIntersectionN<N, T> rayIntersection(RayN<N, T> const & ray, T max_time = -1) const
    {
      return aab.rayIntersection(ray.toObjectSpace(frame), max_time);
    }

  private:
    AxisAlignedBoxT aab;     ///< The axis-aligned box in the local frame.
    CoordinateFrameT frame;  ///< The local frame for the box.

}; // class BoxN

} // namespace Thea

#include "Box3.hpp"

#endif
