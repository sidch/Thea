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

#ifndef __Thea_BoxN_hpp__
#define __Thea_BoxN_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"

namespace Thea {

/** An arbitrarily oriented box, implemented as an axis-aligned box in a coordinate frame. */
template <int N, typename T = Real>
class /* THEA_API */ BoxN : public RayIntersectableN<N, T>
{
  public:
    typedef Vector<N, T>            VectorT;           ///< N-dimensional vector.
    typedef AxisAlignedBoxN<N, T>   AxisAlignedBoxT;   ///< N-dimensional axis-aligned box.
    typedef CoordinateFrameN<N, T>  CoordinateFrameT;  ///< N-dimensional coordinate frame.

    THEA_DECL_SMART_POINTERS(BoxN)

    /** Default constructor. Creates a null box. */
    BoxN() {}

    /** Initialize to an axis-aligned box in a coordinate frame. */
    BoxN(AxisAlignedBoxT const & aab_, CoordinateFrameT const & frame_ = CoordinateFrameT::identity())
    : aab(aab_), frame(frame_) {}

    /** Cast the box to a different scalar type. */
    template <typename U> BoxN<N, U> cast() const { return BoxN<N, U>(aab.template cast<U>(), frame.template cast<U>()); }

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
     * N <= sizeof(uintx).
     */
    VectorT getCorner(uintx i) const
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
      // FIXME: Currently works only for N < sizeof(uintx)
      for (uintx i = 0; i < (1 << N); ++i)
        if (!contains(aab_.getCorner(i)))
          return false;

      return true;
    }

    /** Check if the box contains another. */
    bool contains(BoxN const & other) const
    {
      CoordinateFrameT combined_transform = frame.inverse() * other.frame;

      // FIXME: Currently works only for N < sizeof(uintx)
      for (uintx i = 0; i < (1 << N); ++i)
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
