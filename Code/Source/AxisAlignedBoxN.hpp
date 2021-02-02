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
// First version: 2009
//
//============================================================================

#ifndef __Thea_AxisAlignedBoxN_hpp__
#define __Thea_AxisAlignedBoxN_hpp__

#include "Common.hpp"
#include "Math.hpp"
#include "LineN.hpp"
#include "LineSegmentN.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"
#include "RayN.hpp"
#include <bitset>

namespace Thea {

// Forward declarations
template <int N, typename T> class AxisAlignedBoxN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for axis-aligned boxes in N-dimensional space, where N is any <b>positive</b> (non-zero) integer
 * and T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ AxisAlignedBoxNBase : public RayIntersectableN<N, T>
{
  public:
    typedef AxisAlignedBoxN<N, T>  AxisAlignedBoxT;  ///< N-dimensional axis-aligned box.
    typedef Vector<N, T>           VectorT;          ///< N-dimensional vector.

    THEA_DECL_SMART_POINTERS(AxisAlignedBoxT)

    /**
     * Default constructor, creates a null box.
     *
     * @see isNull()
     */
    AxisAlignedBoxNBase() : lo(VectorT::Zero()), hi(VectorT::Zero()), is_null(true) {}

    /** Constructor. Sets the box to be a single point. */
    AxisAlignedBoxNBase(VectorT const & v) : lo(v), hi(v), is_null(false) {}

    /** Constructor. Sets the extents of the box. */
    AxisAlignedBoxNBase(VectorT const & lo_, VectorT const & hi_) : lo(lo_), hi(hi_), is_null(false) {}

    /** Cast the box to a different scalar type. */
    template <typename U> AxisAlignedBoxN<N, U> cast() const
    {
      return AxisAlignedBoxN<N, U>(lo.template cast<U>(), hi.template cast<U>());
    }

    /** Check if two boxes are equal. Null boxes are considered to be equal to each other. */
    bool operator==(AxisAlignedBoxNBase const & other) const
    {
      return (is_null && other.is_null) || (!is_null && !other.is_null && lo == other.lo && hi == other.hi);
    }

    /**
     * Check if the box is null. A null box has no extent or location in space. It is <i>not</i> equivalent to a box with zero
     * extent. All null boxes are considered equal.
     */
    bool isNull() const { return is_null; }

    /**
     * Set the box to be null.
     *
     * @see isNull()
     */
    void setNull()
    {
      lo = hi = VectorT::Zero();
      is_null = true;
    }

    /** Get the minimum corner of the box. The return value is zero if the box is null. */
    VectorT const & getLow() const { return lo; }

    /** Get the maximum corner of the box. The return value is zero if the box is null. */
    VectorT const & getHigh() const { return hi; }

    /** Set the corners of the box. */
    void set(VectorT const & lo_, VectorT const & hi_)
    {
      lo = lo_;
      hi = hi_;
      is_null = false;
    }

    /**
     * Get the i'th corner of the box, for \a i in the range [0, 2^N - 1]. This function works as expected only if
     * N <= sizeof(uintx). The d'th coordinate of the returned point is assigned the d'th coordinate of the maximum corner of
     * the box if the d'th bit of i (where the 0'th bit is the least significant) is 1, and the minimum corner if it is 0.
     */
    VectorT getCorner(uintx i) const
    {
      VectorT ret;
      for (intx j = 0; j < N; ++j)
        ret[j] = ((i >> j) & 0x01) ? hi[j] : lo[j];

      return ret;
    }

    /**
     * Transform the box and return a new axis-aligned box which tightly encloses the result.
     *
     * Algorithm derived from the Ogre source code, http://www.ogre3d.org.
     */
    template <typename TransformT> AxisAlignedBoxT transformAndBound(TransformT const & tr) const
    {
      AxisAlignedBoxT ret;     // null box
      std::bitset<N> counter;  // all-zero bitset
      VectorT corner;
      do
      {
        // Get the current corner
        for (intx j = 0; j < N; ++j)
          corner[j] = counter.test((size_t)j) ? hi[j] : lo[j];

        // Transform and merge it
        ret.merge(tr * corner);

        // Get the index of the next corner by adding 1 to the current index
        for (size_t j = 0; j < N; j++)
        {
          if (counter.test(j))
            counter.set(j, false);
          else
          {
            counter.set(j, true);
            break;
          }
        }

      } while (counter.any());

      return ret;
    }

    /** Grow to include a point. Equivalent to assignment if this box is initially null. */
    void merge(VectorT const & p)
    {
      if (is_null)
      {
        lo = hi = p;
        is_null = false;
      }
      else
      {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
      }
    }

    /** Grow to include a point. Alias for merge(VectorT const &), for use with PointCollectorN. */
    void addPoint(VectorT const & p) { merge(p); }

    /** Grow to include a box. Equivalent to assignment if this box is initially null. No-op if the argument is a null box. */
    void merge(AxisAlignedBoxT const & aab)
    {
      if (!aab.is_null)
      {
        if (is_null)
        {
          lo = aab.lo;
          hi = aab.hi;
          is_null = false;
        }
        else
        {
          lo = lo.cwiseMin(aab.lo);
          hi = hi.cwiseMax(aab.hi);
        }
      }
    }

    /** Get the center of the box. The return value is zero for null boxes. */
    VectorT getCenter() const { return static_cast<T>(0.5) * (hi + lo); }

    /** Set the center of the box, without changing its size. No-op if the box is null. */
    void setCenter(VectorT const & center)
    {
      if (!is_null)
      {
        VectorT shift = center - getCenter();
        hi += shift;
        lo += shift;
      }
    }

    /** Get the extent of the box. The return value is zero for null boxes. */
    VectorT getExtent() const { return (hi - lo); }

    /** Set the extent (size) of the box, without changing its center. No-op if the box is null. */
    void setExtent(VectorT const & ext)
    {
      if (!is_null)
      {
        VectorT half_ext = static_cast<T>(0.5) * ext;
        VectorT center = getCenter();
        lo = center - half_ext;
        hi = center + half_ext;
      }
    }

    /** Get the volume of the box. A null box has zero volume. */
    T volume() const
    {
      if (is_null) return 0;

      VectorT ext = getExtent();
      T vol = ext[0];
      for (intx i = 1; i < N; ++i)
        vol *= ext[i];

      return vol;
    }

    /** Scale the box by a linear factor relative to its center. No-op if the box is null. */
    void scaleCentered(T const & scale)
    {
      if (!is_null)
      {
        VectorT c = getCenter();
        lo = scale * (lo - c) + c;
        hi = scale * (hi - c) + c;
      }
    }

    /** Scale the box relative to its center. No-op if the box is null. */
    void scaleCentered(VectorT const & scale)
    {
      if (!is_null)
      {
        VectorT c = getCenter();
        lo = scale * (lo - c) + c;
        hi = scale * (hi - c) + c;
      }
    }

    /** Return a box that is a copy of this, scaled by a linear factor relative to its center. */
    AxisAlignedBoxT scaleCenteredCopy(T const & scale) const
    {
      if (is_null) return AxisAlignedBoxT();

      VectorT c = getCenter();
      return AxisAlignedBoxT(scale * (lo - c) + c, scale * (hi - c) + c);
    }

    /** Return a box that is a copy of this, scaled relative to its center. */
    AxisAlignedBoxT scaleCenteredCopy(VectorT const & scale) const
    {
      if (is_null) return AxisAlignedBoxT();

      VectorT c = getCenter();
      return AxisAlignedBoxT(scale * (lo - c) + c, scale * (hi - c) + c);
    }

    /** Test if this box intersects (i.e. contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Test if this box intersects another. Null boxes do not intersect each other. */
    bool intersects(AxisAlignedBoxT const & other) const
    {
      if (is_null || other.is_null) return false;

      // Look for the absence of a separating plane
      for (intx i = 0; i < N; ++i)
        if (lo[i] > other.hi[i] || hi[i] < other.lo[i])
          return false;

      return true;
    }

    /** Test if this box contains a point. */
    bool contains(VectorT const & p) const
    {
      if (is_null) return false;

      for (intx i = 0; i < N; ++i)
      {
        if (p[i] < lo[i]) return false;
        if (p[i] > hi[i]) return false;
      }

      return true;
    }

    /** Test if this box contains another. Null boxes do not contain or are contained by other boxes. */
    bool contains(AxisAlignedBoxT const & other) const
    {
      if (is_null || other.is_null) return false;

      for (intx i = 0; i < N; ++i)
      {
        if (other.lo[i] < lo[i]) return false;
        if (other.hi[i] > hi[i]) return false;
      }

      return true;
    }

    /** Get the square of the closest distance between this box and another one. Returns zero if either box is null. */
    T squaredDistance(AxisAlignedBoxT const & other) const
    {
      if (is_null || other.is_null) return 0;

      VectorT vmax = hi.cwiseMin(other.hi);
      VectorT vmin = lo.cwiseMax(other.lo);

      // Each coord of vmax is greater than the corresponding coord of vmin iff the ranges intersect

      VectorT vdiff = vmin - vmax;         // any axes with overlap have negative values here
      vdiff = vdiff.cwiseMax(VectorT::Zero());  // overlap axes have zero separation

      return vdiff.squaredNorm();
    }

    /** Get the closest distance between this box and another one. Returns zero if either box is null. */
    T distance(AxisAlignedBoxT const & other) const { return static_cast<T>(std::sqrt(squaredDistance(other))); }

    /** Get the square of the closest distance to a point. Returns zero if the box is null. */
    T squaredDistance(VectorT const & point) const
    {
      if (is_null) return 0;

      VectorT vmax = point.cwiseMin(hi);
      VectorT vmin = point.cwiseMax(lo);

      return (vmin - vmax).squaredNorm();
    }

    /** Get the point in this box closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      return p.cwiseMin(hi).cwiseMax(lo);
    }

    /** Get the closest distance to a point. Returns zero if the box is null. */
    T distance(VectorT const & point) const { return std::sqrt(squaredDistance(point)); }

    /** Get the largest distance between points of this box and a given point. Returns zero if the box is null. */
    T maxDistance(VectorT const & point) const { return std::sqrt(squaredMaxDistance(point)); }

    /** Get the square of the largest distance between points of this box and a given point. Returns zero if the box is null. */
    T squaredMaxDistance(VectorT const & point) const
    {
      if (is_null) return 0;

      VectorT vdiff = (hi - point).cwiseMax(point - lo);
      return vdiff.squaredNorm();
    }

    /** Get the distance between this box and an infinite line. */
    T distance(LineN<N, T> const & line) const { return std::sqrt(squaredDistance(line)); }

    /** Get the squared distance between this box and an infinite line, and optionally return the closest pair of points. */
    T squaredDistance(LineN<N, T> const & line, VectorT * this_pt = nullptr, VectorT * line_pt = nullptr) const
    {
      throw FatalError(format("AxisAlignedBoxN: Distance between line and box not implemented in %ld dimension(s)", N));
    }

    /** Get the distance between this box and a line segment. */
    T distance(LineSegmentN<N, T> const & seg) const { return std::sqrt(squaredDistance(seg)); }

    /** Get the squared distance between this box and a line segment, and optionally return the closest pair of points. */
    T squaredDistance(LineSegmentN<N, T> const & seg, VectorT * this_pt = nullptr, VectorT * seg_pt = nullptr) const
    {
      throw FatalError(format("AxisAlignedBoxN: Distance between line segment and box not implemented in %ld dimension(s)", N));
    }

    /** Get the distance between this box and a ray. */
    T distance(RayN<N, T> const & ray) const { return std::sqrt(squaredDistance(ray)); }

    /** Get the squared distance between this box and a ray, and optionally return the closest pair of points. */
    T squaredDistance(RayN<N, T> const & ray, VectorT * this_pt = nullptr, VectorT * ray_pt = nullptr) const
    {
      throw FatalError(format("AxisAlignedBoxN: Distance between ray and box not implemented in %ld dimension(s)", N));
    }

    /** Get a string representing the box as "[low, high]", or "[null]" if the box is null. */
    std::string toString() const
    {
      if (is_null)
        return "[null]";

      std::ostringstream oss;
      oss << "[" << Thea::toString(lo) << ", " << Thea::toString(hi) << "]";
      return oss.str();
    }

    bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    {
      // TODO: Speed this up: see G3D::Intersect::rayAABox
      return (rayIntersectionTime(ray) > -0.5);
    }

    T rayIntersectionTime(RayN<N, T> const & ray, T max_time = -1) const
    {
      // TODO: Speed this up: see G3D::Intersect::rayAABox
      RayIntersectionN<N, T> isec = rayIntersection(ray, max_time);
      return isec.getTime();
    }

    RayIntersectionN<N, T> rayIntersection(RayN<N, T> const & ray, T max_time = -1) const
    {
      // Early exit
      if (is_null)
        return RayIntersectionN<N, T>(-1);
      else if (max_time >= 0)
      {
        T ray_sqlen = ray.getDirection().squaredNorm();
        T origin_sqdist = squaredDistance(ray.getOrigin());  // fast operation
        if (origin_sqdist > max_time * max_time * ray_sqlen)
          return RayIntersectionN<N, T>(-1);
      }

      VectorT max_t(-1), location;
      VectorT const & origin = ray.getOrigin();
      VectorT const & dir = ray.getDirection();
      bool inside = true;

      // Find candidate planes.
      for (intx i = 0; i < N; ++i)
      {
        if (origin[i] < lo[i])
        {
          location[i] = lo[i];
          inside = false;

          // Calculate distance to candidate plane along this coordinate
          if (Math::fuzzyNe(dir[i], static_cast<T>(0)))
            max_t[i] = (lo[i] - origin[i]) / dir[i];
        }
        else if (origin[i] > hi[i])
        {
          location[i] = hi[i];
          inside = false;

          // Calculate distance to candidate plane along this coordinate
          if (Math::fuzzyNe(dir[i], static_cast<T>(0)))
            max_t[i] = (hi[i] - origin[i]) / dir[i];
        }
      }

      if (inside)  // ray origin inside bounding box
        return RayIntersectionN<N, T>(0);  // undefined normal inside box

      // Get largest of the max_t's for final choice of intersection
      int which_plane = 0;
      for (intx i = 1; i < N; ++i)
        if (max_t[i] > max_t[which_plane])
          which_plane = i;

      // Check final candidate actually inside box
      if (max_t[which_plane] < static_cast<T>(-0.5))  // miss the box
        return RayIntersectionN<N, T>(-1);

      for (intx i = 0; i < N; ++i)
      {
        if (i != which_plane)
        {
          location[i] = origin[i] + max_t[which_plane] * dir[i];
          if (location[i] < lo[i] || location[i] > hi[i])  // on this plane we're outside the box extents, so we miss the box
            return RayIntersectionN<N, T>(-1);
        }
      }

      // Choose the normal to be the plane normal facing into the ray
      VectorT normal = VectorT::Zero();
      normal[which_plane] = (dir[which_plane] > 0) ? -1 : 1;

      return RayIntersectionN<N, T>(max_t[which_plane], &normal);
    }

  private:
    VectorT lo, hi;
    bool is_null;

}; // class AxisAlignedBoxNBase

} // namespace Internal

/** An axis-aligned box in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <int N, typename T = Real>
class /* THEA_API */ AxisAlignedBoxN : public Internal::AxisAlignedBoxNBase<N, T>
{
  private:
    typedef Internal::AxisAlignedBoxNBase<N, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    /** Default constructor, creates a null box. */
    AxisAlignedBoxN() : BaseT() {}

    /** Constructor. Sets the box to be a single point. */
    AxisAlignedBoxN(VectorT const & v) : BaseT(v) {}

    /** Constructor. Sets the extents of the box. */
    AxisAlignedBoxN(VectorT const & lo_, VectorT const & hi_) : BaseT(lo_, hi_) {}

}; // class AxisAlignedBoxN

} // namespace Thea

#include "AxisAlignedBox3.hpp"

#endif
