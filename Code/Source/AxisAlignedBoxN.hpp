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
//============== License for line-AAB distance calculation ===================
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
//
//============================================================================

#ifndef __Thea_AxisAlignedBoxN_hpp__
#define __Thea_AxisAlignedBoxN_hpp__

#include "Common.hpp"
#include "LineN.hpp"
#include "LineSegmentN.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"
#include "RayN.hpp"
#include <bitset>

namespace Thea {

// Forward declarations
template <int N, typename T> class AxisAlignedBoxN;

namespace Internal {

template <typename T>
struct LineAABDistanceResult
{
  T sqrDistance;
  T lineParameter;
};

template <int N, typename T> class LineAABDistance
{
  public:
    typedef Vector<N, T> VectorT;
    typedef LineAABDistanceResult<T> Result;

    // Compute the distance and closest point between a line and an axis-aligned box whose center is the origin.  On input,
    // 'point' is the line origin and 'direction' is the line direction.  On output, 'point' is the point on the box closest to
    // the line.  The 'direction' is non-const to allow transforming the problem into the first octant.
    static void DoQuery(VectorT & point, VectorT & direction, VectorT const & extent, Result & result)
    {
      throw FatalError(format("AxisAlignedBoxN: Distance between line and box not implemented in %ld dimension(s)", N));
    }

    static bool Intersect(VectorT const & point, VectorT const & direction, VectorT const & extent)
    {
      throw FatalError(format("AxisAlignedBoxN: Line-box intersection not implemented in %ld dimension(s)", N));
    }
};

// From https://github.com/davideberly/GeometricTools/blob/master/GTE/Mathematics/DistLine2AlignedBox2.h
//  and https://github.com/davideberly/GeometricTools/blob/master/GTE/Mathematics/IntrLine2AlignedBox2.h
template <typename T>
class LineAABDistance<2, T>
{
  public:
    typedef Vector<2, T> VectorT;
    typedef LineAABDistanceResult<T> Result;

    static void DoQuery(VectorT & point, VectorT & direction, VectorT const & extent, Result & result);
    static bool Intersect(VectorT const & point, VectorT const & direction, VectorT const & extent);

private:
    static void DoQuery2D(VectorT & origin, VectorT const & direction, VectorT const & extent, Result & result);
    static void DoQuery1D(int i0, int i1, VectorT & origin, VectorT const & direction, VectorT const & extent,
                          Result & result);
    static void DoQuery0D(VectorT & origin, VectorT const & extent, Result& result);

    // Compute Dot((x0,x1),Perp(y0,y1)) = x0*y1 - x1*y0, where v0 = (x0,x1) and v1 = (y0,y1).
    static T DotPerp(VectorT const & v0, VectorT const & v1) { return v0[0] * v1[1] - v0[1] * v1[0]; }
};

// From https://github.com/davideberly/GeometricTools/blob/master/GTE/Mathematics/DistLine3CanonicalBox3.h
//  and https://github.com/davideberly/GeometricTools/blob/master/GTE/Mathematics/IntrLine3AlignedBox3.h
template <typename T>
class LineAABDistance<3, T>
{
  public:
    typedef Vector<3, T> VectorT;
    typedef LineAABDistanceResult<T> Result;

    static void DoQuery(VectorT & point, VectorT & direction, VectorT const & extent, Result & result);
    static bool Intersect(VectorT const & point, VectorT const & direction, VectorT const & extent);

  private:
    static void Face(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & PmE,
                     VectorT const & extent, Result & result);
    static void DoQuery3D(VectorT & origin, VectorT const & direction, VectorT const & extent, Result & result);
    static void DoQuery2D(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & extent,
                          Result & result);
    static void DoQuery1D(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & extent,
                          Result & result);
    static void DoQuery0D(VectorT & origin, VectorT const & extent, Result & result);
};

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
     * Get the i'th corner vertex of the box, for \a i in the range [0, 2^N - 1]. This function works as expected only if
     * N <= sizeof(uintx). The d'th coordinate of the returned point is assigned the d'th coordinate of the maximum corner of
     * the box if the d'th bit of i (where the 0'th bit is the least significant) is 1, and of the minimum corner if it is 0.
     */
    VectorT getVertex(uintx i) const
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

    /** Check if the box intersects a line. */
    bool intersects(LineN<N, T> const & line) const
    {
      // Translate the line and box so that the box has center at the origin
      VectorT box_center = this->getCenter(), box_extent = static_cast<T>(0.5) * this->getExtent();
      VectorT point = line.getPoint() - box_center;
      VectorT direction = line.getDirection();

      return Internal::LineAABDistance<N, T>::Intersect(point, direction, box_extent);
    }

    /** Check if the box intersects a line segment. */
    bool intersects(LineSegmentN<N, T> const & seg) const
    {
      // Translate the segment and box so that the box has center at the origin
      VectorT box_center = this->getCenter(), box_extent = static_cast<T>(0.5) * this->getExtent();
      VectorT seg_center = seg.getPoint(static_cast<T>(0.5)) - box_center;
      T seg_extent = seg.length();
      VectorT seg_dir = (seg_extent > std::numeric_limits<T>::min() ? VectorT(seg.getDirection() / seg_extent)
                                                                    : VectorT(VectorT::UnitX()));
      seg_extent /= 2;

      for (int i = 0; i < N; ++i)
        if (std::fabs(seg_center[i]) > box_extent[i] + seg_extent * std::fabs(seg_dir[i]))
          return false;

      return Internal::LineAABDistance<N, T>::Intersect(seg_center, seg_dir, box_extent);
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

      VectorT vdiff = vmin - vmax;              // any axes with overlap have negative values here
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
      // Translate the line and box so that the box has center at the origin
      VectorT box_center = this->getCenter(), box_extent = static_cast<T>(0.5) * this->getExtent();
      VectorT point = line.getPoint() - box_center;
      VectorT direction = line.getDirection();

      typename Internal::LineAABDistance<N, T>::Result result;
      Internal::LineAABDistance<N, T>::DoQuery(point, direction, box_extent, result);

      if (this_pt) *this_pt = box_center + point;
      if (line_pt) *line_pt = line.getPoint() + result.lineParameter * line.getDirection();

      return result.sqrDistance;
    }

    /** Get the distance between this box and a line segment. */
    T distance(LineSegmentN<N, T> const & seg) const { return std::sqrt(squaredDistance(seg)); }

    /** Get the squared distance between this box and a line segment, and optionally return the closest pair of points. */
    T squaredDistance(LineSegmentN<N, T> const & seg, VectorT * this_pt = nullptr, VectorT * seg_pt = nullptr) const
    {
      // Translate the segment and box so that the box has center at the origin
      VectorT box_center = this->getCenter(), box_extent = static_cast<T>(0.5) * this->getExtent();
      VectorT e0 = seg.getEndpoint(0);
      VectorT point = e0 - box_center;
      VectorT seg_dir = seg.getDirection();  // no need for normalization
      VectorT direction = seg_dir;  // can be overwritten by DoQuery

      typename Internal::LineAABDistance<N, T>::Result result;
      Internal::LineAABDistance<N, T>::DoQuery(point, direction, box_extent, result);

      if (result.lineParameter >= 0 && result.lineParameter <= 1)
      {
        if (this_pt) *this_pt = box_center + point;
        if (seg_pt)  *seg_pt  = e0 + result.lineParameter * seg_dir;
      }
      else
      {
        VectorT e1 = seg.getEndpoint(1);
        VectorT c0 = this->closestPoint(e0);
        VectorT c1 = this->closestPoint(e1);
        Real sqdist0 = (c0 - e0).squaredNorm();
        Real sqdist1 = (c1 - e1).squaredNorm();
        if (sqdist0 < sqdist1)
        {
          if (this_pt) *this_pt = c0;
          if (seg_pt)  *seg_pt  = e0;
          result.sqrDistance = sqdist0;
        }
        else
        {
          if (this_pt) *this_pt = c1;
          if (seg_pt)  *seg_pt  = e1;
          result.sqrDistance = sqdist1;
        }
      }

      return result.sqrDistance;
    }

    /** Get the distance between this box and a ray. */
    T distance(RayN<N, T> const & ray) const { return std::sqrt(squaredDistance(ray)); }

    /** Get the squared distance between this box and a ray, and optionally return the closest pair of points. */
    T squaredDistance(RayN<N, T> const & ray, VectorT * this_pt = nullptr, VectorT * ray_pt = nullptr) const
    {
      // Translate the ray and box so that the box has center at the origin
      VectorT box_center = this->getCenter(), box_extent = static_cast<T>(0.5) * this->getExtent();
      VectorT point = ray.getOrigin() - box_center;
      VectorT ray_dir = ray.getDirection();  // no need for normalization
      VectorT direction = ray_dir;  // can be overwritten by DoQuery

      typename Internal::LineAABDistance<N, T>::Result result;
      Internal::LineAABDistance<N, T>::DoQuery(point, direction, box_extent, result);

      if (result.lineParameter >= 0)
      {
        if (this_pt) *this_pt = box_center + point;
        if (ray_pt)  *ray_pt  = ray.getOrigin() + result.lineParameter * ray_dir;
      }
      else
      {
        VectorT c = this->closestPoint(ray.getOrigin());
        result.sqrDistance = (c - ray.getOrigin()).squaredNorm();

        if (this_pt) *this_pt = c;
        if (ray_pt)  *ray_pt  = ray.getOrigin();
      }

      return result.sqrDistance;
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

      VectorT max_t, location; max_t.fill(-1);
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

    /** Get a string representing the box as "[low, high]", or "[null]" if the box is null. */
    std::string toString() const
    {
      if (is_null)
        return "[null]";

      std::ostringstream oss;
      oss << "[" << Thea::toString(lo) << ", " << Thea::toString(hi) << "]";
      return oss.str();
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

namespace Internal {

//=============================================================================================================================
//
// Implementation of line-AAB distance in 2D
//
//=============================================================================================================================

template <typename T>
void
LineAABDistance<2, T>::DoQuery(VectorT & point, VectorT & direction, VectorT const & extent, Result & result)
{
  // Make copies of the origin line origin and direction because we'll need them to compute the distance
  VectorT orig_point = point, orig_dir = direction;

  // Apply reflections so that the direction has nonnegative
  // components.
  T const zero = static_cast<T>(0);
  std::array<bool, 2> reflect{ false, false };
  for (int i = 0; i < 2; ++i)
  {
    if (direction[i] < zero)
    {
      point[i] = -point[i];
      direction[i] = -direction[i];
      reflect[i] = true;
    }
  }

  // Compute the line-box distance and closest points. The DoQueryND
  // calls compute result.lineParameter and put the box's closest
  // point in 'point'.
  if (direction[0] > zero)
  {
    if (direction[1] > zero)
    {
      // The direction signs are (+,+). If the line does not
      // intersect the box, the only possible closest box points
      // are K[0] = (-e0,e1) or K[1] = (e0,-e1). If the line
      // intersects the box, the closest points are the same and
      // chosen to be the intersection with box edge x0 = e0 or
      // x1 = e1. For the remaining discussion, define K[2] =
      // (e0,e1).
      //
      // Test where the candidate corners are relative to the
      // line. If D = (d0,d1), then Perp(D) = (d1,-d0). The
      // corner K[i] = P + t[i] * D + s[i] * Perp(D), where
      // s[i] = Dot(K[i]-P,Perp(D))/|D|^2. K[0] is closest when
      // s[0] >= 0 or K[1] is closest when s[1] <= 0. Otherwise,
      // the line intersects the box. If s[2] >= 0, the common
      // closest point is chosen to be (p0+(e1-p1)*d0/d1,e1). If
      // s[2] < 0, the common closest point is chosen to be
      // (e0,p1+(e0-p0)*d1/d0).
      //
      // It is sufficient to test the signs of Dot(K[i],Perp(D))
      // and defer the division by |D|^2 until needed for
      // computing the closest point.
      DoQuery2D(point, direction, extent, result);
    }
    else
    {
      // The direction signs are (+,0). The parameter is the
      // value of t for which P + t * D = (e0, p1).
      DoQuery1D(0, 1, point, direction, extent, result);
    }
  }
  else
  {
    if (direction[1] > zero)
    {
      // The direction signs are (0,+). The parameter is the
      // value of t for which P + t * D = (p0, e1).
      DoQuery1D(1, 0, point, direction, extent, result);
    }
    else
    {
      // The direction signs are (0,0). The line is degenerate
      // to a point (its point). Clamp the point to the box
      // to obtain the closest point.
      DoQuery0D(point, extent, result);
    }
  }

  // Undo the reflections applied previously.
  for (int i = 0; i < 3; ++i)
  {
    if (reflect[i])
    {
      point[i] = -point[i];
    }
  }

  result.sqrDistance = (orig_point + result.lineParameter * orig_dir - point).squaredNorm();
}

template <typename T>
bool
LineAABDistance<2, T>::Intersect(VectorT const & point, VectorT const & direction, VectorT const & extent)
{
  T LHS = std::fabs(DotPerp(direction, point));
  T RHS = extent[0] * std::fabs(direction[1]) +
          extent[1] * std::fabs(direction[0]);
  return (LHS <= RHS);
}

template <typename T>
void
LineAABDistance<2, T>::DoQuery2D(VectorT & origin, VectorT const & direction, VectorT const & extent, Result & result)
{
  T const zero = static_cast<T>(0);
  VectorT K0{ -extent[0], extent[1] };
  VectorT delta = K0 - origin;
  T K0dotPerpD = DotPerp(delta, direction);
  if (K0dotPerpD >= zero)
  {
    result.lineParameter = delta.dot(direction) / direction.dot(direction);
    origin = K0;
  }
  else
  {
    VectorT K1{ extent[0], -extent[1] };
    delta = K1 - origin;
    T K1dotPerpD = DotPerp(delta, direction);
    if (K1dotPerpD <= zero)
    {
      result.lineParameter = delta.dot(direction) / direction.dot(direction);
      origin = K1;
    }
    else
    {
      VectorT K2{ extent[0], extent[1] };
      delta = K2 - origin;
      T K2dotPerpD = DotPerp(delta, direction);
      if (K2dotPerpD >= zero)
      {
        result.lineParameter = (extent[1] - origin[1]) / direction[1];
        origin[0] += result.lineParameter * direction[0];
        origin[1] =  extent[1];
      }
      else
      {
        result.lineParameter = (extent[0] - origin[0]) / direction[0];
        origin[0] =  extent[0];
        origin[1] += result.lineParameter * direction[1];
      }
    }
  }
}

template <typename T>
void
LineAABDistance<2, T>::DoQuery1D(int i0, int i1, VectorT & origin, VectorT const & direction, VectorT const & extent,
                                 Result & result)
{
  result.lineParameter = (extent[i0] - origin[i0]) / direction[i0];
  origin[i0] = extent[i0];
  origin[i1] = Math::clamp(origin[i1], -extent[i1], extent[i1]);
}

template <typename T>
void
LineAABDistance<2, T>::DoQuery0D(VectorT & origin, VectorT const & extent, Result & result)
{
  result.lineParameter = static_cast<T>(0);
  origin[0] = Math::clamp(origin[0], -extent[0], extent[0]);
  origin[1] = Math::clamp(origin[1], -extent[1], extent[1]);
}

//=============================================================================================================================
//
// Implementation of line-AAB distance in 3D
//
//=============================================================================================================================

template <typename T>
void
LineAABDistance<3, T>::DoQuery(VectorT & point, VectorT & direction, VectorT const & extent, Result & result)
{
  // Copies are made so that we can transform the line direction to
  // the first octant (nonnegative components) using reflections.
  T const zero = static_cast<T>(0);
  std::array<bool, 3> reflect{ false, false, false };
  for (int i = 0; i < 3; ++i)
  {
    if (direction[i] < zero)
    {
      point[i] = -point[i];
      direction[i] = -direction[i];
      reflect[i] = true;
    }
  }

  // Compute the line-box distance and closest points. The DoQueryND
  // calls compute result.lineParameter and update result.sqrDistance.
  // The result.distance is computed after the specialized queries.
  // The closest point on the box is computed afterwards.
  result.sqrDistance = zero;
  result.lineParameter = zero;
  if (direction[0] > zero)
  {
    if (direction[1] > zero)
    {
      if (direction[2] > zero)  // (+,+,+)
      {
        DoQuery3D(point, direction, extent, result);
      }
      else  // (+,+,0)
      {
        DoQuery2D(0, 1, 2, point, direction, extent, result);
      }
    }
    else
    {
      if (direction[2] > zero)  // (+,0,+)
      {
        DoQuery2D(0, 2, 1, point, direction, extent, result);
      }
      else  // (+,0,0)
      {
        DoQuery1D(0, 1, 2, point, direction, extent, result);
      }
    }
  }
  else
  {
    if (direction[1] > zero)
    {
      if (direction[2] > zero)  // (0,+,+)
      {
        DoQuery2D(1, 2, 0, point, direction, extent, result);
      }
      else  // (0,+,0)
      {
        DoQuery1D(1, 0, 2, point, direction, extent, result);
      }
    }
    else
    {
      if (direction[2] > zero)  // (0,0,+)
      {
        DoQuery1D(2, 0, 1, point, direction, extent, result);
      }
      else  // (0,0,0)
      {
        DoQuery0D(point, extent, result);
      }
    }
  }

  // Undo the reflections applied previously.
  for (int i = 0; i < 3; ++i)
  {
    if (reflect[i])
    {
      point[i] = -point[i];
    }
  }
}

template <typename T>
bool
LineAABDistance<3, T>::Intersect(VectorT const & point, VectorT const & direction, VectorT const & extent)
{
  VectorT WxD = direction.cross(point);
  std::array<T, 3> absWdU
  {
    std::fabs(direction[0]),
    std::fabs(direction[1]),
    std::fabs(direction[2])
  };

  return std::fabs(WxD[0]) <= extent[1] * absWdU[2] + extent[2] * absWdU[1]
      && std::fabs(WxD[1]) <= extent[0] * absWdU[2] + extent[2] * absWdU[0]
      && std::fabs(WxD[2]) <= extent[0] * absWdU[1] + extent[1] * absWdU[0];
}

template <typename T>
void
LineAABDistance<3, T>::Face(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & PmE,
                            VectorT const & extent, Result & result)
{
  T const zero = static_cast<T>(0);
  T const two = static_cast<T>(2);

  VectorT PpE = origin + extent;

  if (direction[i0] * PpE[i1] >= direction[i1] * PmE[i0])
  {
    if (direction[i0] * PpE[i2] >= direction[i2] * PmE[i0])
    {
      // v[i1] >= -e[i1], v[i2] >= -e[i2] (distance = 0)
      origin[i0] = extent[i0];
      origin[i1] -= direction[i1] * PmE[i0] / direction[i0];
      origin[i2] -= direction[i2] * PmE[i0] / direction[i0];
      result.lineParameter = -PmE[i0] / direction[i0];
    }
    else
    {
      // v[i1] >= -e[i1], v[i2] < -e[i2]
      T lenSqr = direction[i0] * direction[i0] +
        direction[i2] * direction[i2];
      T tmp = lenSqr * PpE[i1] - direction[i1] * (direction[i0] * PmE[i0] +
        direction[i2] * PpE[i2]);
      if (tmp <= two * lenSqr * extent[i1])
      {
        T t = tmp / lenSqr;
        lenSqr += direction[i1] * direction[i1];
        tmp = PpE[i1] - t;
        T delta = direction[i0] * PmE[i0] + direction[i1] * tmp +
          direction[i2] * PpE[i2];
        result.lineParameter = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + tmp * tmp +
          PpE[i2] * PpE[i2] + delta * result.lineParameter;

        origin[i0] = extent[i0];
        origin[i1] = t - extent[i1];
        origin[i2] = -extent[i2];
      }
      else
      {
        lenSqr += direction[i1] * direction[i1];
        T delta = direction[i0] * PmE[i0] + direction[i1] * PmE[i1] +
          direction[i2] * PpE[i2];
        result.lineParameter = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1]
          + PpE[i2] * PpE[i2] + delta * result.lineParameter;

        origin[i0] = extent[i0];
        origin[i1] = extent[i1];
        origin[i2] = -extent[i2];
      }
    }
  }
  else
  {
    if (direction[i0] * PpE[i2] >= direction[i2] * PmE[i0])
    {
      // v[i1] < -e[i1], v[i2] >= -e[i2]
      T lenSqr = direction[i0] * direction[i0] +
        direction[i1] * direction[i1];
      T tmp = lenSqr * PpE[i2] - direction[i2] * (direction[i0] * PmE[i0] +
        direction[i1] * PpE[i1]);
      if (tmp <= two * lenSqr * extent[i2])
      {
        T t = tmp / lenSqr;
        lenSqr += direction[i2] * direction[i2];
        tmp = PpE[i2] - t;
        T delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] +
          direction[i2] * tmp;
        result.lineParameter = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
          tmp * tmp + delta * result.lineParameter;

        origin[i0] = extent[i0];
        origin[i1] = -extent[i1];
        origin[i2] = t - extent[i2];
      }
      else
      {
        lenSqr += direction[i2] * direction[i2];
        T delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] +
          direction[i2] * PmE[i2];
        result.lineParameter = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
          PmE[i2] * PmE[i2] + delta * result.lineParameter;

        origin[i0] = extent[i0];
        origin[i1] = -extent[i1];
        origin[i2] = extent[i2];
      }
    }
    else
    {
      // v[i1] < -e[i1], v[i2] < -e[i2]
      T lenSqr = direction[i0] * direction[i0] +
        direction[i2] * direction[i2];
      T tmp = lenSqr * PpE[i1] - direction[i1] * (direction[i0] * PmE[i0] +
        direction[i2] * PpE[i2]);
      if (tmp >= zero)
      {
        // v[i1]-edge is closest
        if (tmp <= two * lenSqr * extent[i1])
        {
          T t = tmp / lenSqr;
          lenSqr += direction[i1] * direction[i1];
          tmp = PpE[i1] - t;
          T delta = direction[i0] * PmE[i0] + direction[i1] * tmp +
            direction[i2] * PpE[i2];
          result.lineParameter = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + tmp * tmp +
            PpE[i2] * PpE[i2] + delta * result.lineParameter;

          origin[i0] = extent[i0];
          origin[i1] = t - extent[i1];
          origin[i2] = -extent[i2];
        }
        else
        {
          lenSqr += direction[i1] * direction[i1];
          T delta = direction[i0] * PmE[i0] + direction[i1] * PmE[i1]
            + direction[i2] * PpE[i2];
          result.lineParameter = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1]
            + PpE[i2] * PpE[i2] + delta * result.lineParameter;

          origin[i0] = extent[i0];
          origin[i1] = extent[i1];
          origin[i2] = -extent[i2];
        }
        return;
      }

      lenSqr = direction[i0] * direction[i0] +
        direction[i1] * direction[i1];
      tmp = lenSqr * PpE[i2] - direction[i2] * (direction[i0] * PmE[i0] +
        direction[i1] * PpE[i1]);
      if (tmp >= zero)
      {
        // v[i2]-edge is closest
        if (tmp <= two * lenSqr * extent[i2])
        {
          T t = tmp / lenSqr;
          lenSqr += direction[i2] * direction[i2];
          tmp = PpE[i2] - t;
          T delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] +
            direction[i2] * tmp;
          result.lineParameter = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
            tmp * tmp + delta * result.lineParameter;

          origin[i0] = extent[i0];
          origin[i1] = -extent[i1];
          origin[i2] = t - extent[i2];
        }
        else
        {
          lenSqr += direction[i2] * direction[i2];
          T delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1]
            + direction[i2] * PmE[i2];
          result.lineParameter = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1]
            + PmE[i2] * PmE[i2] + delta * result.lineParameter;

          origin[i0] = extent[i0];
          origin[i1] = -extent[i1];
          origin[i2] = extent[i2];
        }
        return;
      }

      // (v[i1],v[i2])-corner is closest
      lenSqr += direction[i2] * direction[i2];
      T delta = direction[i0] * PmE[i0] + direction[i1] * PpE[i1] +
        direction[i2] * PpE[i2];
      result.lineParameter = -delta / lenSqr;
      result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1]
        + PpE[i2] * PpE[i2] + delta * result.lineParameter;

      origin[i0] = extent[i0];
      origin[i1] = -extent[i1];
      origin[i2] = -extent[i2];
    }
  }
}

template <typename T>
void
LineAABDistance<3, T>::DoQuery3D(VectorT & origin, VectorT const & direction, VectorT const & extent, Result & result)
{
  VectorT PmE = origin - extent;
  T prodDxPy = direction[0] * PmE[1];
  T prodDyPx = direction[1] * PmE[0];

  if (prodDyPx >= prodDxPy)
  {
    T prodDzPx = direction[2] * PmE[0];
    T prodDxPz = direction[0] * PmE[2];
    if (prodDzPx >= prodDxPz)
    {
      // line intersects x = e0
      Face(0, 1, 2, origin, direction, PmE, extent, result);
    }
    else
    {
      // line intersects z = e2
      Face(2, 0, 1, origin, direction, PmE, extent, result);
    }
  }
  else
  {
    T prodDzPy = direction[2] * PmE[1];
    T prodDyPz = direction[1] * PmE[2];
    if (prodDzPy >= prodDyPz)
    {
      // line intersects y = e1
      Face(1, 2, 0, origin, direction, PmE, extent, result);
    }
    else
    {
      // line intersects z = e2
      Face(2, 0, 1, origin, direction, PmE, extent, result);
    }
  }
}

template <typename T>
void
LineAABDistance<3, T>::DoQuery2D(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & extent,
                                 Result & result)
{
  T const zero = static_cast<T>(0);

  T PmE0 = origin[i0] - extent[i0];
  T PmE1 = origin[i1] - extent[i1];
  T prod0 = direction[i1] * PmE0;
  T prod1 = direction[i0] * PmE1;

  if (prod0 >= prod1)
  {
    // line intersects P[i0] = e[i0]
    origin[i0] = extent[i0];

    T PpE1 = origin[i1] + extent[i1];
    T delta = prod0 - direction[i0] * PpE1;
    if (delta >= zero)
    {
      T lenSqr = direction[i0] * direction[i0] +
        direction[i1] * direction[i1];
      result.sqrDistance += delta * delta / lenSqr;
      origin[i1] = -extent[i1];
      result.lineParameter = -(direction[i0] * PmE0 +
        direction[i1] * PpE1) / lenSqr;
    }
    else
    {
      origin[i1] -= prod0 / direction[i0];
      result.lineParameter = -PmE0 / direction[i0];
    }
  }
  else
  {
    // line intersects P[i1] = e[i1]
    origin[i1] = extent[i1];

    T PpE0 = origin[i0] + extent[i0];
    T delta = prod1 - direction[i1] * PpE0;
    if (delta >= zero)
    {
      T lenSqr = direction[i0] * direction[i0] +
        direction[i1] * direction[i1];
      result.sqrDistance += delta * delta / lenSqr;
      origin[i0] = -extent[i0];
      result.lineParameter = -(direction[i0] * PpE0 +
        direction[i1] * PmE1) / lenSqr;
    }
    else
    {
      origin[i0] -= prod1 / direction[i1];
      result.lineParameter = -PmE1 / direction[i1];
    }
  }

  if (origin[i2] < -extent[i2])
  {
    T delta = origin[i2] + extent[i2];
    result.sqrDistance += delta * delta;
    origin[i2] = -extent[i2];
  }
  else if (origin[i2] > extent[i2])
  {
    T delta = origin[i2] - extent[i2];
    result.sqrDistance += delta * delta;
    origin[i2] = extent[i2];
  }
}

template <typename T>
void
LineAABDistance<3, T>::DoQuery1D(int i0, int i1, int i2, VectorT & origin, VectorT const & direction, VectorT const & extent,
                                 Result & result)
{
  result.lineParameter = (extent[i0] - origin[i0]) / direction[i0];

  origin[i0] = extent[i0];
  for (auto i : { i1, i2 })
  {
    if (origin[i] < -extent[i])
    {
      T delta = origin[i] + extent[i];
      result.sqrDistance += delta * delta;
      origin[i] = -extent[i];
    }
    else if (origin[i] > extent[i])
    {
      T delta = origin[i] - extent[i];
      result.sqrDistance += delta * delta;
      origin[i] = extent[i];
    }
  }
}

template <typename T>
void
LineAABDistance<3, T>::DoQuery0D(VectorT & origin, VectorT const & extent, Result & result)
{
  for (int i = 0; i < 3; ++i)
  {
    if (origin[i] < -extent[i])
    {
      T delta = origin[i] + extent[i];
      result.sqrDistance += delta * delta;
      origin[i] = -extent[i];
    }
    else if (origin[i] > extent[i])
    {
      T delta = origin[i] - extent[i];
      result.sqrDistance += delta * delta;
      origin[i] = extent[i];
    }
  }
}

} // namespace Internal

} // namespace Thea

#include "AxisAlignedBox3.hpp"

#endif
