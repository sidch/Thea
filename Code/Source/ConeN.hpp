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
// First version: 2021
//
//============================================================================

#ifndef __Thea_ConeN_hpp__
#define __Thea_ConeN_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"
#include "AxisAlignedBoxN.hpp"
#include "LineSegmentN.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "HyperplaneN.hpp"
#include "RayIntersectableN.hpp"
#include <sstream>

namespace Thea {

/**
 * A solid cone in N dimensions, defined as the convex hull of an (N - 1)-dimensional ball and a point. Currently, this is only
 * implemented for N = 3.
 */
template <int N, typename T = Real> class /* THEA_API */ ConeN;

/** A solid 3D cone, whose base is a circular disk and apex is a point perpendicularly above the disk centre. */
template <typename T>
class /* THEA_API */ ConeN<3, T> : public RayIntersectableN<3, T>
{
  private:
    typedef LineSegmentN<3, T>      LineSegmentT;      ///< Compatible line segment type.
    typedef AffineTransformN<3, T>  AffineTransformT;  ///< Compatible affine transform type.

  public:
    THEA_DECL_SMART_POINTERS(ConeN)

    typedef Vector<3, T> VectorT;  ///< Compatible vector type.

    /** Default constructor. Does not initialize anything. */
    ConeN() {}

    /** Initialize with an axis (line segment from the center of the base to the apex) and a base radius. */
    ConeN(VectorT const & base_center_, VectorT const & apex_, T base_radius_)
    : segment(base_center_, apex_), base_radius(base_radius_)
    {
      updateTransform();
    }

    /** Get the center of the base disk. */
    VectorT getBaseCenter() const { return segment.getEndpoint(0); }

    /** Get the apex point. */
    VectorT getApex() const { return segment.getEndpoint(1); }

    /** Get the unnormalized direction vector of the cone from the center of the base to the apex. */
    VectorT const & getDirection() const { return segment.getDirection(); }

    /** Get the length of the cone. */
    T length() const { return len; }

    /** Get the square of the length of the cone. */
    T squaredLength() const { return len * len; }

    /** Set the axis of the cone (line from the center of the base to the apex). */
    void setAxis(VectorT const & base_center_, VectorT const & apex_)
    {
      segment = LineSegmentT(base_center_, apex_);
      updateTransform();
    }

    /** Get the base radius of the cone. */
    T getBaseRadius() const { return base_radius; }

    /** Set the base radius of the cone. */
    void setBaseRadius(T base_radius_) { base_radius = base_radius_; updateTransform(); }

    /** Test if this cone intersects (contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the cone contains a point. */
    bool contains(VectorT const & p) const
    {
      VectorT local_p = world2local * p;
      return local_p.z() >= 0 && local_p.z() <= 1 && local_p.template head<2>().squaredNorm() <= Math::square(1 - local_p.z());
    }

    /** Get a bounding box for the cone. */
    AxisAlignedBoxN<3, T> getBounds() const
    {
      return AxisAlignedBoxN<3, T>(VectorT(-base_radius, -base_radius, 0),
                                   VectorT( base_radius,  base_radius, len)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the cone. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[axis: " << segment.toString() << ", base_radius: " << base_radius << ']';
      return oss.str();
    }

    bool rayIntersects(RayN<3, T> const & ray, T max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time) >= 0;
    }

    T rayIntersectionTime(RayN<3, T> const & ray, T max_time = -1) const
    {
      int surface;
      return rayIntersectionTime(ray, max_time, surface);
    }

    RayIntersectionN<3, T> rayIntersection(RayN<3, T> const & ray, T max_time = -1) const
    {
      int surface;
      T t = rayIntersectionTime(ray, max_time, surface);
      if (surface >= 0)
      {
        VectorT n;
        switch (surface)
        {
          case 0: n = -unit_dir; break;
          case 1:
          {
            VectorT p = ray.getPoint(t);
            VectorT cp = segment.closestPoint(p);
            VectorT h = getApex() - cp;
            VectorT d = p - cp;
            T hlen2 = h.squaredNorm();
            T dlen2 = d.squaredNorm();
            if (dlen2 > 0)
              n = (h + (hlen2 / dlen2) * d).stableNormalized();  // Exchange lengths of h and d. FIXME: Avoid the final sqrt?
            else
              n = VectorT::Zero();
            break;
          }
          default: n = VectorT::Zero();
        }

        return RayIntersectionN<3, T>(t, &n);
      }
      else
        return RayIntersectionN<3, T>(-1);
    }

  private:
    /** Compute the transform that maps world coordinates to normalized cone coordinates. */
    void updateTransform()
    {
      len = segment.length();
      if (len > 0 && base_radius > 0)
      {
        unit_dir = getDirection().stableNormalized();
        Matrix<3, 3, T> lin = Math::scaling(VectorT(1 / base_radius, 1 / base_radius, 1 / len))
                            * Math::rotationArc(unit_dir, VectorT(0, 0, 1), /* normalize_dirs = */ false);
        world2local = AffineTransformT(lin, -(lin * getBaseCenter()));
      }
      else
      {
        unit_dir = VectorT::Zero();
        world2local.setIdentity();
      }
    }

    /**
     * Helper function for rayIntersectionTime(RayN<3, T> const &, T) that also returns the hit surface: 0 for bottom, 1 for
     * sides, 2 for interior, negative for no hit.
     */
    T rayIntersectionTime(RayN<3, T> const & ray, T max_time, int & surface) const
    {
      if (contains(ray.getOrigin()))
      {
        surface = 2;
        return 0;
      }

      surface = -1;
      RayN<3, T> local_ray = ray.transform(world2local);
      double min_root = -1;

      // Base
      auto plane = HyperplaneN<3, T>::fromPointAndNormal(VectorT::Zero(), VectorT::UnitZ());
      auto plane_time = plane.rayIntersectionTime(local_ray, max_time);
      if (plane_time >= 0 && (min_root < 0 || plane_time < min_root))
      {
        VectorT p = local_ray.getPoint(plane_time);
        if (p.head(2).squaredNorm() <= 1)
        {
          min_root = plane_time;
          surface = 0;
        }
      }

      // Sides
      auto lrd = local_ray.template cast<double>();
      Vector<3, double> p2 = lrd.getOrigin().cwiseProduct(lrd.getOrigin());
      Vector<3, double> pu = lrd.getOrigin().cwiseProduct(lrd.getDirection());
      Vector<3, double> u2 = lrd.getDirection().cwiseProduct(lrd.getDirection());
      double diff_z = lrd.getOrigin().z() - 1;

      double s[3];
      s[2] = u2[0] + u2[1] - u2[2];
      s[1] = 2 * (pu[0] + pu[1] - diff_z * local_ray.getDirection()[2]);
      s[0] = p2[0] + p2[1] - diff_z * diff_z;

      double roots[2];
      int num_roots = Math::solveQuadratic(s[0], s[1], s[2], roots);

      for (int i = 0; i < num_roots; ++i)
        if (roots[i] >= 0 && (max_time < 0 || roots[i] <= max_time) && (min_root < 0 || roots[i] < min_root))
        {
          VectorT point = local_ray.getPoint((T)roots[i]);
          if (point[2] >= 0 && point[2] <= 1)
          {
            min_root = roots[i];
            surface = 1;
          }
        }

      return (T)min_root;
    }

    LineSegmentT      segment;      ///< Line segment defining the axis of the cone.
    T                 base_radius;  ///< Radius of the base of the cone.
    VectorT           unit_dir;     ///< Unit vector along the axis of the cone.
    T                 len;          ///< Length of the cone (center of base to apex).
    AffineTransformT  world2local;  ///< Transform from world coordinates to normalized cone coordinates.

}; // class ConeN<3, T>

} // namespace Thea

#include "Cone3.hpp"

#endif
