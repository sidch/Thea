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

#ifndef __Thea_Cone3_hpp__
#define __Thea_Cone3_hpp__

#include "Common.hpp"
#include "AxisAlignedBox3.hpp"
#include "LineSegment3.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "Plane3.hpp"
#include "RayIntersectable3.hpp"
#include <sstream>

namespace Thea {

/** A solid 3D cone, whose base is a circular disk and apex is a point perpendicularly above the disk centre. */
class /* THEA_API */ Cone3 : public RayIntersectable3
{
  public:
    THEA_DECL_SMART_POINTERS(Cone3)

    /** Default constructor. Does not initialize anything. */
    Cone3() {}

    /** Initialize with an axis (line segment from the center of the base to the apex) and a base radius. */
    Cone3(Vector3 const & base_center_, Vector3 const & apex_, Real base_radius_)
    : segment(base_center_, apex_), base_radius(base_radius_)
    {
      updateTransform();
    }

    /** Get the center of the base disk. */
    Vector3 getBaseCenter() const { return segment.getEndpoint(0); }

    /** Get the apex point. */
    Vector3 getApex() const { return segment.getEndpoint(1); }

    /** Get the unnormalized direction vector of the cone from the center of the base to the apex. */
    Vector3 const & getDirection() const { return segment.getDirection(); }

    /** Get the length of the cone. */
    Real length() const { return len; }

    /** Get the square of the length of the cone. */
    Real squaredLength() const { return len * len; }

    /** Set the axis of the cone (line from the center of the base to the apex). */
    void setAxis(Vector3 const & base_center_, Vector3 const & apex_)
    {
      segment = LineSegment3(base_center_, apex_);
      updateTransform();
    }

    /** Get the base radius of the cone. */
    Real getBaseRadius() const { return base_radius; }

    /** Set the base radius of the cone. */
    void setBaseRadius(Real base_radius_) { base_radius = base_radius_; updateTransform(); }

    /** Test if this cone intersects (contains) a point. */
    bool intersects(Vector3 const & p) const { return contains(p); }

    /** Check if the cone contains a point. */
    bool contains(Vector3 const & p) const
    {
      Vector3 local_p = world2local * p;
      return local_p.z() >= 0 && local_p.z() <= 1 && local_p.head<2>().squaredNorm() <= Math::square(1 - local_p.z());
    }

    /** Get a bounding box for the cone. */
    AxisAlignedBox3 getBounds() const
    {
      return AxisAlignedBox3(Vector3(-base_radius, -base_radius, 0),
                             Vector3( base_radius,  base_radius, len)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the cone. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[axis: " << segment.toString() << ", base_radius: " << base_radius << ']';
      return oss.str();
    }

    bool rayIntersects(Ray3 const & ray, Real max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time) >= 0;
    }

    Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const
    {
      int surface;
      return rayIntersectionTime(ray, max_time, surface);
    }

    RayIntersection3 rayIntersection(Ray3 const & ray, Real max_time = -1) const
    {
      int surface;
      Real t = rayIntersectionTime(ray, max_time, surface);
      if (surface >= 0)
      {
        Vector3 n;
        switch (surface)
        {
          case 0: n = -unit_dir; break;
          case 1:
          {
            Vector3 p = ray.getPoint(t);
            Vector3 cp = segment.closestPoint(p);
            Vector3 h = getApex() - cp;
            Vector3 d = p - cp;
            Real hlen2 = h.squaredNorm();
            Real dlen2 = d.squaredNorm();
            if (dlen2 > 0)
              n = (h + (hlen2 / dlen2) * d).stableNormalized();  // Exchange lengths of h and d. FIXME: Avoid the final sqrt?
            else
              n = Vector3::Zero();
            break;
          }
          default: n = Vector3::Zero();
        }

        return RayIntersection3(t, &n);
      }
      else
        return RayIntersection3(-1);
    }

  private:
    /** Compute the transform that maps world coordinates to normalized cone coordinates. */
    void updateTransform()
    {
      len = segment.length();
      if (len > 0 && base_radius > 0)
      {
        unit_dir = getDirection().stableNormalized();
        Matrix3 lin = Math::scaling(Vector3(1 / base_radius, 1 / base_radius, 1 / len))
                    * Math::rotationArc(unit_dir, Vector3(0, 0, 1), /* normalize_dirs = */ false);
        world2local = AffineTransform3(lin, -(lin * getBaseCenter()));
      }
      else
      {
        unit_dir = Vector3::Zero();
        world2local.setIdentity();
      }
    }

    /**
     * Helper function for rayIntersectionTime(Ray3 const &, Real) that also returns the hit surface: 0 for bottom, 1 for sides,
     * 2 for interior, negative for no hit.
     */
    Real rayIntersectionTime(Ray3 const & ray, Real max_time, int & surface) const
    {
      if (contains(ray.getOrigin()))
      {
        surface = 2;
        return 0;
      }

      surface = -1;
      Ray3 local_ray = ray.transform(world2local);
      double min_root = -1;

      // Base
      Plane3 plane = Plane3::fromPointAndNormal(Vector3::Zero(), Vector3::UnitZ());
      auto plane_time = plane.rayIntersectionTime(local_ray, max_time);
      if (plane_time >= 0 && (min_root < 0 || plane_time < min_root))
      {
        Vector3 p = local_ray.getPoint(plane_time);
        if (p.head(2).squaredNorm() <= 1)
        {
          min_root = plane_time;
          surface = 0;
        }
      }

      // Sides
      Vector3 p2  = local_ray.getOrigin().cwiseProduct(local_ray.getOrigin());
      Vector3 pu  = local_ray.getOrigin().cwiseProduct(local_ray.getDirection());
      Vector3 u2  = local_ray.getDirection().cwiseProduct(local_ray.getDirection());
      Real diff_z = local_ray.getOrigin().z() - 1;

      double s[3];
      s[2] = u2[0] + u2[1] - u2[2];
      s[1] = 2 * (pu[0] + pu[1] - diff_z * local_ray.getDirection()[2]);
      s[0] = p2[0] + p2[1] - diff_z * diff_z;

      double roots[2];
      int num_roots = Math::solveQuadratic(s[0], s[1], s[2], roots);

      for (int i = 0; i < num_roots; ++i)
        if (roots[i] >= 0 && (max_time < 0 || roots[i] <= max_time) && (min_root < 0 || roots[i] < min_root))
        {
          Vector3 point = local_ray.getPoint((Real)roots[i]);
          if (point[2] >= 0 && point[2] <= 1)
          {
            min_root = roots[i];
            surface = 1;
          }
        }

      return (Real)min_root;
    }

    LineSegment3      segment;      ///< Line segment defining the axis of the cone.
    Real              base_radius;  ///< Radius of the base of the cone.
    Vector3           unit_dir;     ///< Unit vector along the axis of the cone.
    Real              len;          ///< Length of the cone (center of base to apex).
    AffineTransform3  world2local;  ///< Transform from world coordinates to normalized cone coordinates.

}; // class Cone3

} // namespace Thea

#endif
