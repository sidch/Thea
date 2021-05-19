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

#ifndef __Thea_Cylinder3_hpp__
#define __Thea_Cylinder3_hpp__

#include "Common.hpp"
#include "AxisAlignedBox3.hpp"
#include "LineSegment3.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "Plane3.hpp"
#include "RayIntersectable3.hpp"
#include <sstream>

namespace Thea {

/** A solid 3D cylinder, traced by a circular disk whose center moves along a line segment. */
class /* THEA_API */ Cylinder3 : public RayIntersectable3
{
  public:
    THEA_DECL_SMART_POINTERS(Cylinder3)

    /** Default constructor. Does not initialize anything. */
    Cylinder3() {}

    /** Initialize with an axis (line segment from the center of one end to the center of the other end) and a radius. */
    Cylinder3(Vector3 const & begin_, Vector3 const & end_, Real radius_) : segment(begin_, end_), radius(radius_)
    {
      updateTransform();
    }

    /** Get the center of one capping disk of the cylinder: 0 returns the first endpoint and 1 returns the second. */
    Vector3 getEndpoint(int i) const { return segment.getEndpoint(i); }

    /** Get the unnormalized direction vector of the cylinder from the center of one end to the center of the other. */
    Vector3 const & getDirection() const { return segment.getDirection(); }

    /** Get the length of the cylinder. */
    Real length() const { return len; }

    /** Get the square of the length of the cylinder. */
    Real squaredLength() const { return len * len; }

    /** Set the axis of the cylinder. */
    void setAxis(Vector3 const & begin_, Vector3 const & end_)
    {
      segment = LineSegment3(begin_, end_);
      updateTransform();
    }

    /** Get the radius of the cylinder. */
    Real getRadius() const { return radius; }

    /** Set the radius of the cylinder. */
    void setRadius(Real radius_) { radius = radius_; updateTransform(); }

    /** Test if this cylinder intersects (contains) a point. */
    bool intersects(Vector3 const & p) const { return contains(p); }

    /** Check if the cylinder contains a point. */
    bool contains(Vector3 const & p) const
    {
      Vector3 local_p = world2local * p;
      return local_p.z() >= 0 && local_p.z() <= 1 && local_p.squaredNorm() <= 1;
    }

    /** Get a bounding box for the cylinder. */
    AxisAlignedBox3 getBounds() const
    {
      return AxisAlignedBox3(Vector3(-radius, -radius, 0),
                             Vector3( radius,  radius, len)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the cylinder. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[axis: " << segment.toString() << ", radius: " << radius << ']';
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
          case 1: n =  unit_dir; break;
          case 2:
          {
            Vector3 p = ray.getPoint(t);
            n = (p - segment.closestPoint(p)) / radius;
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
    /** Compute the transform that maps world coordinates to normalized cylinder coordinates. */
    void updateTransform()
    {
      len = segment.length();
      if (len > 0 && radius > 0)
      {
        unit_dir = getDirection().stableNormalized();
        Matrix3 lin = Math::scaling(Vector3(1 / radius, 1 / radius, 1 / len))
                    * Math::rotationArc(unit_dir, Vector3(0, 0, 1), /* normalize_dirs = */ false);
        world2local = AffineTransform3(lin, -(lin * getEndpoint(0)));
      }
      else
      {
        unit_dir = Vector3::Zero();
        world2local.setIdentity();
      }
    }

    /**
     * Helper function for rayIntersectionTime(Ray3 const &, Real) that also returns the hit surface: 0 for bottom, 1 for top,
     * 2 for sides, 3 for interior, negative for no hit.
     */
    Real rayIntersectionTime(Ray3 const & ray, Real max_time, int & surface) const
    {
      if (contains(ray.getOrigin()))
      {
        surface = 3;
        return 0;
      }

      surface = -1;
      Ray3 local_ray = ray.transform(world2local);
      double min_root = -1;

      // Bottom and top
      for (int i = 0; i <= 1; ++i)
      {
        Plane3 plane = Plane3::fromPointAndNormal(Vector3(0, 0, i), Vector3::UnitZ());
        auto plane_time = plane.rayIntersectionTime(local_ray, max_time);
        if (plane_time >= 0 && (min_root < 0 || plane_time < min_root))
        {
          Vector3 p = local_ray.getPoint(plane_time);
          if (p.head(2).squaredNorm() <= 1)
          {
            min_root = plane_time;
            surface = i;
          }
        }
      }

      // Sides
      Vector3 p2 = local_ray.getOrigin().cwiseProduct(local_ray.getOrigin());
      Vector3 pu = local_ray.getOrigin().cwiseProduct(local_ray.getDirection());
      Vector3 u2 = local_ray.getDirection().cwiseProduct(local_ray.getDirection());

      double s[3];
      s[2] = u2[0] + u2[1];
      s[1] = 2 * (pu[0] + pu[1]);
      s[0] = p2[0] + p2[1] - 1;

      double roots[2];
      int num_roots = Math::solveQuadratic(s[0], s[1], s[2], roots);

      for (int i = 0; i < num_roots; ++i)
        if (roots[i] >= 0 && (max_time < 0 || roots[i] <= max_time) && (min_root < 0 || roots[i] < min_root))
        {
          Vector3 point = local_ray.getPoint((Real)roots[i]);
          if (point[2] >= 0 && point[2] <= 1)
          {
            min_root = roots[i];
            surface = 2;
          }
        }

      return (Real)min_root;
    }

    LineSegment3      segment;      ///< Line segment defining the axis of the cylinder.
    Real              radius;       ///< Radius of the cylinder.
    Vector3           unit_dir;     ///< Unit vector along the axis of the cylinder.
    Real              len;          ///< Length of the cylinder.
    AffineTransform3  world2local;  ///< Transform from world coordinates to normalized cylinder coordinates.

}; // class Cylinder3

} // namespace Thea

#endif
