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

#ifndef __Thea_Torus3_hpp__
#define __Thea_Torus3_hpp__

#include "Common.hpp"
#include "AxisAlignedBox3.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayIntersectable3.hpp"
#include <sstream>

namespace Thea {

/**
 * A solid 3D torus, traced by a circular disk whose radius is the inner radius of the torus, and whose center moves along a
 * larger circle whose radius is the outer radius of the torus.
 */
class /* THEA_API */ Torus3 : public RayIntersectable3
{
  public:
    THEA_DECL_SMART_POINTERS(Torus3)

    /** Default constructor. Does not initialize anything. */
    Torus3() {}

    /**
     * Initialize with the center point, the axis (unnormalized direction vector orthogonal to the plane of the torus), and
     * outer and inner radii.
     */
    Torus3(Vector3 const & center_, Vector3 const & axis_, Real outer_radius_, Real inner_radius_)
    : center(center_), axis(axis_.stableNormalized()), outer_radius(outer_radius_), inner_radius(inner_radius_)
    {
      updateTransform();
    }

    /** Get the center of the torus. */
    Vector3 const & getCenter() const { return center; }

    /** Set the center of the torus. */
    void setCenter(Vector3 const & center_) { center = center_; updateTransform(); }

    /** Get the axis of the torus, normalized to unit length. */
    Vector3 const & getAxis() const { return axis; }

    /** Set the axis of the torus from an unnormalized direction vector. */
    void setAxis(Vector3 const & axis_) { axis = axis_.stableNormalized(); updateTransform(); }

    /** Get the outer radius of the torus. */
    Real getOuterRadius() const { return outer_radius; }

    /** Set the outer radius of the torus. */
    void setOuterRadius(Real outer_radius_) { outer_radius = outer_radius_; }

    /** Get the inner radius of the torus. */
    Real getInnerRadius() const { return inner_radius; }

    /** Set the inner radius of the torus. */
    void setInnerRadius(Real inner_radius_) { inner_radius = inner_radius_; }

    /** Test if this torus intersects (contains) a point. */
    bool intersects(Vector3 const & p) const { return contains(p); }

    /** Check if the torus contains a point. */
    bool contains(Vector3 const & p) const
    {
      Vector3 local_p = world2local * p;
      Real a2 = Math::square(outer_radius);
      Real b2 = Math::square(inner_radius);
      return Math::square(local_p.squaredNorm() + a2 - b2) < 4 * a2 * (local_p.head<2>().squaredNorm());
    }

    /** Get a bounding box for the torus. */
    AxisAlignedBox3 getBounds() const
    {
      Real s = outer_radius + inner_radius;
      return AxisAlignedBox3(Vector3(-s, -s, -inner_radius),
                             Vector3( s,  s,  inner_radius)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the torus. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[center: " << Thea::toString(center) << ", axis: " << Thea::toString(axis)
          << ", outer_radius: " << outer_radius << ", inner_radius: " << inner_radius << ']';
      return oss.str();
    }

    bool rayIntersects(Ray3 const & ray, Real max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time) >= 0;
    }

    Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time, nullptr);
    }

    RayIntersection3 rayIntersection(Ray3 const & ray, Real max_time = -1) const
    {
      Vector3 n;
      Real t = rayIntersectionTime(ray, max_time, &n);
      if (t >= 0)
        return RayIntersection3(t, &n);
      else
        return RayIntersection3(-1);
    }

  private:
    /** Compute the transform that maps world coordinates to local torus coordinates. */
    void updateTransform()
    {
      // Transform is orthonormal -- we exploit this in other places
      Matrix3 lin = Math::rotationArc(axis, Vector3(0, 0, 1), /* normalize_dirs = */ false);
      world2local = AffineTransform3(lin, -(lin * center));
    }

    /**
     * Helper function for rayIntersectionTime(Ray3 const &, Real) that also optionally returns the surface normal at the hit
     * point.
     */
    Real rayIntersectionTime(Ray3 const & ray, Real max_time, Vector3 * normal) const
    {
      if (contains(ray.getOrigin()))
      {
        if (normal) *normal = Vector3::Zero();
        return 0;
      }

      Ray3 local_ray = ray.transform(world2local);
      double min_root = -1;

      double r2pw2 = outer_radius * outer_radius + inner_radius * inner_radius;
      double r2mw2 = r2pw2 - 2 * inner_radius * inner_radius;

      Vector3 p2 = local_ray.getOrigin().cwiseProduct(local_ray.getOrigin());
      Vector3 pu = local_ray.getOrigin().cwiseProduct(local_ray.getDirection());
      Vector3 u2 = local_ray.getDirection().cwiseProduct(local_ray.getDirection());

      double s[5];
      s[4] = u2[0] * (u2[0] + 2 * u2[1]) + u2[1] * (u2[1] + 2 * u2[2]) + u2[2] * (u2[2] + 2 * u2[0]);
      s[3] = 4 * (pu[0] + pu[1] + pu[2]) * (u2[0] + u2[1] + u2[2]);
      s[2] = 2 * (r2mw2 * u2[2] - r2pw2 * (u2[0] + u2[1]))
           + 8 * (pu[0] * pu[1] + pu[1] * pu[2] + pu[2] * pu[0])
           + 6 * (pu[0] * pu[0] + pu[1] * pu[1] + pu[2] * pu[2])
           + 2 * (p2[0] * (u2[1] + u2[2]) + p2[1] * (u2[2] + u2[0]) + p2[2] * (u2[0] + u2[1]));
      s[1] = 4 * (r2mw2 * pu[2] - r2pw2 * (pu[0] + pu[1]) + (p2[0] + p2[1] + p2[2]) * (pu[0] + pu[1] + pu[2]));
      s[0] = 2 * (r2mw2 * p2[2] - r2pw2 * (p2[0] + p2[1]) + p2[0] * p2[1] + p2[1] * p2[2] + p2[2] * p2[0])
           + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2] + r2mw2 * r2mw2;

      double roots[4];
      int num_roots = Math::solveQuartic(s[0], s[1], s[2], s[3], s[4], roots);

      for (int i = 0; i < num_roots; ++i)
        if (roots[i] >= 0 && (max_time < 0 || roots[i] <= max_time) && (min_root < 0 || roots[i] < min_root))
          min_root = roots[i];

      if (normal && min_root >= 0)
      {
        // Exploit orthonormality of world2local
        Vector3 local_p = local_ray.getPoint((Real)min_root);
        Vector2 ring_point = outer_radius * local_p.head<2>().stableNormalized();
        *normal = (world2local.getLinear() * (local_p - Vector3(ring_point.x(), ring_point.y(), 0))).stableNormalized();
      }

      return (Real)min_root;
    }

    Vector3           center;         ///< The center of the torus.
    Vector3           axis;           ///< Unit vector along the axis of the torus (perpendicular to main disk).
    Real              outer_radius;   ///< Outer radius of the torus.
    Real              inner_radius;   ///< Inner radius of the torus.
    AffineTransform3  world2local;    ///< Transform from world coordinates to local torus coordinates (axis along +Z).

}; // class Torus3

} // namespace Thea

#endif
