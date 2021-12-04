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

#ifndef __Thea_TorusN_hpp__
#define __Thea_TorusN_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"
#include "AxisAlignedBoxN.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"
#include <sstream>

namespace Thea {

/**
 * A solid torus in N dimensions, defined as the volume enclosed by the Cartesian product of N - 1 circles. Currently, this is
 * only implemented for N = 3.
 */
template <int N, typename T = Real> class /* THEA_API */ TorusN;

/**
 * A solid 3D torus, traced by a circular disk whose radius is the inner radius of the torus, and whose center moves along a
 * larger circle whose radius is the outer radius of the torus.
 */
template <typename T>
class /* THEA_API */ TorusN<3, T> : public RayIntersectableN<3, T>
{
  private:
    typedef AffineTransformN<3, T> AffineTransformT;  ///< Compatible affine transform type.

  public:
    THEA_DECL_SMART_POINTERS(TorusN)

    typedef Vector<3, T> VectorT;  ///< Compatible vector type.

    /** Default constructor. Does not initialize anything. */
    TorusN() {}

    /**
     * Initialize with the center point, the axis (unnormalized direction vector orthogonal to the plane of the torus), and
     * outer and inner radii.
     */
    TorusN(VectorT const & center_, VectorT const & axis_, T outer_radius_, T inner_radius_)
    : center(center_), axis(axis_.stableNormalized()), outer_radius(outer_radius_), inner_radius(inner_radius_)
    {
      updateTransform();
    }

    /** Get the center of the torus. */
    VectorT const & getCenter() const { return center; }

    /** Set the center of the torus. */
    void setCenter(VectorT const & center_) { center = center_; updateTransform(); }

    /** Get the axis of the torus, normalized to unit length. */
    VectorT const & getAxis() const { return axis; }

    /** Set the axis of the torus from an unnormalized direction vector. */
    void setAxis(VectorT const & axis_) { axis = axis_.stableNormalized(); updateTransform(); }

    /** Get the outer radius of the torus. */
    T getOuterRadius() const { return outer_radius; }

    /** Set the outer radius of the torus. */
    void setOuterRadius(T outer_radius_) { outer_radius = outer_radius_; }

    /** Get the inner radius of the torus. */
    T getInnerRadius() const { return inner_radius; }

    /** Set the inner radius of the torus. */
    void setInnerRadius(T inner_radius_) { inner_radius = inner_radius_; }

    /** Test if this torus intersects (contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the torus contains a point. */
    bool contains(VectorT const & p) const
    {
      VectorT local_p = world2local * p;
      T a2 = Math::square(outer_radius);
      T b2 = Math::square(inner_radius);
      return Math::square(local_p.squaredNorm() + a2 - b2) < 4 * a2 * (local_p.template head<2>().squaredNorm());
    }

    /** Get a bounding box for the torus. */
    AxisAlignedBoxN<3, T> getBounds() const
    {
      T s = outer_radius + inner_radius;
      return AxisAlignedBoxN<3, T>(VectorT(-s, -s, -inner_radius),
                                   VectorT( s,  s,  inner_radius)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the torus. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[center: " << Thea::toString(center) << ", axis: " << Thea::toString(axis)
          << ", outer_radius: " << outer_radius << ", inner_radius: " << inner_radius << ']';
      return oss.str();
    }

    bool rayIntersects(RayN<3, T> const & ray, T max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time) >= 0;
    }

    T rayIntersectionTime(RayN<3, T> const & ray, T max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time, nullptr);
    }

    RayIntersectionN<3, T> rayIntersection(RayN<3, T> const & ray, T max_time = -1) const
    {
      VectorT n;
      T t = rayIntersectionTime(ray, max_time, &n);
      if (t >= 0)
        return RayIntersectionN<3, T>(t, &n);
      else
        return RayIntersectionN<3, T>(-1);
    }

  private:
    /** Compute the transform that maps world coordinates to local torus coordinates. */
    void updateTransform()
    {
      // Transform is orthonormal -- we exploit this in other places
      Matrix<3, 3, T> lin = Math::rotationArc(axis, VectorT(0, 0, 1), /* normalize_dirs = */ false);
      world2local = AffineTransformT(lin, -(lin * center));
    }

    /**
     * Helper function for rayIntersectionTime(RayN<3, T> const &, T) that also optionally returns the surface normal at the hit
     * point.
     */
    T rayIntersectionTime(RayN<3, T> const & ray, T max_time, VectorT * normal) const
    {
      if (contains(ray.getOrigin()))
      {
        if (normal) *normal = VectorT::Zero();
        return 0;
      }

      auto local_ray = ray.transform(world2local);
      double min_root = -1;

      double r2 = Math::square((double)outer_radius);
      double w2 = Math::square((double)inner_radius);
      double r2pw2 = r2 + w2;
      double r2mw2 = r2 - w2;

      auto lrd = local_ray.template cast<double>();
      Vector<3, double> p2 = lrd.getOrigin().cwiseProduct(lrd.getOrigin());
      Vector<3, double> pu = lrd.getOrigin().cwiseProduct(lrd.getDirection());
      Vector<3, double> u2 = lrd.getDirection().cwiseProduct(lrd.getDirection());

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
        VectorT local_p = local_ray.getPoint((T)min_root);
        Vector2 ring_point = outer_radius * local_p.template head<2>().stableNormalized();
        *normal = (world2local.getLinear() * (local_p - VectorT(ring_point.x(), ring_point.y(), 0))).stableNormalized();
      }

      return (T)min_root;
    }

    VectorT           center;         ///< The center of the torus.
    VectorT           axis;           ///< Unit vector along the axis of the torus (perpendicular to main disk).
    T                 outer_radius;   ///< Outer radius of the torus.
    T                 inner_radius;   ///< Inner radius of the torus.
    AffineTransformT  world2local;    ///< Transform from world coordinates to local torus coordinates (axis along +Z).

}; // class TorusN<3, T>

} // namespace Thea

#include "Torus3.hpp"

#endif
