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

#ifndef __Thea_CylinderN_hpp__
#define __Thea_CylinderN_hpp__

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
 * A solid cylinder in N dimensions, defined as the Cartesian product of an (N - 1)-dimensional ball and a line segment.
 * Currently, this is only implemented for N = 3.
 */
template <int N, typename T = Real> class /* THEA_API */ CylinderN;

/** A solid 3D cylinder, traced by a circular disk whose center moves along an orthogonal line segment. */
template <typename T>
class /* THEA_API */ CylinderN<3, T> : public RayIntersectableN<3, T>
{
  private:
    typedef LineSegmentN<3, T>      LineSegmentT;      ///< Compatible line segment type.
    typedef AffineTransformN<3, T>  AffineTransformT;  ///< Compatible affine transform type.

  public:
    THEA_DECL_SMART_POINTERS(CylinderN)

    typedef Vector<3, T> VectorT;  ///< Compatible vector type.

    /** Default constructor. Does not initialize anything. */
    CylinderN() {}

    /** Initialize with an axis (line segment from the center of one end to the center of the other end) and a radius. */
    CylinderN(VectorT const & begin_, VectorT const & end_, T radius_) : segment(begin_, end_), radius(radius_)
    {
      updateTransform();
    }

    /** Get the center of one capping disk of the cylinder: 0 returns the first endpoint and 1 returns the second. */
    VectorT getEndpoint(int i) const { return segment.getEndpoint(i); }

    /** Get the unnormalized direction vector of the cylinder from the center of one end to the center of the other. */
    VectorT const & getDirection() const { return segment.getDirection(); }

    /** Get the length of the cylinder. */
    T length() const { return len; }

    /** Get the square of the length of the cylinder. */
    T squaredLength() const { return len * len; }

    /** Set the axis of the cylinder. */
    void setAxis(VectorT const & begin_, VectorT const & end_)
    {
      segment = LineSegmentT(begin_, end_);
      updateTransform();
    }

    /** Get the radius of the cylinder. */
    T getRadius() const { return radius; }

    /** Set the radius of the cylinder. */
    void setRadius(T radius_) { radius = radius_; updateTransform(); }

    /** Test if this cylinder intersects (contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the cylinder contains a point. */
    bool contains(VectorT const & p) const
    {
      VectorT local_p = world2local * p;
      return local_p.z() >= 0 && local_p.z() <= 1 && local_p.squaredNorm() <= 1;
    }

    /** Get a bounding box for the cylinder. */
    AxisAlignedBoxN<3, T> getBounds() const
    {
      return AxisAlignedBoxN<3, T>(VectorT(-radius, -radius, 0),
                                   VectorT( radius,  radius, len)).transformAndBound(world2local.inverse());
    }

    /** Get a textual representation of the cylinder. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[axis: " << segment.toString() << ", radius: " << radius << ']';
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
          case 1: n =  unit_dir; break;
          case 2:
          {
            VectorT p = ray.getPoint(t);
            n = (p - segment.closestPoint(p)) / radius;
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
    /** Compute the transform that maps world coordinates to normalized cylinder coordinates. */
    void updateTransform()
    {
      len = segment.length();
      if (len > 0 && radius > 0)
      {
        unit_dir = getDirection().stableNormalized();
        Matrix<3, 3, T> lin = Math::scaling(VectorT(1 / radius, 1 / radius, 1 / len))
                            * Math::rotationArc(unit_dir, VectorT(0, 0, 1), /* normalize_dirs = */ false);
        world2local = AffineTransformT(lin, -(lin * getEndpoint(0)));
      }
      else
      {
        unit_dir = VectorT::Zero();
        world2local.setIdentity();
      }
    }

    /**
     * Helper function for rayIntersectionTime(RayN<3, T> const &, T) that also returns the hit surface: 0 for bottom, 1 for
     * top, 2 for sides, 3 for interior, negative for no hit.
     */
    T rayIntersectionTime(RayN<3, T> const & ray, T max_time, int & surface) const
    {
      if (contains(ray.getOrigin()))
      {
        surface = 3;
        return 0;
      }

      surface = -1;
      auto local_ray = ray.transform(world2local);
      double min_root = -1;

      // Bottom and top
      for (int i = 0; i <= 1; ++i)
      {
        auto plane = HyperplaneN<3, T>::fromPointAndNormal(VectorT(0, 0, i), VectorT::UnitZ());
        auto plane_time = plane.rayIntersectionTime(local_ray, max_time);
        if (plane_time >= 0 && (min_root < 0 || plane_time < min_root))
        {
          VectorT p = local_ray.getPoint(plane_time);
          if (p.head(2).squaredNorm() <= 1)
          {
            min_root = plane_time;
            surface = i;
          }
        }
      }

      // Sides
      auto lrd = local_ray.template cast<double>();
      Vector<3, double> p2 = lrd.getOrigin().cwiseProduct(lrd.getOrigin());
      Vector<3, double> pu = lrd.getOrigin().cwiseProduct(lrd.getDirection());
      Vector<3, double> u2 = lrd.getDirection().cwiseProduct(lrd.getDirection());

      double s[3];
      s[2] = u2[0] + u2[1];
      s[1] = 2 * (pu[0] + pu[1]);
      s[0] = p2[0] + p2[1] - 1;

      double roots[2];
      int num_roots = Math::solveQuadratic(s[0], s[1], s[2], roots);

      for (int i = 0; i < num_roots; ++i)
        if (roots[i] >= 0 && (max_time < 0 || roots[i] <= max_time) && (min_root < 0 || roots[i] < min_root))
        {
          VectorT point = local_ray.getPoint((T)roots[i]);
          if (point[2] >= 0 && point[2] <= 1)
          {
            min_root = roots[i];
            surface = 2;
          }
        }

      return (T)min_root;
    }

    LineSegmentT      segment;      ///< Line segment defining the axis of the cylinder.
    T                 radius;       ///< Radius of the cylinder.
    VectorT           unit_dir;     ///< Unit vector along the axis of the cylinder.
    T                 len;          ///< Length of the cylinder.
    AffineTransformT  world2local;  ///< Transform from world coordinates to normalized cylinder coordinates.

}; // class CylinderN<3, T>

} // namespace Thea

#include "Cylinder3.hpp"

#endif
