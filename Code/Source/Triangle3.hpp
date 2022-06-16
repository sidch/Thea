//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the THEA_TRI3_BSD license, as detailed in the
// accompanying THEA_TRI3_LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// THEA_TRI3_LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
//
//============== License for line-triangle distance calculation ==============
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/THEA_TRI3_LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/THEA_TRI3_LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
//
//============================================================================

#ifndef __Thea_Triangle3_hpp__
#define __Thea_Triangle3_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include "BallN.hpp"
#include "BoxN.hpp"
#include "HyperplaneN.hpp"
#include "LineSegmentN.hpp"
#include "Math.hpp"
#include "RayIntersectableN.hpp"
#include "TriangleN.hpp"
#include <limits>

namespace Thea {

// Internal functions
namespace Triangle3Internal {

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * int tri_tri_intersect(float const V0[3],float const V1[3],float const V2[3],
 *                       float const U0[3],float const U1[3],float const U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDivTriTriIsect(float const V0[3],float const V1[3],float const V2[3],
 *                      float const U0[3],float const U1[3],float const U2[3]);
 *
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(float const V0[3],float const V1[3],float const V2[3],
 *                                      float const U0[3],float const U1[3],float const U2[3],int *coplanar,
 *                                      float isectpt1[3],float isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */

template <typename T>
int tri_tri_intersect(T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3]);

template <typename T>
int NoDivTriTriIsect(T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3]);

template <typename T>
int tri_tri_intersect_with_isectline(T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3],
                                     int * coplanar, T isectpt1[3], T isectpt2[3]);

// Check if a point is inside a triangle.
template <typename T>
bool isPointInsideTriangle(Vector<3, T> const & v0, Vector<3, T> const & v1, Vector<3, T> const & v2, int primary_axis,
                           Vector<3, T> const &  p);

// Closest point on the perimeter of a triangle.
template <typename T>
Vector<3, T> closestPointOnTrianglePerimeter(Vector<3, T> const & v0, Vector<3, T> const & v1, Vector<3, T> const & v2,
                                             Vector<3, T> const & point);

// Intersection time of a ray with a triangle. Returns a negative value if the ray does not intersect the triangle.
template <typename T>
T rayTriangleIntersectionTime(RayN<3, T> const & ray, Vector<3, T> const & v0,
                              Vector<3, T> const & edge01, Vector<3, T> const & edge02);

} // namespace Triangle3Internal

/**
 * A triangle (convex hull of 3 points) in 3-space, with precomputed properties for fast access. To account for the fact that
 * the triangle vertices may be stored in different ways (e.g. internally within the object, or as indices into an external
 * vertex pool), the class is parametrized on the way the vertices are stored. <code>VertexTripleT</code> must be
 * default-constructible, copy-constructible and assignable, and provide an efficient member function with the signature:
 * \code
 * VectorT const & getVertex(int i) const
 * \endcode
 */
template <typename VertexTripleT, typename T>
class /* THEA_DLL_LOCAL */ TriangleN<3, VertexTripleT, T>
: public TriangleNBase<3, VertexTripleT, T>, public RayIntersectableN<3, T>
{
  private:
    typedef TriangleNBase<3, VertexTripleT, T> BaseT;  ///< Base class.

  public:
    typedef typename BaseT::VectorT       VectorT;       ///< A vector in N-D space.
    typedef typename BaseT::VertexTriple  VertexTriple;  ///< Stores and provides access to the triangle's vertices.

    /** Default constructor. Does not initialize anything. */
    TriangleN() {}

    /** Construct from a set of vertices. */
    explicit TriangleN(VertexTripleT const & vertices_) : BaseT(vertices_) { update(); }

    /** Initializing constructor, callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or derived from it. */
    TriangleN(VectorT const & v0, VectorT const & v1, VectorT const & v2) : BaseT(v0, v1, v2) { update(); }

    /** Initialize the triangle from its vertex triple. */
    void set(VertexTriple const & vertices_) { BaseT::set(vertices_); update(); }

    /**
     * Set the i'th vertex. This function is callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or derived from
     * it.
     */
    void setVertex(int i, VectorT const & v) { BaseT::set(i, v); update(); }

    /**
     * Set all three vertices at once. This function is callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or
     * derived from it.
     */
    void set(VectorT const & v0, VectorT const & v1, VectorT const & v2) { BaseT::set(v0, v1, v2); update(); }

    /**
     * Update the properties of the triangle, assuming the positions of its three corners have changed. Useful for external
     * callers if the mutable positions are not locally stored within the triangle.
     */
    void update() const
    {
      VectorT v0 = BaseT::getVertex(0), v1 = BaseT::getVertex(1), v2 = BaseT::getVertex(2);

      plane = HyperplaneN<3, T>::fromThreePoints(v0, v1, v2);
      primary_axis = (int)Math::maxAbsAxis(plane.getNormal());

      centroid = (v0 + v1 + v2) / 3;
      edge01   = v1 - v0;
      edge02   = v2 - v0;

      area = (T)0.5 * std::fabs(plane.getNormal().dot(edge01.cross(edge02)));  // exploit precomputed normal to avoid sqrt
    }

    /** Get the plane of the triangle. */
    HyperplaneN<3, T> const & getPlane() const { return plane; }

    /** Get the primary axis of the triangle (closest to normal). */
    intx getPrimaryAxis() const { return primary_axis; }

    /** Get the normal of the triangle (right-hand rule, going round vertices in order 0, 1, 2). */
    VectorT const & getNormal() const { return plane.getNormal(); }

    /** Get the centroid of the triangle. */
    VectorT const & getCentroid() const { return centroid; }

    /** Get the edge vector corresponding to getVertex(1) - getVertex(0). */
    VectorT const & getEdge01() const { return edge01; }

    /** Get the edge vector corresponding to getVertex(2) - getVertex(0). */
    VectorT const & getEdge02() const { return edge02; }

    /** Get the area of the triangle. */
    T getArea() const { return area; }

    /** Get a uniformly distributed random sample from the triangle. */
    VectorT randomPoint() const
    {
      // From G3D::TriangleN

      // Choose a random point in the parallelogram
      T s = (T)Random::common().uniform01();
      T t = (T)Random::common().uniform01();

      if (s + t > 1)
      {
        // Outside the triangle; reflect about the diagonal of the parallelogram
        t = 1 - t;
        s = 1 - s;
      }

      return s * edge01 + t * edge02 + BaseT::getVertex(0);
    }

    /** Get the barycentric coordinates of a point assumed to be in the plane of the triangle. */
    VectorT barycentricCoordinates(VectorT const & p) const
    {
      if (area <= 0)
        return VectorT::Zero();

      VectorT const & n = getNormal();
      VectorT p0 = BaseT::getVertex(0) - p;
      VectorT p1 = BaseT::getVertex(1) - p;
      VectorT p2 = BaseT::getVertex(2) - p;

      T b0 = n.dot(p1.cross(p2)) / (2 * area);
      T b1 = n.dot(p2.cross(p0)) / (2 * area);

      return VectorT(b0, b1, 1 - b0 - b1);
    }

    /** Get a bounding box for the triangle. */
    AxisAlignedBoxN<3, T> getBounds() const
    {
      VectorT v0 = BaseT::getVertex(0), v1 = BaseT::getVertex(1), v2 = BaseT::getVertex(2);
      return AxisAlignedBoxN<3, T>(v0.cwiseMin(v1.cwiseMin(v2)), v0.cwiseMax(v1.cwiseMax(v2)));
    }

    /** Check if the triangle intersects (that is, contains) a point. */
    bool intersects(VectorT const & p) const { return contains(p); }

    /** Check if the triangle intersects another triangle. */
    template <typename OtherVertexTripleT> bool intersects(TriangleN<3, OtherVertexTripleT, T> const & other) const
    {
      VectorT p0 = BaseT::getVertex(0);  VectorT p1 = BaseT::getVertex(1);  VectorT p2 = BaseT::getVertex(2);
      VectorT q0 = other.getVertex(0);   VectorT q1 = other.getVertex(1);   VectorT q2 = other.getVertex(2);

      T v0[3] = { p0.x(), p0.y(), p0.z() };
      T v1[3] = { p1.x(), p1.y(), p1.z() };
      T v2[3] = { p2.x(), p2.y(), p2.z() };

      T u0[3] = { q0.x(), q0.y(), q0.z() };
      T u1[3] = { q1.x(), q1.y(), q1.z() };
      T u2[3] = { q2.x(), q2.y(), q2.z() };

      return Triangle3Internal::NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
    }

    /**
     * Check if the triangle intersects another triangle. If they do intersect, check if they are coplanar. If they are not
     * coplanar, compute the line of intersection. The intersection test is somewhat slower than intersects().
     */
    template <typename OtherVertexTripleT>
    bool intersects(TriangleN<3, OtherVertexTripleT, T> const & other, bool & coplanar, LineSegmentN<3, T> & seg) const
    {
      VectorT p0 = BaseT::getVertex(0);  VectorT p1 = BaseT::getVertex(1);  VectorT p2 = BaseT::getVertex(2);
      VectorT q0 = other.getVertex(0);   VectorT q1 = other.getVertex(1);   VectorT q2 = other.getVertex(2);

      int i_coplanar;
      VectorT isectpt1, isectpt2;
      int isec = Triangle3Internal::tri_tri_intersect_with_isectline(&p0[0], &p1[0], &p2[0], &q0[0], &q1[0], &q2[0],
                                                                     &i_coplanar, &isectpt1[0], &isectpt2[0]);
      if (isec)
      {
        coplanar = (bool)i_coplanar;
        if (!coplanar)
          seg = LineSegmentN<3, T>(isectpt1, isectpt2);

        return true;
      }

      return false;
    }

    /** Check if the triangle intersects a ball. */
    bool intersects(BallN<3, T> const & ball) const { throw Error("TriangleN<3>: Intersection with ball not implemented"); }

    /** Check if the triangle intersects an axis-aligned box. */
    bool intersects(AxisAlignedBoxN<3, T> const & aab) const
    { throw Error("TriangleN<3>: Intersection with AAB not implemented"); }

    /** Check if the triangle intersects an oriented box. */
    bool intersects(Box3 const & box) const { throw Error("TriangleN<3>: Intersection with oriented box not implemented"); }

    /** Check if the triangle contains a point. */
    bool contains(VectorT const & p) const
    {
      return Triangle3Internal::isPointInsideTriangle(BaseT::getVertex(0), BaseT::getVertex(1), BaseT::getVertex(2),
                                                      primary_axis, p);
    }

    /** Get the distance of the triangle from a point. */
    T distance(VectorT const & p) const { return std::sqrt(squaredDistance(p)); }

    /** Get the squared distance of the triangle from a point. */
    T squaredDistance(VectorT const & p) const
    {
      return (closestPoint(p) - p).squaredNorm();
    }

    /** Get the point on this triangle closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      // Project the point onto the plane of the triangle
      VectorT proj = plane.closestPoint(p);

      if (contains(proj))
        return proj;
      else  // the closest point is on the perimeter instead
        return Triangle3Internal::closestPointOnTrianglePerimeter(BaseT::getVertex(0), BaseT::getVertex(1), BaseT::getVertex(2),
                                                                  p);
    }

    /** Get the distance of the triangle from another triangle. */
    template <typename OtherVertexTripleT>
    T distance(TriangleN<3, OtherVertexTripleT, T> const & other) const { return std::sqrt(squaredDistance(other)); }

    /**
     * Get the squared distance between this triangle and another triangle, and optionally return the closest pair of points.
     */
    template <typename OtherVertexTripleT>
    T squaredDistance(TriangleN<3, OtherVertexTripleT, T> const & other, VectorT * this_pt = nullptr,
                      VectorT * other_pt = nullptr) const
    {
      // From Christer Ericson, "T-Time Collision Detection", Morgan-Kaufman, 2005.

      VectorT c1(0, 0, 0), c2(0, 0, 0);  // squash an uninitialized variable warning
      T min_sqdist = std::numeric_limits<T>::max();

      // First test for intersection
      bool coplanar;
      LineSegmentN<3, T> seg;
      if (intersects(other, coplanar, seg))
      {
        if (coplanar)
          c1 = c2 = (T)0.5 * (getCentroid() + other.getCentroid());  // THEA_TRI3_FIXME: This needn't be in the intersection
        else
          c1 = c2 = (T)0.5 * (seg.getPoint(0) + seg.getPoint(1));

        min_sqdist = 0;
      }
      else
      {
        // Edge-edge distances
        VectorT p, q;
        T d2, s, t;
        for (int i = 0; i < 3; ++i)
        {
          int i2 = (i + 1) % 3;

          for (int j = 0; j < 3; ++j)
          {
            int j2 = (j + 1) % 3;
            d2 = Internal::closestPtSegmentSegment<3, T>(BaseT::getVertex(i), BaseT::getVertex(i2), false,
                                                         other.getVertex(j), other.getVertex(j2), false,
                                                         s, t, p, q);
            if (d2 < min_sqdist)
            {
              min_sqdist = d2;
              c1 = p;
              c2 = q;
            }
          }
        }

        // Distance from vertex of triangle 2 to triangle 1, if the former projects inside the latter
        for (int i = 0; i < 3; ++i)
        {
          q = other.getVertex(i);
          p = getPlane().closestPoint(q);
          if (contains(p))
          {
            d2 = (p - q).squaredNorm();
            if (d2 < min_sqdist)
            {
              min_sqdist = d2;
              c1 = p;
              c2 = q;
            }
          }
        }

        // Distance from vertex of triangle 1 to triangle 2, if the former projects inside the latter
        for (int i = 0; i < 3; ++i)
        {
          p = BaseT::getVertex(i);
          q = other.getPlane().closestPoint(p);
          if (other.contains(q))
          {
            d2 = (p - q).squaredNorm();
            if (d2 < min_sqdist)
            {
              min_sqdist = d2;
              c1 = p;
              c2 = q;
            }
          }
        }
      }

      if (this_pt)  *this_pt  = c1;
      if (other_pt) *other_pt = c2;

      return min_sqdist;
    }

    /** Get the distance of this triangle from a ball. */
    T distance(BallN<3, T> const & ball) const
    {
      return std::max(distance(ball.getCenter()) - ball.getRadius(), static_cast<T>(0));
    }

    /** Get the squared distance between this triangle and a ball, and optionally return the closest pair of points. */
    T squaredDistance(BallN<3, T> const & ball, VectorT * this_pt = nullptr, VectorT * ball_pt = nullptr) const
    {
      if (!this_pt && !ball_pt)
      {
        T x = distance(ball);
        return x * x;
      }

      VectorT c1 = closestPoint(ball.getCenter());
      VectorT c2;

      VectorT diff = c1 - ball.getCenter();
      T d2 = diff.squaredNorm();
      T r2 = ball.getRadius() * ball.getRadius();
      if (d2 < r2)  // point inside ball
      {
        c2 = c1;
        d2 = 0;
      }
      else
      {
        if (r2 <= std::numeric_limits<T>::min())
          c2 = ball.getCenter();
        else
        {
          c2 = ball.getCenter() + std::sqrt(r2 / d2) * diff;
          d2 = (c1 - c2).squaredNorm();
        }
      }

      if (this_pt) *this_pt = c1;
      if (ball_pt) *ball_pt = c2;

      return d2;
    }

    /** Get the distance of this triangle from an infinite line. */
    T distance(LineN<3, T> const & line) const
    {
      return std::sqrt(squaredDistance(line));
    }

    /**
     * Get the squared distance between this triangle and an infinite line, and optionally return the closest pair of points.
     */
    T squaredDistance(LineN<3, T> const & line, VectorT * this_pt = nullptr, VectorT * line_pt = nullptr) const
    {
      // Two tests where one would suffice, but for now it avoids having to code a new function or modify an existing one.
      // Shift the ray origins to ensure the rays overlap and there are no errors at line.getPoint().

      // Forward half
      {
        RayN<3, T> ray(line.getPoint() - line.getDirection(), line.getDirection());
        T t = rayIntersectionTime(ray);
        if (t >= 0)
        {
          VectorT c = ray.getPoint(t);
          if (this_pt) *this_pt = c;
          if (line_pt) *line_pt = c;
          return 0;
        }
      }

      // Backward half
      {
        RayN<3, T> ray(line.getPoint() + line.getDirection(), -line.getDirection());
        T t = rayIntersectionTime(ray);
        if (t >= 0)
        {
          VectorT c = ray.getPoint(t);
          if (this_pt) *this_pt = c;
          if (line_pt) *line_pt = c;
          return 0;
        }
      }

      return squaredDistanceToPerimeter(line, this_pt, line_pt);
    }

    /** Get the distance of this triangle from a line segment. */
    T distance(LineSegmentN<3, T> const & seg) const
    {
      return std::sqrt(squaredDistance(seg));
    }

    /**
     * Get the squared distance between this triangle and a line segment, and optionally return the closest pair of points.
     */
    T squaredDistance(LineSegmentN<3, T> const & seg, VectorT * this_pt = nullptr, VectorT * seg_pt = nullptr) const
    {
      // From https://www.geometrictools.com/THEA_TRI3_GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      RayN<3, T> ray(seg.getEndpoint(0), seg.getDirection());
      T t = rayIntersectionTime(ray);
      if (t >= 0 && t <= 1)
      {
        VectorT c = ray.getPoint(t);
        if (this_pt) *this_pt = c;
        if (seg_pt)  *seg_pt  = c;
        return 0;
      }

      return squaredDistanceToPerimeter(seg, this_pt, seg_pt);
    }

    /** Get the distance of this triangle from a ray. */
    T distance(RayN<3, T> const & ray) const
    {
      return std::sqrt(squaredDistance(ray));
    }

    /**
     * Get the squared distance between this triangle and a ray, and optionally return the closest pair of points.
     */
    T squaredDistance(RayN<3, T> const & ray, VectorT * this_pt = nullptr, VectorT * ray_pt = nullptr) const
    {
      // From https://www.geometrictools.com/THEA_TRI3_GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      T t = rayIntersectionTime(ray);
      if (t >= 0)
      {
        VectorT c = ray.getPoint(t);
        if (this_pt) *this_pt = c;
        if (ray_pt)  *ray_pt  = c;
        return 0;
      }

      return squaredDistanceToPerimeter(ray, this_pt, ray_pt);
    }

    T rayIntersectionTime(RayN<3, T> const & ray, T max_time = -1) const
    {
      T t = Triangle3Internal::rayTriangleIntersectionTime(ray, BaseT::getVertex(0), getEdge01(), getEdge02());
      return (max_time >= 0 && t > max_time) ? -1 : t;
    }

    RayIntersectionN<3, T> rayIntersection(RayN<3, T> const & ray, T max_time = -1) const
    {
      T t = Triangle3Internal::rayTriangleIntersectionTime(ray, BaseT::getVertex(0), getEdge01(), getEdge02());
      if (t >= 0 && (max_time < 0 || t <= max_time))
      {
        VectorT n = getNormal();
        return RayIntersectionN<3, T>(t, &n);
      }

      return RayIntersectionN<3, T>(-1);
    }

  protected:
    /**
     * Get the squared distance between the perimeter of this triangle and a line-like object (LineN, RayN or LineSegmentN), and
     * optionally return the closest pair of points.
     */
    template <typename LineLikeT>
    T squaredDistanceToPerimeter(LineLikeT const & line, VectorT * this_pt = nullptr, VectorT * line_pt = nullptr) const
    {
      // From https://www.geometrictools.com/THEA_TRI3_GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      T d2 = -1;
      VectorT c1, c2;
      for (int i = 0; i < 3; ++i)
      {
        LineSegmentN<3, T> edge(BaseT::getVertex(i), BaseT::getVertex((i + 1) % 3));
        T edge_d2 = edge.squaredDistance(line, (this_pt ? &c1 : nullptr), (line_pt ? &c2 : nullptr));
        if (d2 < 0 || edge_d2 < d2)
        {
          d2 = edge_d2;
          if (this_pt) *this_pt = c1;
          if (line_pt) *line_pt = c2;
        }
      }

      return d2;
    }

    // Cached properties computed from the vertices
    mutable HyperplaneN<3, T>  plane;         ///< Plane of the triangle.
    mutable int                primary_axis;  ///< Primary axis (closest to normal).
    mutable VectorT            centroid;      ///< Centroid of the triangle (mean of three vertices).
    mutable VectorT            edge01;        ///< vertices[1] - vertices[0]
    mutable VectorT            edge02;        ///< vertices[2] - vertices[0]
    mutable T                  area;          ///< TriangleN area.

}; // class TriangleN<3, VertexTripleT, T>

/**
 * The default triangle class in 3-dimensional real space, whose vertex positions are stored locally in the TriangleN object
 * itself.
 */
typedef LocalTriangleN<3, Real> Triangle3;

namespace Triangle3Internal {

/* if THEA_TRI3_USE_EPSILON_TEST is true then we do a check:
         if |dv|<THEA_TRI3_EPSILON then dv=0;
   else no check is done (which is less robust)
*/
#define THEA_TRI3_USE_EPSILON_TEST true
#define THEA_TRI3_EPSILON Math::eps<T>()

/* some macros */
#define THEA_TRI3_CROSS(dest,v1,v2) \
  dest[0]=v1[1]*v2[2]-v1[2]*v2[1];  \
  dest[1]=v1[2]*v2[0]-v1[0]*v2[2];  \
  dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define THEA_TRI3_DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define THEA_TRI3_SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2];

#define THEA_TRI3_ADD(dest,v1,v2) dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2];

#define THEA_TRI3_MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2];

#define THEA_TRI3_SET(dest,src) dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2];

/* sort so that a<=b */
#define THEA_TRI3_SORT(a,b) \
  if(a>b)    \
  {          \
    T c;     \
    c=a;     \
    a=b;     \
    b=c;     \
  }

#define THEA_TRI3_ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
  isect0=VV0+(VV1-VV0)*D0/(D0-D1);                          \
  isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define THEA_TRI3_COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0)                                                                      \
  {                                                                               \
    /* here we know that D0D2<=0 */                                               \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */    \
    THEA_TRI3_ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);                          \
  }                                                                               \
  else if(D0D2>0)                                                                 \
  {                                                                               \
    /* here we know that d0d1<=0 */                                               \
    THEA_TRI3_ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);                          \
  }                                                                               \
  else if(D1*D2>0 || D0!=0)                                                       \
  {                                                                               \
    /* here we know that d0d1<=0 or that D0!=0 */                                 \
    THEA_TRI3_ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);                          \
  }                                                                               \
  else if(D1!=0)                                                                  \
  {                                                                               \
    THEA_TRI3_ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);                          \
  }                                                                               \
  else if(D2!=0)                                                                  \
  {                                                                               \
    THEA_TRI3_ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);                          \
  }                                                                               \
  else                                                                            \
  {                                                                               \
    /* triangles are coplanar */                                                  \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);                                \
  }

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems THEA_TRI3_III,
   pp. 199-202 */
#define THEA_TRI3_EDGE_EDGE_TEST(V0,U0,U1)            \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }

#define THEA_TRI3_EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2)   \
  {                                                        \
    T Ax,Ay,Bx,By,Cx,Cy,e,d,f;                             \
    Ax=V1[i0]-V0[i0];                                      \
    Ay=V1[i1]-V0[i1];                                      \
    /* test edge U0,U1 against V0,V1 */                    \
    THEA_TRI3_EDGE_EDGE_TEST(V0,U0,U1);                    \
    /* test edge U1,U2 against V0,V1 */                    \
    THEA_TRI3_EDGE_EDGE_TEST(V0,U1,U2);                    \
    /* test edge U2,U1 against V0,V1 */                    \
    THEA_TRI3_EDGE_EDGE_TEST(V0,U2,U0);                    \
  }

#define THEA_TRI3_POINT_IN_TRI(V0,U0,U1,U2)   \
  {                                           \
    T a,b,c,d0,d1,d2;                         \
    /* is T1 completly inside T2? */          \
    /* check if V0 is inside tri(U0,U1,U2) */ \
    a=U1[i1]-U0[i1];                          \
    b=-(U1[i0]-U0[i0]);                       \
    c=-a*U0[i0]-b*U0[i1];                     \
    d0=a*V0[i0]+b*V0[i1]+c;                   \
                                              \
    a=U2[i1]-U1[i1];                          \
    b=-(U2[i0]-U1[i0]);                       \
    c=-a*U1[i0]-b*U1[i1];                     \
    d1=a*V0[i0]+b*V0[i1]+c;                   \
                                              \
    a=U0[i1]-U2[i1];                          \
    b=-(U0[i0]-U2[i0]);                       \
    c=-a*U2[i0]-b*U2[i1];                     \
    d2=a*V0[i0]+b*V0[i1]+c;                   \
    if(d0*d1>0)                               \
    {                                         \
      if(d0*d2>0) return 1;                   \
    }                                         \
  }

template <typename T>
int coplanar_tri_tri(T const N[3], T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3])
{
  T A[3];
  short i0, i1;
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  A[0] = std::fabs(N[0]);
  A[1] = std::fabs(N[1]);
  A[2] = std::fabs(N[2]);

  if (A[0] > A[1])
  {
    if (A[0] > A[2])
    {
      i0 = 1;    /* A[0] is greatest */
      i1 = 2;
    }
    else
    {
      i0 = 0;    /* A[2] is greatest */
      i1 = 1;
    }
  }
  else   /* A[0]<=A[1] */
  {
    if (A[2] > A[1])
    {
      i0 = 0;    /* A[2] is greatest */
      i1 = 1;
    }
    else
    {
      i0 = 0;    /* A[1] is greatest */
      i1 = 2;
    }
  }

  /* test all edges of triangle 1 against the edges of triangle 2 */
  THEA_TRI3_EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
  THEA_TRI3_EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
  THEA_TRI3_EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);
  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  THEA_TRI3_POINT_IN_TRI(V0, U0, U1, U2);
  THEA_TRI3_POINT_IN_TRI(U0, V0, V1, V2);
  return 0;
}

template <typename T>
int tri_tri_intersect(T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3])
{
  T E1[3], E2[3];
  T N1[3], N2[3], d1, d2;
  T du0, du1, du2, dv0, dv1, dv2;
  T D[3];
  T isect1[2], isect2[2];
  T du0du1, du0du2, dv0dv1, dv0dv2;
  short index;
  T vp0, vp1, vp2;
  T up0, up1, up2;
  T b, c, max;
  /* compute plane equation of triangle(V0,V1,V2) */
  THEA_TRI3_SUB(E1, V1, V0);
  THEA_TRI3_SUB(E2, V2, V0);
  THEA_TRI3_CROSS(N1, E1, E2);
  d1 = -THEA_TRI3_DOT(N1, V0);
  /* plane equation 1: N1.X+d1=0 */
  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = THEA_TRI3_DOT(N1, U0) + d1;
  du1 = THEA_TRI3_DOT(N1, U1) + d1;
  du2 = THEA_TRI3_DOT(N1, U2) + d1;
  /* coplanarity robustness check */
#if THEA_TRI3_USE_EPSILON_TEST

  if (std::fabs(du0) < THEA_TRI3_EPSILON) du0 = 0;

  if (std::fabs(du1) < THEA_TRI3_EPSILON) du1 = 0;

  if (std::fabs(du2) < THEA_TRI3_EPSILON) du2 = 0;

#endif
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;

  if (du0du1 > 0 && du0du2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;                   /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  THEA_TRI3_SUB(E1, U1, U0);
  THEA_TRI3_SUB(E2, U2, U0);
  THEA_TRI3_CROSS(N2, E1, E2);
  d2 = -THEA_TRI3_DOT(N2, U0);
  /* plane equation 2: N2.X+d2=0 */
  /* put V0,V1,V2 into plane equation 2 */
  dv0 = THEA_TRI3_DOT(N2, V0) + d2;
  dv1 = THEA_TRI3_DOT(N2, V1) + d2;
  dv2 = THEA_TRI3_DOT(N2, V2) + d2;
#if THEA_TRI3_USE_EPSILON_TEST

  if (std::fabs(dv0) < THEA_TRI3_EPSILON) dv0 = 0;

  if (std::fabs(dv1) < THEA_TRI3_EPSILON) dv1 = 0;

  if (std::fabs(dv2) < THEA_TRI3_EPSILON) dv2 = 0;

#endif
  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;

  if (dv0dv1 > 0 && dv0dv2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;                   /* no intersection occurs */

  /* compute direction of intersection line */
  THEA_TRI3_CROSS(D, N1, N2);
  /* compute and index to the largest component of D */
  max = std::fabs(D[0]);
  index = 0;
  b = std::fabs(D[1]);
  c = std::fabs(D[2]);

  if (b > max) max = b, index = 1;

  if (c > max) max = c, index = 2;

  /* this is the simplified projection onto L*/
  vp0 = V0[index];
  vp1 = V1[index];
  vp2 = V2[index];
  up0 = U0[index];
  up1 = U1[index];
  up2 = U2[index];
  /* compute interval for triangle 1 */
  THEA_TRI3_COMPUTE_INTERVALS(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, isect1[0], isect1[1]);
  /* compute interval for triangle 2 */
  THEA_TRI3_COMPUTE_INTERVALS(up0, up1, up2, du0, du1, du2, du0du1, du0du2, isect2[0], isect2[1]);
  THEA_TRI3_SORT(isect1[0], isect1[1]);
  THEA_TRI3_SORT(isect2[0], isect2[1]);

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;

  return 1;
}

#define THEA_TRI3_NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
  { \
    if(D0D1>0) \
    { \
      /* here we know that D0D2<=0 */ \
      /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
      A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
    } \
    else if(D0D2>0)\
    { \
      /* here we know that d0d1<=0.0f */ \
      A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
    } \
    else if(D1*D2>0 || D0!=0) \
    { \
      /* here we know that d0d1<=0.0f or that D0!=0.0f */ \
      A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
    } \
    else if(D1!=0) \
    { \
      A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
    } \
    else if(D2!=0) \
    { \
      A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
    } \
    else \
    { \
      /* triangles are coplanar */ \
      return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
    } \
  }

template <typename T>
int NoDivTriTriIsect(T const V0[3], T const V1[3], T const V2[3], T const U0[3], T const U1[3], T const U2[3])
{
  T E1[3], E2[3];
  T N1[3], N2[3], d1, d2;
  T du0, du1, du2, dv0, dv1, dv2;
  T D[3];
  T isect1[2], isect2[2];
  T du0du1, du0du2, dv0dv1, dv0dv2;
  short index;
  T vp0, vp1, vp2;
  T up0, up1, up2;
  T bb, cc, max;
  T a, b, c, x0, x1;
  T d, e, f, y0, y1;
  T xx, yy, xxyy, tmp;
  /* compute plane equation of triangle(V0,V1,V2) */
  THEA_TRI3_SUB(E1, V1, V0);
  THEA_TRI3_SUB(E2, V2, V0);
  THEA_TRI3_CROSS(N1, E1, E2);
  d1 = -THEA_TRI3_DOT(N1, V0);
  /* plane equation 1: N1.X+d1=0 */
  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = THEA_TRI3_DOT(N1, U0) + d1;
  du1 = THEA_TRI3_DOT(N1, U1) + d1;
  du2 = THEA_TRI3_DOT(N1, U2) + d1;
  /* coplanarity robustness check */
#if THEA_TRI3_USE_EPSILON_TEST

  if (THEA_TRI3_FABS(du0) < THEA_TRI3_EPSILON) du0 = 0;

  if (THEA_TRI3_FABS(du1) < THEA_TRI3_EPSILON) du1 = 0;

  if (THEA_TRI3_FABS(du2) < THEA_TRI3_EPSILON) du2 = 0;

#endif
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;

  if (du0du1 > 0 && du0du2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;                   /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  THEA_TRI3_SUB(E1, U1, U0);
  THEA_TRI3_SUB(E2, U2, U0);
  THEA_TRI3_CROSS(N2, E1, E2);
  d2 = -THEA_TRI3_DOT(N2, U0);
  /* plane equation 2: N2.X+d2=0 */
  /* put V0,V1,V2 into plane equation 2 */
  dv0 = THEA_TRI3_DOT(N2, V0) + d2;
  dv1 = THEA_TRI3_DOT(N2, V1) + d2;
  dv2 = THEA_TRI3_DOT(N2, V2) + d2;
#if THEA_TRI3_USE_EPSILON_TEST

  if (THEA_TRI3_FABS(dv0) < THEA_TRI3_EPSILON) dv0 = 0;

  if (THEA_TRI3_FABS(dv1) < THEA_TRI3_EPSILON) dv1 = 0;

  if (THEA_TRI3_FABS(dv2) < THEA_TRI3_EPSILON) dv2 = 0;

#endif
  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;

  if (dv0dv1 > 0 && dv0dv2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;                   /* no intersection occurs */

  /* compute direction of intersection line */
  THEA_TRI3_CROSS(D, N1, N2);
  /* compute and index to the largest component of D */
  max = (T)THEA_TRI3_FABS(D[0]);
  index = 0;
  bb = (T)THEA_TRI3_FABS(D[1]);
  cc = (T)THEA_TRI3_FABS(D[2]);

  if (bb > max) max = bb, index = 1;

  if (cc > max) max = cc, index = 2;

  /* this is the simplified projection onto L*/
  vp0 = V0[index];
  vp1 = V1[index];
  vp2 = V2[index];
  up0 = U0[index];
  up1 = U1[index];
  up2 = U2[index];
  /* compute interval for triangle 1 */
  THEA_TRI3_NEWCOMPUTE_INTERVALS(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, a, b, c, x0, x1);
  /* compute interval for triangle 2 */
  THEA_TRI3_NEWCOMPUTE_INTERVALS(up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f, y0, y1);
  xx = x0 * x1;
  yy = y0 * y1;
  xxyy = xx * yy;
  tmp = a * xxyy;
  isect1[0] = tmp + b * x1 * yy;
  isect1[1] = tmp + c * x0 * yy;
  tmp = d * xxyy;
  isect2[0] = tmp + e * xx * y1;
  isect2[1] = tmp + f * xx * y0;
  THEA_TRI3_SORT(isect1[0], isect1[1]);
  THEA_TRI3_SORT(isect2[0], isect2[1]);

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;

  return 1;
}

/* sort so that a<=b */
#define THEA_TRI3_SORT2(a,b,smallest) \
  if(a>b)                             \
  {                                   \
    T c;                              \
    c=a;                              \
    a=b;                              \
    b=c;                              \
    smallest=1;                       \
  }                                   \
  else smallest=0;

template <typename T>
inline void isect2(T const VTX0[3], T const VTX1[3], T const VTX2[3], T VV0, T VV1, T VV2,
                   T D0, T D1, T D2, T * isect0, T * isect1, T isectpoint0[3], T isectpoint1[3])
{
  T tmp = D0 / (D0 - D1);
  T diff[3];
  *isect0 = VV0 + (VV1 - VV0) * tmp;
  THEA_TRI3_SUB(diff, VTX1, VTX0);
  THEA_TRI3_MULT(diff, diff, tmp);
  THEA_TRI3_ADD(isectpoint0, diff, VTX0);
  tmp = D0 / (D0 - D2);
  *isect1 = VV0 + (VV2 - VV0) * tmp;
  THEA_TRI3_SUB(diff, VTX2, VTX0);
  THEA_TRI3_MULT(diff, diff, tmp);
  THEA_TRI3_ADD(isectpoint1, VTX0, diff);
}

template <typename T>
inline int compute_intervals_isectline(T const VERT0[3], T const VERT1[3], T const VERT2[3],
                                       T VV0, T VV1, T VV2, T D0, T D1, T D2,
                                       T D0D1, T D0D2, T * isect0, T * isect1,
                                       T isectpoint0[3], T isectpoint1[3])
{
  if (D0D1 > 0)
  {
    /* here we know that D0D2<=0 */
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
  }
  else if (D0D2 > 0)
  {
    /* here we know that d0d1<=0 */
    isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
  }
  else if (D1 * D2 > 0 || D0 != 0)
  {
    /* here we know that d0d1<=0 or that D0!=0 */
    isect2(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, isect0, isect1, isectpoint0, isectpoint1);
  }
  else if (D1 != 0)
  {
    isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, isect0, isect1, isectpoint0, isectpoint1);
  }
  else if (D2 != 0)
  {
    isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, isect0, isect1, isectpoint0, isectpoint1);
  }
  else
  {
    /* triangles are coplanar */
    return 1;
  }

  return 0;
}

#define THEA_TRI3_COMPUTE_INTERVALS_ISECTLINE(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2, \
                                              isect0,isect1,isectpoint0,isectpoint1) \
  if(D0D1>0) \
  { \
    /* here we know that D0D2<=0 */ \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1); \
  }
#if 0
  else if (D0D2 > 0) \
  { \
    /* here we know that d0d1<=0 */ \
    isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, &isect0, &isect1, isectpoint0, isectpoint1); \
  } \
  else if (D1 * D2 > 0 || D0 != 0) \
  { \
    /* here we know that d0d1<=0 or that D0!=0 */ \
    isect2(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, &isect0, &isect1, isectpoint0, isectpoint1); \
  } \
  else if (D1 != 0) \
  { \
    isect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, &isect0, &isect1, isectpoint0, isectpoint1); \
  } \
  else if (D2 != 0) \
  { \
    isect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, &isect0, &isect1, isectpoint0, isectpoint1); \
  } \
  else \
  { \
    /* triangles are coplanar */ \
    coplanar = 1; \
    return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2); \
  }
#endif

template <typename T>
int tri_tri_intersect_with_isectline(T const V0[3], T const V1[3], T const V2[3],
                                     T const U0[3], T const U1[3], T const U2[3],
                                     int * coplanar, T isectpt1[3], T isectpt2[3])
{
  T E1[3], E2[3];
  T N1[3], N2[3], d1, d2;
  T du0, du1, du2, dv0, dv1, dv2;
  T D[3];
  T isect1[2] = {0}, isect2[2] = {0};
  T isectpointA1[3] = {0}, isectpointA2[3] = {0};
  T isectpointB1[3] = {0}, isectpointB2[3] = {0};
  T du0du1, du0du2, dv0dv1, dv0dv2;
  short index;
  T vp0, vp1, vp2;
  T up0, up1, up2;
  T b, c, max;
  // T tmp,diff[3];
  int smallest1, smallest2;
  /* compute plane equation of triangle(V0,V1,V2) */
  THEA_TRI3_SUB(E1, V1, V0);
  THEA_TRI3_SUB(E2, V2, V0);
  THEA_TRI3_CROSS(N1, E1, E2);
  d1 = -THEA_TRI3_DOT(N1, V0);
  /* plane equation 1: N1.X+d1=0 */
  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = THEA_TRI3_DOT(N1, U0) + d1;
  du1 = THEA_TRI3_DOT(N1, U1) + d1;
  du2 = THEA_TRI3_DOT(N1, U2) + d1;
  /* coplanarity robustness check */
#if THEA_TRI3_USE_EPSILON_TEST

  if (std::fabs(du0) < THEA_TRI3_EPSILON) du0 = 0;

  if (std::fabs(du1) < THEA_TRI3_EPSILON) du1 = 0;

  if (std::fabs(du2) < THEA_TRI3_EPSILON) du2 = 0;

#endif
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;

  if (du0du1 > 0 && du0du2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;              /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  THEA_TRI3_SUB(E1, U1, U0);
  THEA_TRI3_SUB(E2, U2, U0);
  THEA_TRI3_CROSS(N2, E1, E2);
  d2 = -THEA_TRI3_DOT(N2, U0);
  /* plane equation 2: N2.X+d2=0 */
  /* put V0,V1,V2 into plane equation 2 */
  dv0 = THEA_TRI3_DOT(N2, V0) + d2;
  dv1 = THEA_TRI3_DOT(N2, V1) + d2;
  dv2 = THEA_TRI3_DOT(N2, V2) + d2;
#if THEA_TRI3_USE_EPSILON_TEST

  if (std::fabs(dv0) < THEA_TRI3_EPSILON) dv0 = 0;

  if (std::fabs(dv1) < THEA_TRI3_EPSILON) dv1 = 0;

  if (std::fabs(dv2) < THEA_TRI3_EPSILON) dv2 = 0;

#endif
  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;

  if (dv0dv1 > 0 && dv0dv2 > 0) /* same sign on all of them + not equal 0 ? */
    return 0;              /* no intersection occurs */

  /* compute direction of intersection line */
  THEA_TRI3_CROSS(D, N1, N2);
  /* compute and index to the largest component of D */
  max = std::fabs(D[0]);
  index = 0;
  b = std::fabs(D[1]);
  c = std::fabs(D[2]);

  if (b > max) max = b, index = 1;

  if (c > max) max = c, index = 2;

  /* this is the simplified projection onto L*/
  vp0 = V0[index];
  vp1 = V1[index];
  vp2 = V2[index];
  up0 = U0[index];
  up1 = U1[index];
  up2 = U2[index];
  /* compute interval for triangle 1 */
  *coplanar = compute_intervals_isectline(V0, V1, V2, vp0, vp1, vp2, dv0, dv1, dv2,
                                          dv0dv1, dv0dv2, &isect1[0], &isect1[1], isectpointA1, isectpointA2);

  if (*coplanar) return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);

  /* compute interval for triangle 2 */
  compute_intervals_isectline(U0, U1, U2, up0, up1, up2, du0, du1, du2,
                              du0du1, du0du2, &isect2[0], &isect2[1], isectpointB1, isectpointB2);
  THEA_TRI3_SORT2(isect1[0], isect1[1], smallest1);
  THEA_TRI3_SORT2(isect2[0], isect2[1], smallest2);

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;

  /* at this point, we know that the triangles intersect */

  if (isect2[0] < isect1[0])
  {
    if (smallest1 == 0)
    {
      THEA_TRI3_SET(isectpt1, isectpointA1);
    }
    else
    {
      THEA_TRI3_SET(isectpt1, isectpointA2);
    }

    if (isect2[1] < isect1[1])
    {
      if (smallest2 == 0)
      {
        THEA_TRI3_SET(isectpt2, isectpointB2);
      }
      else
      {
        THEA_TRI3_SET(isectpt2, isectpointB1);
      }
    }
    else
    {
      if (smallest1 == 0)
      {
        THEA_TRI3_SET(isectpt2, isectpointA2);
      }
      else
      {
        THEA_TRI3_SET(isectpt2, isectpointA1);
      }
    }
  }
  else
  {
    if (smallest2 == 0)
    {
      THEA_TRI3_SET(isectpt1, isectpointB1);
    }
    else
    {
      THEA_TRI3_SET(isectpt1, isectpointB2);
    }

    if (isect2[1] > isect1[1])
    {
      if (smallest1 == 0)
      {
        THEA_TRI3_SET(isectpt2, isectpointA2);
      }
      else
      {
        THEA_TRI3_SET(isectpt2, isectpointA1);
      }
    }
    else
    {
      if (smallest2 == 0)
      {
        THEA_TRI3_SET(isectpt2, isectpointB2);
      }
      else
      {
        THEA_TRI3_SET(isectpt2, isectpointB1);
      }
    }
  }

  return 1;
}

#undef THEA_TRI3_USE_EPSILON_TEST
#undef THEA_TRI3_EPSILON
#undef THEA_TRI3_CROSS
#undef THEA_TRI3_DOT
#undef THEA_TRI3_SUB
#undef THEA_TRI3_ADD
#undef THEA_TRI3_MULT
#undef THEA_TRI3_SET
#undef THEA_TRI3_SORT
#undef THEA_TRI3_ISECT
#undef THEA_TRI3_COMPUTE_INTERVALS
#undef THEA_TRI3_EDGE_EDGE_TEST
#undef THEA_TRI3_EDGE_AGAINST_TRI_EDGES
#undef THEA_TRI3_POINT_IN_TRI
#undef THEA_TRI3_NEWCOMPUTE_INTERVALS
#undef THEA_TRI3_SORT2
#undef THEA_TRI3_COMPUTE_INTERVALS_ISECTLINE

template <typename T>
Vector<3, T> closestPointOnLineSegment(
  Vector<3, T> const  &  v0,
  Vector<3, T> const  &  v1,
  Vector<3, T> const  &  edgeDirection,
  T const                edgeLength,
  Vector<3, T> const  &  point)
{
  // Vector towards the point
  Vector<3, T> c = point - v0;
  // Projected onto the edge itself
  T t = edgeDirection.dot(c);

  if (t <= 0)
  {
    // Before the start
    return v0;
  }
  else if (t >= edgeLength)
  {
    // After the end
    return v1;
  }
  else
  {
    // At distance t along the edge
    return v0 + edgeDirection * t;
  }
}

template <typename T>
Vector<3, T>
closestPointOnTrianglePerimeter(
  Vector<3, T> const     v[3],
  Vector<3, T> const     edgeDirection[3],
  T const                edgeLength[3],
  Vector<3, T> const  &  point,
  int                 &  edgeIndex)
{
  // Closest point on segment from v[i] to v[i + 1]
  Vector<3, T> r[3];
  // Distance squared from r[i] to point
  T d[3];
  // Index of the next point
  static const int next[] = {1, 2, 0};

  for (int i = 0; i < 3; ++i)
  {
    r[i] = closestPointOnLineSegment(v[i], v[next[i]], edgeDirection[i], edgeLength[i], point);
    d[i] = (r[i] - point).squaredNorm();
  }

  if (d[0] < d[1])
  {
    if (d[0] < d[2])
    {
      // Between v0 and v1
      edgeIndex = 0;
    }
    else
    {
      // Between v2 and v0
      edgeIndex = 2;
    }
  }
  else
  {
    if (d[1] < d[2])
    {
      // Between v1 and v2
      edgeIndex = 1;
    }
    else
    {
      // Between v2 and v0
      edgeIndex = 2;
    }
  }

  return r[edgeIndex];
}

template <typename T>
Vector<3, T>
closestPointOnTrianglePerimeter(
  Vector<3, T> const & v0,
  Vector<3, T> const & v1,
  Vector<3, T> const & v2,
  Vector<3, T> const & point)
{
  Vector<3, T> v[3] = {v0, v1, v2};
  Vector<3, T> edgeDirection[3] = {(v1 - v0), (v2 - v1), (v0 - v2)};
  T edgeLength[3];

  for (int i = 0; i < 3; ++i)
  {
    edgeLength[i] = edgeDirection[i].norm();
    edgeDirection[i] /= edgeLength[i];
  }

  int edgeIndex;
  return closestPointOnTrianglePerimeter(v, edgeDirection, edgeLength, point, edgeIndex);
}

template <typename T>
bool
isPointInsideTriangle(
  Vector<3, T> const  &  v0,
  Vector<3, T> const  &  v1,
  Vector<3, T> const  &  v2,
  int                    primary_axis,
  Vector<3, T> const  &  p)
{
  // Check that the point is within the triangle using a Barycentric coordinate test on a two dimensional plane.
  int i, j;
  switch (primary_axis)
  {
    case 1          : i = 2; j = 0; break;
    case 2          : i = 0; j = 1; break;
    default /* 0 */ : i = 1; j = 2;
  }

  // See if all barycentric coordinates are non-negative

  // 2D area via cross product
#define THEA_TRI3_AREA2(d, e, f)  (((e)[i] - (d)[i]) * ((f)[j] - (d)[j]) - ((f)[i] - (d)[i]) * ((e)[j] - (d)[j]))

  // Area of the polygon
  T area = THEA_TRI3_AREA2(v0, v1, v2);
  if (area == 0)
  {
    // This triangle has zero area, so the point must not be in it unless the triangle point is the test point.
    return (v0 == p);
  }

  T inv_area = 1.0f / area;
  T b[3];

  // (Avoid normalization until absolutely necessary)
  b[0] = THEA_TRI3_AREA2(p, v1, v2) * inv_area;
  if ((b[0] < 0) || (b[0] > 1))
    return false;

  b[1] = THEA_TRI3_AREA2(v0, p, v2) * inv_area;
  if ((b[1] < 0) || (b[1] > 1))
    return false;

  b[2] = 1.0f - b[0] - b[1];

#undef THEA_TRI3_AREA2

  return ((b[2] >= 0) && (b[2] <= 1));
}

template <typename T>
T
rayTriangleIntersectionTime(RayN<3, T> const & ray, Vector<3, T> const & v0,
                            Vector<3, T> const & edge01, Vector<3, T> const & edge02)
{
  // The code is taken from Dave Eberly's Wild Magic library, v5.3, released under the Boost license:
  // http://www.boost.org/THEA_TRI3_LICENSE_1_0.txt .
  static T const EPSILON = Math::eps<T>();
  Vector<3, T> diff = ray.getOrigin() - v0;
  Vector<3, T> normal = edge01.cross(edge02);
  // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = ray direction, E1 = edge01, E2 = edge02, N = Cross(E1,E2)) by
  //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
  //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
  //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
  T DdN = ray.getDirection().dot(normal);
  int sign;

  if (DdN > EPSILON)
    sign = 1;
  else if (DdN < -EPSILON)
  {
    sign = -1;
    DdN = -DdN;
  }
  else
  {
    // Ray and triangle are parallel, call it a "no intersection" even if the ray does intersect
    return -1;
  }

  T DdQxE2 = sign * ray.getDirection().dot(diff.cross(edge02));

  if (DdQxE2 >= 0)
  {
    T DdE1xQ = sign * ray.getDirection().dot(edge01.cross(diff));

    if (DdE1xQ >= 0)
    {
      if (DdQxE2 + DdE1xQ <= DdN)
      {
        // Line intersects triangle, check if ray does
        T QdN = -sign * diff.dot(normal);

        if (QdN >= 0)
        {
          // Ray intersects triangle.
          return QdN / DdN;
        }

        // else: t < 0, no intersection
      }

      // else: b1 + b2 > 1, no intersection
    }

    // else: b2 < 0, no intersection
  }

  // else: b1 < 0, no intersection
  return -1;
}

} // namespace Triangle3Internal
} // namespace Thea

#endif
