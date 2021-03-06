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
//============== License for line-triangle distance calculation ==============
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
//
//============================================================================

#ifndef __Thea_Triangle3_hpp__
#define __Thea_Triangle3_hpp__

#include "Common.hpp"
#include "AxisAlignedBox3.hpp"
#include "Ball3.hpp"
#include "Box3.hpp"
#include "LineSegment3.hpp"
#include "Math.hpp"
#include "Plane3.hpp"
#include "RayIntersectable3.hpp"
#include <limits>

namespace Thea {

/**
 * Stores the positions of a triangle's three vertices locally, and provides access to them.
 *
 * @see Triangle3
 */
class THEA_API TriangleLocalVertexTriple3
{
  public:
    /** Default constructor. */
    TriangleLocalVertexTriple3() {}

    /** Initializing constructor. */
    TriangleLocalVertexTriple3(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
    }

    /** Get the i'th vertex. */
    Vector3 const & getVertex(int i) const { return vertices[i]; }

  private:
    Vector3 vertices[3];  ///< Vertex positions.

}; // TriangleLocalVertexTriple3

// Internal functions
namespace Triangle3Internal {

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * int tri_tri_intersect(Real const V0[3],Real const V1[3],Real const V2[3],
 *                       Real const U0[3],Real const U1[3],Real const U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDivTriTriIsect(Real const V0[3],Real const V1[3],Real const V2[3],
 *                      Real const U0[3],Real const U1[3],Real const U2[3]);
 *
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(Real const V0[3],Real const V1[3],Real const V2[3],
 *                                      Real const U0[3],Real const U1[3],Real const U2[3],int *coplanar,
 *                                      Real isectpt1[3],Real isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */

THEA_API int tri_tri_intersect(Real const V0[3],Real const V1[3],Real const V2[3],
                               Real const U0[3],Real const U1[3],Real const U2[3]);

THEA_API int NoDivTriTriIsect(Real const V0[3],Real const V1[3],Real const V2[3],
                              Real const U0[3],Real const U1[3],Real const U2[3]);

THEA_API int tri_tri_intersect_with_isectline(Real const V0[3],Real const V1[3],Real const V2[3],
                                              Real const U0[3],Real const U1[3],Real const U2[3],int *coplanar,
                                              Real isectpt1[3],Real isectpt2[3]);

// Check if a point is inside a triangle.
THEA_API bool isPointInsideTriangle(Vector3 const &  v0,
                                    Vector3 const &  v1,
                                    Vector3 const &  v2,
                                    int              primary_axis,
                                    Vector3 const &  p);

// Closest point on the perimeter of a triangle.
THEA_API Vector3 closestPointOnTrianglePerimeter(Vector3 const & v0,
                                                 Vector3 const & v1,
                                                 Vector3 const & v2,
                                                 Vector3 const & point);

// Intersection time of a ray with a triangle. Returns a negative value if the ray does not intersect the triangle.
THEA_API Real rayTriangleIntersectionTime(Ray3 const & ray, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02);

} // namespace Triangle3Internal

/**
 * Base class for a triangle in 3-space, with precomputed properties for fast access. To account for the fact that the triangle
 * vertices may be stored in different ways (e.g. internally within the object, or as indices into an external vertex pool), the
 * class is parametrized on the way the vertices are stored. <code>VertexTripleT</code> must be default-constructible and
 * provide an efficient member function with the signature
 *
 * \code
 * Vector3 [const &] getVertex(int i) const
 * \endcode
 *
 * @note This class cannot be used directly (it has protected constructors). Use Triangle3 or LocalTriangle3 instead.
 */
template <typename VertexTripleT = TriangleLocalVertexTriple3>
class /* THEA_DLL_LOCAL */ Triangle3Base : public RayIntersectable3
{
  public:
    typedef VertexTripleT VertexTriple;  ///< Stores and provides access to the triangle's vertices.

  protected:
    /** Default constructor. Does not initialize anything. */
    Triangle3Base() {}

    /** Construct from a set of vertices. */
    explicit Triangle3Base(VertexTriple const & vertices_) : vertices(vertices_) { update(); }

    /** Copy constructor. */
    Triangle3Base(Triangle3Base const & src)
    : vertices(src.vertices), plane(src.plane), primary_axis(src.primary_axis), centroid(src.centroid), edge01(src.edge01),
      edge02(src.edge02), area(src.area)
    {}

  public:
    /** Initialize the triangle from its vertices. */
    void set(VertexTriple const & vertices_) { vertices = vertices_; update(); }

    /**
     * Update the properties of the triangle, assuming the positions of its three corners have changed. Useful for external
     * callers if the mutable positions are not locally stored within the triangle.
     */
    void update() const
    {
      Vector3 v0 = getVertex(0), v1 = getVertex(1), v2 = getVertex(2);

      plane = Plane3::fromThreePoints(v0, v1, v2);
      primary_axis = (int)Math::maxAbsAxis(plane.getNormal());

      centroid = (v0 + v1 + v2) / 3;
      edge01   = v1 - v0;
      edge02   = v2 - v0;

      area = 0.5f * std::fabs(plane.getNormal().dot(edge01.cross(edge02)));  // exploit precomputed normal to avoid sqrt
    }

    /** Get a vertex of the triangle. */
    Vector3 getVertex(int i) const { return vertices.getVertex(i); }

    /** Get the vertices of the triangle. */
    VertexTriple const & getVertices() const { return vertices; }

    /** Get the plane of the triangle. */
    Plane3 const & getPlane() const { return plane; }

    /** Get the primary axis of the triangle (closest to normal). */
    intx getPrimaryAxis() const { return primary_axis; }

    /** Get the normal of the triangle (right-hand rule, going round vertices in order 0, 1, 2). */
    Vector3 const & getNormal() const { return plane.getNormal(); }

    /** Get the centroid of the triangle. */
    Vector3 const & getCentroid() const { return centroid; }

    /** Get the edge vector corresponding to getVertex(1) - getVertex(0). */
    Vector3 const & getEdge01() const { return edge01; }

    /** Get the edge vector corresponding to getVertex(2) - getVertex(0). */
    Vector3 const & getEdge02() const { return edge02; }

    /** Get the area of the triangle. */
    Real getArea() const { return area; }

    /** Get a uniformly distributed random sample from the triangle. */
    Vector3 randomPoint() const
    {
      // From G3D::Triangle

      // Choose a random point in the parallelogram
      float s = Random::common().uniform01();
      float t = Random::common().uniform01();

      if (s + t > 1.0f)
      {
        // Outside the triangle; reflect about the diagonal of the parallelogram
        t = 1.0f - t;
        s = 1.0f - s;
      }

      return s * edge01 + t * edge02 + getVertex(0);
    }

    /** Get the barycentric coordinates of a point assumed to be in the plane of the triangle. */
    Vector3 barycentricCoordinates(Vector3 const & p) const
    {
      if (area <= 0)
        return Vector3::Zero();

      Vector3 const & n = getNormal();
      Vector3 p0 = getVertex(0) - p;
      Vector3 p1 = getVertex(1) - p;
      Vector3 p2 = getVertex(2) - p;

      Real b0 = n.dot(p1.cross(p2)) / (2 * area);
      Real b1 = n.dot(p2.cross(p0)) / (2 * area);

      return Vector3(b0, b1, 1 - b0 - b1);
    }

    /** Get a bounding box for the triangle. */
    AxisAlignedBox3 getBounds() const
    {
      Vector3 v0 = getVertex(0), v1 = getVertex(1), v2 = getVertex(2);
      return AxisAlignedBox3(v0.cwiseMin(v1.cwiseMin(v2)), v0.cwiseMax(v1.cwiseMax(v2)));
    }

    /** Check if the triangle intersects (that is, contains) a point. */
    bool intersects(Vector3 const & p) const { return contains(p); }

    /** Check if the triangle intersects another triangle. */
    template <typename OtherVertexTripleT> bool intersects(Triangle3Base<OtherVertexTripleT> const & other) const
    {
      Vector3 p0 = getVertex(0);        Vector3 p1 = getVertex(1);        Vector3 p2 = getVertex(2);
      Vector3 q0 = other.getVertex(0);  Vector3 q1 = other.getVertex(1);  Vector3 q2 = other.getVertex(2);

      Real v0[3] = { p0.x(), p0.y(), p0.z() };
      Real v1[3] = { p1.x(), p1.y(), p1.z() };
      Real v2[3] = { p2.x(), p2.y(), p2.z() };

      Real u0[3] = { q0.x(), q0.y(), q0.z() };
      Real u1[3] = { q1.x(), q1.y(), q1.z() };
      Real u2[3] = { q2.x(), q2.y(), q2.z() };

      return Triangle3Internal::NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
    }

    /**
     * Check if the triangle intersects another triangle. If they do intersect, check if they are coplanar. If they are not
     * coplanar, compute the line of intersection. The intersection test is somewhat slower than intersects().
     */
    template <typename OtherVertexTripleT>
    bool intersects(Triangle3Base<OtherVertexTripleT> const & other, bool & coplanar, LineSegment3 & seg) const
    {
      Vector3 p0 = getVertex(0);        Vector3 p1 = getVertex(1);        Vector3 p2 = getVertex(2);
      Vector3 q0 = other.getVertex(0);  Vector3 q1 = other.getVertex(1);  Vector3 q2 = other.getVertex(2);

      int i_coplanar;
      Vector3 isectpt1, isectpt2;
      int isec = Triangle3Internal::tri_tri_intersect_with_isectline(&p0[0], &p1[0], &p2[0], &q0[0], &q1[0], &q2[0],
                                                                     &i_coplanar, &isectpt1[0], &isectpt2[0]);
      if (isec)
      {
        coplanar = (bool)i_coplanar;
        if (!coplanar)
          seg = LineSegment3(isectpt1, isectpt2);

        return true;
      }

      return false;
    }

    /** Check if the triangle intersects a ball. */
    bool intersects(Ball3 const & ball) const { throw Error("Triangle3: Intersection with ball not implemented"); }

    /** Check if the triangle intersects an axis-aligned box. */
    bool intersects(AxisAlignedBox3 const & aab) const { throw Error("Triangle3: Intersection with AAB not implemented"); }

    /** Check if the triangle intersects an oriented box. */
    bool intersects(Box3 const & box) const { throw Error("Triangle3: Intersection with oriented box not implemented"); }

    /** Check if the triangle contains a point. */
    bool contains(Vector3 const & p) const
    {
      return Triangle3Internal::isPointInsideTriangle(getVertex(0), getVertex(1), getVertex(2), primary_axis, p);
    }

    /** Get the distance of the triangle from a point. */
    Real distance(Vector3 const & p) const { return std::sqrt(squaredDistance(p)); }

    /** Get the squared distance of the triangle from a point. */
    Real squaredDistance(Vector3 const & p) const
    {
      return (closestPoint(p) - p).squaredNorm();
    }

    /** Get the point on this triangle closest to a given point. */
    Vector3 closestPoint(Vector3 const & p) const
    {
      // Project the point onto the plane of the triangle
      Vector3 proj = plane.closestPoint(p);

      if (contains(proj))
        return proj;
      else  // the closest point is on the perimeter instead
        return Triangle3Internal::closestPointOnTrianglePerimeter(getVertex(0), getVertex(1), getVertex(2), p);
    }

    /** Get the distance of the triangle from another triangle. */
    template <typename OtherVertexTripleT>
    Real distance(Triangle3Base<OtherVertexTripleT> const & other) const { return std::sqrt(squaredDistance(other)); }

    /**
     * Get the squared distance between this triangle and another triangle, and optionally return the closest pair of points.
     */
    template <typename OtherVertexTripleT>
    Real squaredDistance(Triangle3Base<OtherVertexTripleT> const & other, Vector3 * this_pt = nullptr,
                         Vector3 * other_pt = nullptr) const
    {
      // From Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.

      Vector3 c1(0, 0, 0), c2(0, 0, 0);  // squash an uninitialized variable warning
      Real min_sqdist = std::numeric_limits<Real>::max();

      // First test for intersection
      bool coplanar;
      LineSegment3 seg;
      if (intersects(other, coplanar, seg))
      {
        if (coplanar)
          c1 = c2 = (Real)0.5 * (getCentroid() + other.getCentroid());  // FIXME: This needn't be in the intersection
        else
          c1 = c2 = (Real)0.5 * (seg.getPoint(0) + seg.getPoint(1));

        min_sqdist = 0;
      }
      else
      {
        // Edge-edge distances
        Vector3 p, q;
        Real d2, s, t;
        for (int i = 0; i < 3; ++i)
        {
          int i2 = (i + 1) % 3;

          for (int j = 0; j < 3; ++j)
          {
            int j2 = (j + 1) % 3;
            d2 = Internal::closestPtSegmentSegment<3, Real>(getVertex(i), getVertex(i2), false,
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
          p = getVertex(i);
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
    Real distance(Ball3 const & ball) const
    {
      return std::max(distance(ball.getCenter()) - ball.getRadius(), static_cast<Real>(0));
    }

    /** Get the squared distance between this triangle and a ball, and optionally return the closest pair of points. */
    Real squaredDistance(Ball3 const & ball, Vector3 * this_pt = nullptr, Vector3 * ball_pt = nullptr) const
    {
      if (!this_pt && !ball_pt)
      {
        Real x = distance(ball);
        return x * x;
      }

      Vector3 c1 = closestPoint(ball.getCenter());
      Vector3 c2;

      Vector3 diff = c1 - ball.getCenter();
      Real d2 = diff.squaredNorm();
      Real r2 = ball.getRadius() * ball.getRadius();
      if (d2 < r2)  // point inside ball
      {
        c2 = c1;
        d2 = 0;
      }
      else
      {
        if (r2 < 1e-30)
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
    Real distance(Line3 const & line) const
    {
      return std::sqrt(squaredDistance(line));
    }

    /**
     * Get the squared distance between this triangle and an infinite line, and optionally return the closest pair of points.
     */
    Real squaredDistance(Line3 const & line, Vector3 * this_pt = nullptr, Vector3 * line_pt = nullptr) const
    {
      // Two tests where one would suffice, but for now it avoids having to code a new function or modify an existing one.
      // Shift the ray origins to ensure the rays overlap and there are no errors at line.getPoint().

      // Forward half
      {
        Ray3 ray(line.getPoint() - line.getDirection(), line.getDirection());
        Real t = rayIntersectionTime(ray);
        if (t >= 0)
        {
          Vector3 c = ray.getPoint(t);
          if (this_pt) *this_pt = c;
          if (line_pt) *line_pt = c;
          return 0;
        }
      }

      // Backward half
      {
        Ray3 ray(line.getPoint() + line.getDirection(), -line.getDirection());
        Real t = rayIntersectionTime(ray);
        if (t >= 0)
        {
          Vector3 c = ray.getPoint(t);
          if (this_pt) *this_pt = c;
          if (line_pt) *line_pt = c;
          return 0;
        }
      }

      return squaredDistanceToPerimeter(line, this_pt, line_pt);
    }

    /** Get the distance of this triangle from a line segment. */
    Real distance(LineSegment3 const & seg) const
    {
      return std::sqrt(squaredDistance(seg));
    }

    /**
     * Get the squared distance between this triangle and a line segment, and optionally return the closest pair of points.
     */
    Real squaredDistance(LineSegment3 const & seg, Vector3 * this_pt = nullptr, Vector3 * seg_pt = nullptr) const
    {
      // From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      Ray3 ray(seg.getEndpoint(0), seg.getDirection());
      Real t = rayIntersectionTime(ray);
      if (t >= 0 && t <= 1)
      {
        Vector3 c = ray.getPoint(t);
        if (this_pt) *this_pt = c;
        if (seg_pt)  *seg_pt  = c;
        return 0;
      }

      return squaredDistanceToPerimeter(seg, this_pt, seg_pt);
    }

    /** Get the distance of this triangle from a ray. */
    Real distance(Ray3 const & ray) const
    {
      return std::sqrt(squaredDistance(ray));
    }

    /**
     * Get the squared distance between this triangle and a ray, and optionally return the closest pair of points.
     */
    Real squaredDistance(Ray3 const & ray, Vector3 * this_pt = nullptr, Vector3 * ray_pt = nullptr) const
    {
      // From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      Real t = rayIntersectionTime(ray);
      if (t >= 0)
      {
        Vector3 c = ray.getPoint(t);
        if (this_pt) *this_pt = c;
        if (ray_pt)  *ray_pt  = c;
        return 0;
      }

      return squaredDistanceToPerimeter(ray, this_pt, ray_pt);
    }

    Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const
    {
      Real t = Triangle3Internal::rayTriangleIntersectionTime(ray, getVertex(0), getEdge01(), getEdge02());
      return (max_time >= 0 && t > max_time) ? -1 : t;
    }

    RayIntersection3 rayIntersection(Ray3 const & ray, Real max_time = -1) const
    {
      Real t = Triangle3Internal::rayTriangleIntersectionTime(ray, getVertex(0), getEdge01(), getEdge02());
      if (t >= 0 && (max_time < 0 || t <= max_time))
      {
        Vector3 n = getNormal();
        return RayIntersection3(t, &n);
      }

      return RayIntersection3(-1);
    }

  protected:
    /**
     * Get the squared distance between the perimiter of this triangle and a line-like object (Line3, Ray3 or LineSegment3), and
     * optionally return the closest pair of points.
     */
    template <typename LineLikeT>
    Real squaredDistanceToPerimeter(LineLikeT const & line, Vector3 * this_pt = nullptr, Vector3 * line_pt = nullptr) const
    {
      // From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistLine3Triangle3.h

      Real d2 = -1;
      Vector3 c1, c2;
      for (int i = 0; i < 3; ++i)
      {
        LineSegment3 edge(getVertex(i), getVertex((i + 1) % 3));
        Real edge_d2 = edge.squaredDistance(line, (this_pt ? &c1 : nullptr), (line_pt ? &c2 : nullptr));
        if (d2 < 0 || edge_d2 < d2)
        {
          d2 = edge_d2;
          if (this_pt) *this_pt = c1;
          if (line_pt) *line_pt = c2;
        }
      }

      return d2;
    }

    VertexTriple  vertices;  ///< The vertices of the triangle.

    // Cached properties computed from the vertices
    mutable Plane3   plane;         ///< Plane of the triangle.
    mutable int      primary_axis;  ///< Primary axis (closest to normal).
    mutable Vector3  centroid;      ///< Centroid of the triangle (mean of three vertices).
    mutable Vector3  edge01;        ///< vertices[1] - vertices[0]
    mutable Vector3  edge02;        ///< vertices[2] - vertices[0]
    mutable Real     area;          ///< Triangle area.

}; // class Triangle3Base

// Forward declaration
template <typename VertexTripleT = TriangleLocalVertexTriple3> class Triangle3;

/**
 * A triangle with three vertex positions stored locally, in the class itself. This class adds a more direct constructor and
 * set() method for convenience to the default Triangle3 template.
 *
 * @see Triangle3
 */
template <>
class THEA_API Triangle3<TriangleLocalVertexTriple3> : public Triangle3Base<TriangleLocalVertexTriple3>
{
  private:
    typedef Triangle3Base<TriangleLocalVertexTriple3> BaseT;

  public:
    THEA_DECL_SMART_POINTERS(Triangle3)

    /** Default constructor. Does not initialize anything. */
    Triangle3() {}

    /** Construct from a set of vertices. */
    explicit Triangle3(TriangleLocalVertexTriple3 const & vertices_) : BaseT(vertices_) {}

    /** Construct from a set of vertices. */
    Triangle3(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2) : BaseT(TriangleLocalVertexTriple3(v0, v1, v2)) {}

    /** Copy constructor. */
    Triangle3(Triangle3 const & src) : BaseT(src) {}

    /** Transform the triangle by a 4x4 matrix and return the result. */
    Triangle3 transform(Matrix4 const & m) const
    {
      return Triangle3(Math::hmul(m, BaseT::getVertex(0)),
                       Math::hmul(m, BaseT::getVertex(1)),
                       Math::hmul(m, BaseT::getVertex(2)));
    }

    /** Initialize the triangle from its vertices. */
    void set(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2)
    {
      BaseT::set(TriangleLocalVertexTriple3(v0, v1, v2));
      update();
    }

    /** Get a copy of this triangle (for consistency with the default template). */
    Triangle3 localClone() const { return *this; }

}; // class Triangle3<TriangleLocalVertexTriple3>

/** A triangle with three vertex positions stored locally, in the class itself. */
typedef Triangle3<TriangleLocalVertexTriple3> LocalTriangle3;

/**
 * A triangle in 3-space, with precomputed properties for fast access. To account for the fact that the triangle vertices may
 * be stored in different ways (e.g. internally within the object, or as indices into an external vertex pool), the class is
 * parametrized on the way the vertices are stored. <code>VertexTripleT</code> must be default-constructible and provide an
 * efficient member function with the signature
 *
 * \code
 * Vector3 const & getVertex(int i) const
 * \endcode
 */
template <typename VertexTripleT>
class /* THEA_API */ Triangle3 : public Triangle3Base<VertexTripleT>
{
  private:
    typedef Triangle3Base<VertexTripleT> BaseT;

  public:
    THEA_DECL_SMART_POINTERS(Triangle3)

    /** Default constructor. Does not initialize anything. */
    Triangle3() {}

    /** Construct from a set of vertices. */
    explicit Triangle3(VertexTripleT const & vertices_) : BaseT(vertices_) {}

    /** Copy constructor. */
    Triangle3(Triangle3 const & src) : BaseT(src) {}

    /** Transform the triangle by a 4x4 matrix and return the result. */
    LocalTriangle3 transform(Matrix4 const & m) const
    {
      return LocalTriangle3(m * BaseT::getVertex(0), m * BaseT::getVertex(1), m * BaseT::getVertex(2));
    }

    /** Get a new triangle that simply stores copies of the vertex positions of this triangle. */
    LocalTriangle3 localClone() const
    {
      return LocalTriangle3(BaseT::getVertex(0), BaseT::getVertex(1), BaseT::getVertex(2));
    }

}; // class Triangle3

} // namespace Thea

#endif
