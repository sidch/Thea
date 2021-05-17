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
//=================================== Original header ========================
//
// COTD Entry submitted by John W. Ratcliff [jratcliff@verant.com]
//
// ** THIS IS A CODE SNIPPET WHICH WILL EFFICIEINTLY TRIANGULATE ANY
// ** POLYGON/CONTOUR (without holes) AS A STATIC CLASS.  THIS SNIPPET
// ** IS COMPRISED OF 3 FILES, TRIANGULATE.H, THE HEADER FILE FOR THE
// ** TRIANGULATE BASE CLASS, TRIANGULATE.CPP, THE IMPLEMENTATION OF
// ** THE TRIANGULATE BASE CLASS, AND TEST.CPP, A SMALL TEST PROGRAM
// ** DEMONSTRATING THE USAGE OF THE TRIANGULATOR.  THE TRIANGULATE
// ** BASE CLASS ALSO PROVIDES TWO USEFUL HELPER METHODS, ONE WHICH
// ** COMPUTES THE AREA OF A POLYGON, AND ANOTHER WHICH DOES AN EFFICENT
// ** POINT IN A TRIANGLE TEST.
// ** SUBMITTED BY JOHN W. RATCLIFF (jratcliff@verant.com) July 22, 2000
//
// Static class to triangulate any contour/polygon efficiently
// You should replace Vector2d with whatever your own Vector
// class might be.  Does not support polygons with holes.
// Uses STL vectors to represent a dynamic array of vertices.
// This code snippet was submitted to FlipCode.com by
// John W. Ratcliff (jratcliff@verant.com) on July 22, 2000
// I did not write the original code/algorithm for this
// this triangulator, in fact, I can't even remember where I
// found it in the first place.  However, I did rework it into
// the following black-box static class so you can make easy
// use of it in your own code.  Simply replace Vector2d with
// whatever your own Vector implementation might be.
//
//============================================================================

#ifndef __Thea_Polygon3_hpp__
#define __Thea_Polygon3_hpp__

#include "Common.hpp"
#include "Algorithms/PointTraitsN.hpp"
#include "Array.hpp"
#include "AxisAlignedBox3.hpp"
#include "Math.hpp"
#include "MatVec.hpp"

namespace Thea {

namespace Polygon3Internal {

/** A vertex plus an index. */
struct THEA_API IndexedVertex
{
  /** Default constructor. */
  IndexedVertex() {}

  /** Initializing constructor. */
  IndexedVertex(Vector3 const & position_, intx index_) : position(position_), index(index_) {}

  Vector3 position;  ///< The position of the vertex.
  intx index;  ///< The index of the vertex.
};

} // namespace Polygon3Internal

namespace Algorithms {

// Specify that a Polygon3 vertex is a logical 3D point. */
template <>
class IsPointN<Polygon3Internal::IndexedVertex, 3>
{
  public:
    static bool const value = true;
};

// Map a Polygon3 vertex to its 3D position. */
template <>
class PointTraitsN<Polygon3Internal::IndexedVertex, 3>
{
  public:
    static Vector3 const & getPosition(Polygon3Internal::IndexedVertex const & t) { return t.position; }
};

} // namespace Algorithms

/** A polygon in 3-space. Original code due to John W. Ratcliff. */
class THEA_API Polygon3
{
  public:
    THEA_DECL_SMART_POINTERS(Polygon3)

    typedef Polygon3Internal::IndexedVertex IndexedVertex;  ///< A vertex plus an index.

    /** Construct an empty polygon. */
    Polygon3();

    /**
     * Add a vertex to the polygon. The vertex is inserted at the end of the current sequence of vertices, and by default is
     * assigned an index that is one more than the maximum index in the polygon so far (or zero if this is the first vertex).
     * For efficiency the polygon <b>is not checked for planarity</b>, the caller should ensure that all vertices are coplanar.
     */
    void addVertex(Vector3 const & p);

    /**
     * Add an indexed vertex to the polygon. The vertex is inserted at the end of the current sequence of vertices. For
     * efficiency the polygon <b>is not checked for planarity</b>, the caller should ensure that all vertices are coplanar.
     */
    void addVertex(Vector3 const & p, intx index);

    /** Get the number of vertices in the polygon. */
    intx numVertices() const;

    /**
     * Get the vertex at position \a poly_index in the sequence of vertices around the polygon boundary.
     *
     * @note \a poly_index is determined by the sequence of addVertex() calls, <b>NOT</b> by the index supplied in
     *   addVertex(Vector2 const &, intx)!
     */
    IndexedVertex const & getVertex(intx poly_index) const;

    /** Delete all vertices from the polygon. */
    void clear();

    /**
     * Triangulate the polygon and return the set of triangle indices (in successive groups of 3). All prior data in the
     * supplied array are cleared.
     *
     * @return The number of triangles created.
     */
    intx triangulate(Array<intx> & tri_indices, Real epsilon = -1) const;

    /** Compute the area of the polygon. */
    Real computeArea() const { return computeArea(vertices.begin(), vertices.end()); }

    /** Compute the area of an arbitrary polygon (does not require explicit creation of a Polygon3). */
    template <typename VertexInputIterator> static Real computeArea(VertexInputIterator vbegin, VertexInputIterator vend)
    {
      typedef typename std::iterator_traits<VertexInputIterator>::value_type VertexT;

      static Real const EPSILON = 1e-10f;

      Vector3 normal = computeNormal(vbegin, vend);
      if (normal.squaredNorm() < EPSILON)
        return 0;

      VertexInputIterator v0 = vbegin;
      VertexInputIterator v1 = incrementIterator(v0, vbegin, vend);

      if (v1 == v0)  // polygon is a point
        return 0;

      Vector3 p0 = Algorithms::PointTraitsN<VertexT, 3>::getPosition(*v0);
      Vector3 c = Vector3::Zero();
      for (VertexInputIterator vi = vbegin; vi != vend; ++vi)
      {
        Vector3 p1 = Algorithms::PointTraitsN<VertexT, 3>::getPosition(*v1);
        c += p0.cross(p1);

        p0 = p1;
        v1 = incrementIterator(v1, vbegin, vend);
      }

      return std::fabs(0.5f * c.dot(normal));
    }

    /** Compute the unit normal of the polygon. */
    Vector3 computeNormal() const { return computeNormal(vertices.begin(), vertices.end()); }

    /**
     * Compute the unit normal of an arbitrary polygon (does not require explicit creation of a Polygon3). Returns the zero
     * vector if the polygon is degenerate or too small. PointTraitsN<T, 3> must be defined for the value_type T of the
     * iterator.
     */
    template <typename VertexInputIterator> static Vector3 computeNormal(VertexInputIterator vbegin, VertexInputIterator vend)
    {
      typedef typename std::iterator_traits<VertexInputIterator>::value_type VertexT;

      static Real const EPSILON = 1e-10f;

      VertexInputIterator v0 = vbegin;
      VertexInputIterator v1 = incrementIterator(v0, vbegin, vend);
      VertexInputIterator v2 = incrementIterator(v1, vbegin, vend);

      if (v1 == v0 || v2 == v0)  // too few vertices
        return Vector3::Zero();

      Vector3 p0 = Algorithms::PointTraitsN<VertexT, 3>::getPosition(*v0);
      Vector3 p1 = Algorithms::PointTraitsN<VertexT, 3>::getPosition(*v1);

      for (VertexInputIterator vi = vbegin; vi != vend; ++vi)
      {
        Vector3 p2 = Algorithms::PointTraitsN<VertexT, 3>::getPosition(*v2);
        Vector3 normal = (p2 - p1).cross(p0 - p1);
        if (normal.squaredNorm() >= EPSILON)
          return normal.normalized();

        p0 = p1;
        p1 = p2;
        v2 = incrementIterator(v2, vbegin, vend);
      }

      return Vector3::Zero();  // degenerate or too small
    }

    /** Get the bounding box of the polygon. */
    AxisAlignedBox3 const & getBounds();

    /**
     * Utility function to split a (possibly non-convex) quadrilateral into a pair of triangles with the same winding direction.
     * The quadrilateral is assumed to be planar. This function guarantees that if two triangles are returned (the
     * non-degenerate case), then (a) the shared diagonal of the parent quad will be the first and third vertices of each
     * triangle (but in opposite orders to preserve orientation), and (b) the vertices of the first triangle will be consecutive
     * vertices of the quad. In other words, the vertices of the original quad, starting at some vertex (not necessarily \a p0)
     * and in the same winding order, can be recovered as \a i0, \a j0, \a k0, \a j1. Note that the generic triangulate() does
     * not have any such guarantee.
     *
     * @param p0 Position of the first vertex of the quadrilateral.
     * @param p1 Position of the second vertex of the quadrilateral.
     * @param p2 Position of the third vertex of the quadrilateral.
     * @param p3 Position of the fourth vertex of the quadrilateral.
     *
     * @param i0 Used to return the index of the first vertex of the first triangle.
     * @param j0 Used to return the index of the second vertex of the first triangle. This vertex is guaranteed to not be part
     *   of the second triangle.
     * @param k0 Used to return the index of the third vertex of the first triangle.
     *
     * @param i1 Used to return the index of the first vertex of the second triangle.
     * @param j1 Used to return the index of the second vertex of the second triangle. This vertex is guaranteed to not be part
     *   of the first triangle.
     * @param k1 Used to return the index of the third vertex of the second triangle.
     *
     * @param epsilon A tolerance threshold for checking degeneracy, negative for a default value.
     *
     * @return The number of triangles produced (can be < 2 if the quadrilateral is degenerate).
     *
     * @note Might break in very very degenerate cases (not fully tested).
     */
    template <typename T>
    static int triangulateQuad(Vector<3, T> const & p0, Vector<3, T> const & p1,
                               Vector<3, T> const & p2, Vector<3, T> const & p3,
                               intx & i0, intx & j0, intx & k0,
                               intx & i1, intx & j1, intx & k1,
                               T const & epsilon = -1)
    {
      typedef Vector<3, T> VectorT;

      // We have two diagonals to split along. Try one and if it produces a triangle outside the polygon, pick the other one

      // First diagonal is p0-p2
      VectorT n0 = (p1 - p0).cross(p2 - p0);
      VectorT n1 = (p2 - p0).cross(p3 - p0);

      if (n0.dot(n1) < 0)
      {
        // Flip to diagonal p1-p3
        n0 = (p2 - p1).cross(p3 - p1);
        n1 = (p3 - p1).cross(p0 - p1);

        i0 = 1; j0 = 2; k0 = 3;
        i1 = 3; j1 = 0; k1 = 1;
      }
      else
      {
        i0 = 0; j0 = 1; k0 = 2;
        i1 = 2; j1 = 3; k1 = 0;
      }

      // Check for degenerate triangles
      T e2 = (epsilon < 0 ? Math::eps<T>() * Math::eps<T>() : epsilon * epsilon);
      if (n0.squaredNorm() < e2)
      {
        if (n1.squaredNorm() < e2)
          return 0;
        else
        {
          i0 = i1; j0 = j1; k0 = k1;  // the second triangle is the only non-degenerate one
          return 1;
        }
      }
      else if (n1.squaredNorm() < e2)
        return 1;
      else
        return 2;
    }

  private:
    /** Signed area of projection onto primary coordinate plane. */
    Real projArea() const;

    /** Check if a triangle can be removed. */
    bool snip(size_t u, size_t v, size_t w, size_t n, Array<size_t> const & indices, Real epsilon) const;

    /** Advance an iterator round the polygon with vertices [vbegin, vend). */
    template <typename VertexInputIterator>
    static VertexInputIterator incrementIterator(VertexInputIterator vi, VertexInputIterator vbegin, VertexInputIterator vend)
    {
      if (++vi != vend)
        return vi;
      else
        return vbegin;
    }

    Array<IndexedVertex> vertices;
    intx max_index;
    AxisAlignedBox3 bounds;
    mutable Array<Vector2> proj_vertices;

    friend class Polygon2;

}; // class Polygon3

} // namespace Thea

#endif
