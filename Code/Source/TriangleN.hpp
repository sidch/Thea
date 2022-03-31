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
//============================================================================

#ifndef __Thea_Triangle_hpp__
#define __Thea_Triangle_hpp__

#include "Common.hpp"
#include "MatVec.hpp"
#include <type_traits>
#include <limits>

namespace Thea {

/**
 * Stores the positions of a triangle's three vertices locally, and provides access to them.
 *
 * @see TriangleN
 */
template <int N, typename T = Real>
class /* THEA_API */ TriangleLocalVertexTripleN
{
  public:
    typedef Vector<N, T> VectorT;  ///< A vector in N-D space.

    /** Default constructor. */
    TriangleLocalVertexTripleN() {}

    /** Initializing constructor. */
    TriangleLocalVertexTripleN(VectorT const & v0, VectorT const & v1, VectorT const & v2) : vertices{ v0, v1, v2 } {}

    /** Get the i'th vertex. */
    VectorT const & getVertex(int i) const { return vertices[i]; }

    /** Set the i'th vertex. */
    void setVertex(int i, VectorT const & v)
    {
      debugAssertM(i >= 0 && i < 3, "TriangleLocalVertexTripleN: Vertex index out of bounds");

      vertices[i] = v;
    }

    /** Set all three vertices at once. */
    void set(VectorT const & v0, VectorT const & v1, VectorT const & v2)
    {
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
    }

  private:
    VectorT vertices[3];  ///< Vertex positions.

}; // TriangleLocalVertexTripleN

// Forward declarations
template < int N, typename T = Real, typename VertexTripleT = TriangleLocalVertexTripleN<N, T> > class /* THEA_API */ TriangleN;

/**
 * Base class for a triangle (convex hull of 3 points) in N-space. To account for the fact that the triangle vertices may be
 * stored in different ways (e.g. internally within the object, or as indices into an external vertex pool), the class is
 * parametrized on the way the vertices are stored. <code>VertexTripleT</code> must be default-constructible, copy-constructible
 * and assignable, and provide an efficient member function with the signature:
 * \code
 * Vector3 [const &] getVertex(int i) const
 * \endcode
 *
 * @warning This class cannot be used directly (it has protected constructors). Use TriangleN instead to access functions of this
 *   class.
 */
template <int N, typename VertexTripleT, typename T>
class /* THEA_DLL_LOCAL */ TriangleNBase
{
  public:
    typedef Vector<N, T>   VectorT;       ///< A vector in N-D space.
    typedef VertexTripleT  VertexTriple;  ///< Stores and provides access to the triangle's vertices.

    /** Get the wrapped vertex triple. */
    VertexTriple const & getVertices() const { return vertices; }

    /** Get the i'th vertex. */
    VectorT const & getVertex(int i) const { return vertices.getVertex(i); }

    /** Transform the triangle and return the result. */
    template < typename TransformT, typename std::enable_if< !Math::IsHomMatrix<TransformT, N>::value >::type * = nullptr >
    TriangleN< N, TriangleLocalVertexTripleN<N, T>, T > transform(TransformT const & tr) const
    {
      return TriangleN< N, TriangleLocalVertexTripleN<N, T>, T >(tr * getVertex(0), tr * getVertex(1), tr * getVertex(2));
    }

    /** Transform the triangle by a homogeneous matrix and return the result. */
    template < typename MatrixT, typename std::enable_if< Math::IsHomMatrix<MatrixT, N>::value >::type * = nullptr >
    TriangleN< N, TriangleLocalVertexTripleN<N, T>, T > transform(MatrixT const & m) const
    {
      return TriangleN< N, TriangleLocalVertexTripleN<N, T>, T >(Math::hmul(m, getVertex(0)),
                                                                 Math::hmul(m, getVertex(1)),
                                                                 Math::hmul(m, getVertex(2)));
    }

    /** Get a new triangle that simply stores copies of the vertex positions of this triangle. */
    TriangleN< N, TriangleLocalVertexTripleN<N, T>, T > localClone() const
    {
      return TriangleN< N, TriangleLocalVertexTripleN<N, T>, T >(getVertex(0), getVertex(1), getVertex(2));
    }

  protected:
    /* Default constructor. */
    TriangleNBase() {}

    /** Construct from a set of vertices. */
    explicit TriangleNBase(VertexTriple const & vertices_) : vertices(vertices_) {}

    /** Initializing constructor, callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or derived from it. */
    TriangleNBase(VectorT const & v0, VectorT const & v1, VectorT const & v2)
    {
      static_assert(std::is_base_of< TriangleLocalVertexTripleN<N, T>, VertexTriple >::value,
                    "TriangleNBase: Cannot construct non-local vertex triple from three points");

      vertices.set(v0, v1, v2);
    }

    /** Copy constructor. */
    TriangleNBase(TriangleNBase const & src) : vertices(src.vertices) {}

    /** Initialize the triangle from its vertex triple. */
    void set(VertexTriple const & vertices_) { vertices = vertices_; }

    /**
     * Set the i'th vertex. This function is callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or derived from
     * it.
     */
    void setVertex(int i, VectorT const & v)
    {
      static_assert(std::is_base_of< TriangleLocalVertexTripleN<N, T>, VertexTriple >::value,
                    "TriangleNBase: Cannot set vertex of non-local vertex triple");

      vertices.setVertex(i, v);
    }

    /**
     * Set all three vertices at once. This function is callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or
     * derived from it.
     */
    void set(VectorT const & v0, VectorT const & v1, VectorT const & v2)
    {
      static_assert(std::is_base_of< TriangleLocalVertexTripleN<N, T>, VertexTriple >::value,
                    "TriangleNBase: Cannot set vertices of non-local vertex triple");

      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
    }

  private:
    VertexTripleT vertices;  ///< The three vertices of the triangle, accessed through a common interface.

}; // class TriangleNBase

/**
 * A triangle (convex hull of 3 points) in N-space. To account for the fact that the triangle vertices may be stored in
 * different ways (e.g. internally within the object, or as indices into an external vertex pool), the class is parametrized on
 * the way the vertices are stored. <code>VertexTripleT</code> must be default-constructible, copy-constructible and assignable,
 * and provide an efficient member function with the signature:
 * \code
 * VectorT const & getVertex(int i) const
 * \endcode
 *
 * The default implementation contains very few functions. The specializations TriangleN<2, ...> and TriangleN<3, ...> contain
 * most of the useful functionality.
 */
template <int N, typename VertexTripleT, typename T>
class /* THEA_API */ TriangleN : public TriangleNBase<N, VertexTripleT, T>
{
  public:
    typedef TriangleNBase<N, VertexTripleT, T> BaseT;  ///< Base class with function implementations.
    typedef typename BaseT::VectorT VectorT;           ///< A vector in N-D space.

    /** Default constructor. Does not initialize anything. */
    TriangleN() {}

    /** Construct from a set of vertices. */
    explicit TriangleN(VertexTripleT const & vertices_) : BaseT(vertices_) {}

    /** Initializing constructor, callable only if VertexTriple is TriangleLocalVertexTripleN<N, T> or derived from it. */
    TriangleN(VectorT const & v0, VectorT const & v1, VectorT const & v2) : BaseT(v0, v1, v2) {}

}; // class TriangleN

/** A convenience typedef for a triangle whose vertex positions are stored locally in the TriangleN object itself. */
template <int N, typename T = Real> using LocalTriangleN = TriangleN< N, TriangleLocalVertexTripleN<N, T>, T >;

} // namespace Thea

#include "Triangle3.hpp"

#endif
