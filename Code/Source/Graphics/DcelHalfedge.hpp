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

#ifndef __Thea_Graphics_DcelHalfedge_hpp__
#define __Thea_Graphics_DcelHalfedge_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "GraphicsAttributes.hpp"
#include <iterator>

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DcelMesh;
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DcelFace;
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DcelVertex;

namespace DCELInternal {

// Iterate over edges in a collection (every other halfedge).
template <typename BaseIterT>
struct /** THEA_API */ BidirEdgeIterator : public BaseIterT
{
  // Default constructor.
  BidirEdgeIterator() {}

  // General copy constructor.
  template <typename T> BidirEdgeIterator(T const & src) : BaseIterT(src) {}

  // Assignment. Use with caution: \a src must have even index.
  BidirEdgeIterator & operator=(BaseIterT const & src) { BaseIterT::operator=(src); }

  // Pre-increment. */
  BidirEdgeIterator & operator++()
  {
    BaseIterT::operator++();
    BaseIterT::operator++();
    return *this;
  }

  // Pre-decrement.
  BidirEdgeIterator & operator--()
  {
    BaseIterT::operator--();
    BaseIterT::operator--();
    return *this;
  }

  // Post-increment.
  BidirEdgeIterator operator++(int)
  {
    BidirEdgeIterator ret = *this;
    BaseIterT::operator++();
    BaseIterT::operator++();
    return ret;
  }

  // Post-decrement.
  BidirEdgeIterator operator--(int)
  {
    BidirEdgeIterator ret = *this;
    BaseIterT::operator--();
    BaseIterT::operator--();
    return ret;
  }

}; // class BidirEdgeIterator

// Iterate over a circular loop of edges, derefencing to some property of the edge accessed by a pointer. EdgeT should be either
// DcelHalfedge or DcelHalfedge const. This tries to map a circular loop to an iterator with a well-defined begin and end, by
// keeping track of where the iteration began, and allowing only forward progress. There are clear caveats with this approach
// but it should fit a majority of use cases. Use with care.
//
// DerefFunctorT should be callable with the signature: ValueT * operator()(EdgeT *)
// IncrementFunctorT should be callable with the signature: void operator()(EdgeT **)
template <typename EdgeT, typename ValueT, typename DerefFunctorT, typename IncrementFunctorT>
class /** THEA_API */ FwdIterator
: public std::iterator<std::forward_iterator_tag, ValueT, std::ptrdiff_t, ValueT *, ValueT &>
{
  public:
    // Constructor.
    FwdIterator(EdgeT * e_ = nullptr, bool first_ = true) : initial(e_), e(e_), first(first_) {}

    // Construct from a compatible iterator, typically a non-const to const conversion. Compatibility is the responsibility of
    // the caller.
    template <typename E2, typename V2, typename D2, typename I2> explicit FwdIterator(FwdIterator<E2, V2, D2, I2> const & src)
    : initial(src.getInitial()), e(src.getCurrent()), first(src.isAtInitial()) {}

    // Assign from a compatible iterator, typically a non-const to const conversion. Compatibility is the responsibility of
    // the caller.
    template <typename OtherFwdIteratorT> FwdIterator & operator=(OtherFwdIteratorT const & src)
    {
      initial = src.initial;
      e = src.e;
      first = src.first;

      return *this;
    }

    // Equality comparison.
    bool operator==(FwdIterator const & other)
    {
      // True if the current edges match AND (both are null OR neither is at its initial position).
      return e == other.e && (!e || ((e != initial || !first) && (other.e != other.initial || !other.first)));
    }

    // Inequality comparison.
    bool operator!=(FwdIterator const & other) { return !operator==(other); }

    // Dereference.
    ValueT * operator*() const { return DerefFunctorT()(e); }

    // Arrow operator.
    ValueT const * const * operator->() const
    { alwaysAssertM(false, "FwdIterator: Can't call '->' on iterator-over-pointers"); }

    // Pre-increment. */
    FwdIterator & operator++()
    {
      if (e)
      {
        IncrementFunctorT()(&e);
        first = false;
      }

      return *this;
    }

    // Post-increment.
    FwdIterator operator++(int)
    {
      FwdIterator ret = *this;
      if (e)
      {
        IncrementFunctorT()(&e);
        first = false;
      }

      return ret;
    }

    // Get the initial edge from which the iterator was created.
    EdgeT * getInitial() const { return initial; }

    // Get the current edge the iterator is at.
    EdgeT * getCurrent() const { return e; }

    // Check if the iterator is at the initial position without ever having moved away from it, or not.
    bool isAtInitial() const { return first; }

  private:
    EdgeT * initial;  // The edge used to construct the iterator.
    EdgeT * e;        // The current edge.
    bool first;       // Is this the the first time we're seeing the edge the iterator was initialized with?

}; // class FwdIterator

} // namespace DCELInternal

/**
 * Halfedge of DcelMesh.
 *
 * Adapted from: DcelHalfedge class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DcelHalfedge : public AttributedObject<HalfedgeAttribute>
{
  public:
    typedef DcelMesh  <VertexAttribute, HalfedgeAttribute, FaceAttribute> Mesh;    ///< Parent mesh class.
    typedef DcelFace  <VertexAttribute, HalfedgeAttribute, FaceAttribute> Face;    ///< Face of the mesh.
    typedef DcelVertex<VertexAttribute, HalfedgeAttribute, FaceAttribute> Vertex;  ///< Vertex of the mesh.

  private:
    /** "Dereference" a halfedge to obtain the associated face. */
    struct FaceDeref { Face * operator()(DcelHalfedge const * e) const { return const_cast<Face *>(e->face); } };

    /** Move to the halfedge that links to the next non-null face associated with this edge. */
    struct FaceIncrement { void operator()(DcelHalfedge const ** e) const { if ((*e)->twin()->face) *e = (*e)->twin(); } };

  public:
    /** An iterator over the faces incident on the edge. */
    typedef DCELInternal::FwdIterator<DcelHalfedge, Face, FaceDeref, FaceIncrement> FaceIterator;

    /** A const iterator over the faces incident on the edge. */
    typedef DCELInternal::FwdIterator<DcelHalfedge const, Face const, FaceDeref, FaceIncrement> FaceConstIterator;

    /** Default constructor. */
    DcelHalfedge(intx index_ = -1) : index(index_), twin_he(nullptr), next_he(nullptr), face(nullptr), origin(nullptr), bits(0)
    {}

    /** Get the vertex from which this halfedge originates. */
    Vertex const * getOrigin() const
    {
      debugAssertM(origin, "DcelHalfedge: Halfedge has no origin");
      return origin;
    }

    /** Get the vertex from which this halfedge originates. */
    Vertex * getOrigin()
    {
      debugAssertM(origin, "DcelHalfedge: Halfedge has no origin");
      return origin;
    }

    /** Get the vertex at which this halfedge ends. */
    Vertex const * getEnd() const { return twin()->getOrigin(); }

    /** Get the vertex at which this halfedge ends. */
    Vertex * getEnd() { return twin()->getOrigin(); }

    /** Get an endpoint of the edge. \a i = 0 returns the origin and \a i = 1 the end. */
    Vertex const * getEndpoint(int i) const
    {
      debugAssertM(i == 0 || i == 1, "DcelHalfedge: Invalid endpoint index");
      return i == 0 ? origin : getEnd();
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a i = 1 the second. */
    Vertex * getEndpoint(int i)
    {
      debugAssertM(i == 0 || i == 1, "DcelHalfedge: Invalid endpoint index");
      return i == 0 ? origin : getEnd();
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex const * getOtherEndpoint(Vertex const * endpoint) const
    {
      debugAssertM(hasEndpoint(endpoint), "DcelHalfedge: Vertex is not an endpoint of the edge");
      return origin == endpoint ? getEnd() : origin;
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex * getOtherEndpoint(Vertex const * endpoint)
    {
      debugAssertM(hasEndpoint(endpoint), "DcelHalfedge: Vertex is not an endpoint of the edge");
      return origin == endpoint ? getEnd() : origin;
    }

    /** Get the index (0 or 1) of an endpoint given a pointer to it, or a negative value if the neither endpoint matches. */
    int getEndpointIndex(Vertex const * endpoint) const
    {
      return origin == endpoint ? 0 : (getEnd() == endpoint ? 1 : -1);
    }

    /** Check if the edge has a given vertex as an endpoint. */
    bool hasEndpoint(Vertex const * v) const { return origin == v || getEnd() == v; }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentFace(Face const * f) const
    {
      return face == f || twin()->face == f;
    }

    /** Get the next halfedge around the face. */
    DcelHalfedge const * next() const
    {
      debugAssertM(next_he, "DcelHalfedge: Halfedge has no successor");
      return next_he;
    }

    /** Get the next halfedge around the face. */
    DcelHalfedge * next()
    {
      debugAssertM(next_he, "DcelHalfedge: Halfedge has no successor");
      return next_he;
    }

    /** Get the halfedge between the same two vertices, but in the opposite direction. */
    DcelHalfedge const * twin() const
    {
      debugAssertM(twin_he, "DcelHalfedge: Halfedge has no twin_he");
      return twin_he;
    }

    /** Get the next halfedge around the face. */
    DcelHalfedge * twin()
    {
      debugAssertM(twin_he, "DcelHalfedge: Halfedge has no twin_he");
      return twin_he;
    }

    /**
     * Get the next halfedge around the originating vertex.
     */
    DcelHalfedge const * nextAroundOrigin() const { return twin()->next(); }

    /**
     * Get the next halfedge around the originating vertex.
     */
    DcelHalfedge * nextAroundOrigin() { return twin()->next(); }

    /**
     * Get the previous halfedge around the originating vertex.
     */
    DcelHalfedge const * prevAroundOrigin() const { return const_cast<DcelHalfedge *>(this)->prevAroundOrigin(); }

    /**
     * Get the previous halfedge around the originating vertex.
     */
    DcelHalfedge * prevAroundOrigin()
    {
      DcelHalfedge * rval = this;
      DcelHalfedge * next_around_vertex = twin()->next();

      while (next_around_vertex != this)
      {
        rval = next_around_vertex;
        next_around_vertex = next_around_vertex->twin()->next();
      }

      return rval;
    }

    /** Get the number of faces incident on the edge. */
    intx numFaces() const { return (face ? 1 : 0) + (twin()->face ? 1 : 0); }

    /** Get the face adjoining this halfedge (or null if this is a border halfedge). */
    Face const * getFace() const { return face; }

    /** Get the face adjoining this halfedge (or null if this is a border halfedge). */
    Face * getFace() { return face; }

    /** Get an iterator pointing to the first face incident on the bidirectional edge (this edge plus its twin). */
    FaceConstIterator facesBegin() const
    {
      return FaceConstIterator(const_cast<DcelHalfedge *>(this)->facesBegin());
    }

    /** Get an iterator pointing to the first face incident on the bidirectional edge (this edge plus its twin). */
    FaceIterator facesBegin()
    {
      if (face) return FaceIterator(this);
      else if (twin()->face) return FaceIterator(twin());
      else return FaceIterator();
    }

    /**
     * Get an iterator pointing to one position beyond the last face incident on the bidirectional edge (this edge plus its
     * twin).
     */
    FaceConstIterator facesEnd() const
    {
      return FaceConstIterator(const_cast<DcelHalfedge *>(this)->facesEnd());
    }

    /**
     * Get an iterator pointing to one position beyond the last face incident on the bidirectional edge (this edge plus its
     * twin).
     */
    FaceIterator facesEnd()
    {
      if (face) return FaceIterator(this, false);
      else if (twin()->face) return FaceIterator(twin(), false);
      else return FaceIterator();
    }

    /** Check if this is a boundary halfedge, i.e. if its face pointer is null. */
    bool isBoundaryHalfedge() const { return !face; }

    /** Check if this is a boundary edge, i.e. if either this halfedge or its twin has a null face pointer. */
    bool isBoundaryEdge() const { return !face || !twin()->face; }

    /** Check if one or more marker bits are set. */
    bool areBitsSet(unsigned char mask) const { return ((bits & mask) == mask); };

    /** Set one or more marker bits on or off. */
    void setBits(unsigned char mask, bool value) { bits = (value ? (bits | mask) : (bits & ~mask)); };

    /** Set all marker bits to off. */
    void clearAllBits() { bits = 0; };

  private:
    template <typename _VertexAttribute, typename _HalfedgeAttribute, typename _FaceAttribute> friend class DcelMesh;

    intx index;
    DcelHalfedge * twin_he;
    DcelHalfedge * next_he;
    Face * face;
    Vertex * origin;
    unsigned char bits;

}; // class DcelHalfedge

} // namespace Graphics
} // namespace Thea

#endif
