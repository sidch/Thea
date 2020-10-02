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

#ifndef __Thea_Graphics_DcelVertex_hpp__
#define __Thea_Graphics_DcelVertex_hpp__

#include "../Common.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../AttributedObject.hpp"
#include "DcelHalfedge.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DcelMesh;

/**
 * Vertex of DcelMesh.
 *
 * Adapted from: DcelVertex class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DcelVertex
: public PositionAttribute<Vector3>,
  public NormalAttribute<Vector3>,
  public AttributedObject<VertexAttribute>
{
  public:
    typedef DcelMesh    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Mesh;      ///< Parent mesh class.
    typedef DcelHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute>  Halfedge;  ///< Halfedge of the mesh.
    typedef Halfedge                                                         Edge;      /**< Typedef for interoperability with
                                                                                             GeneralMesh. */
    typedef DcelFace    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Face;      ///< Face of the mesh.

  private:
    typedef PositionAttribute<Vector3>  PositionBaseType;
    typedef NormalAttribute<Vector3>    NormalBaseType;

  private:
    /** "Dereference" a halfedge and return the edge itself. */
    struct EdgeDeref { Edge * operator()(Halfedge const * e) const { return const_cast<Edge *>(e); } };

    /** Move to the halfedge that links to the next non-null face associated with this vertex. */
    struct EdgeIncrement { void operator()(Halfedge const ** e) const { (*e) = (*e)->nextAroundOrigin(); } };

    /** "Dereference" a halfedge to obtain the associated face. */
    struct FaceDeref { Face * operator()(Halfedge const * e) const { return const_cast<Face *>(e->getFace()); } };

    /** Move to the halfedge that links to the next non-null face associated with this vertex. */
    struct FaceIncrement
    {
      // Termination assumes the iteration started with at least one non-null face
      void operator()(Halfedge const ** e) const { do { (*e) = (*e)->nextAroundOrigin(); } while (!(*e)->getFace()); }
    };

  public:
    /** An iterator over the edges incident on the vertex. */
    typedef DCELInternal::FwdIterator<Halfedge, Edge, EdgeDeref, EdgeIncrement> EdgeIterator;

    /** A const iterator over the edges incident on the vertex. */
    typedef DCELInternal::FwdIterator<Halfedge const, Edge const, EdgeDeref, EdgeIncrement> EdgeConstIterator;

    /** An iterator over the faces incident on the vertex. */
    typedef DCELInternal::FwdIterator<Halfedge, Face, FaceDeref, FaceIncrement> FaceIterator;

    /** A const iterator over the faces incident on the vertex. */
    typedef DCELInternal::FwdIterator<Halfedge const, Face const, FaceDeref, FaceIncrement> FaceConstIterator;

    /** Default constructor. */
    DcelVertex()
    : NormalBaseType(Vector3::Zero()), leaving(nullptr), index(-1), has_precomputed_normal(false),
      normal_normalization_factor(0)
    {}

    /** Sets the vertex to have a location. */
    explicit DcelVertex(Vector3 const & p)
    : PositionBaseType(p), NormalBaseType(Vector3::Zero()), leaving(nullptr), index(-1), has_precomputed_normal(false),
      normal_normalization_factor(0)
    {}

    /**
     * Sets the vertex to have a location and a precomputed normal. Adding faces will <b>not</b> change the vertex normal if
     * this constructor is used.
     */
    DcelVertex(Vector3 const & p, Vector3 const & n)
    : PositionBaseType(p), NormalBaseType(n), leaving(nullptr), index(-1), has_precomputed_normal(true),
      normal_normalization_factor(0)
    {}

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge const * getEdgeTo(DcelVertex const * v) const  { return const_cast<DcelVertex *>(this)->getEdgeTo(v); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge * getEdgeTo(DcelVertex const * v)
    {
      Halfedge const * rval = nullptr;

      if (leaving)
      {
        if (leaving->twin()->getOrigin() == v)
          rval = leaving;
        else
        {
          Halfedge const * test = leaving->twin()->next();
          while (rval == nullptr && test != leaving)
          {
            Halfedge const * twin = test->twin();
            if (twin->getOrigin() == v)
              rval = test;
            else
              test = twin->next();
          }
        }
      }

      return const_cast<Halfedge *>(rval);
    }

    /** Get a canonical halfedge leaving this vertex, or null if the vertex is isolated. */
    Halfedge const * getHalfedge() const { return leaving; }

    /** Get a canonical halfedge leaving this vertex, or null if the vertex is isolated. */
    Halfedge * getHalfedge() { return leaving; }

    /** Get an iterator pointing to the first edge incident on the vertex. */
    EdgeConstIterator edgesBegin() const { return EdgeConstIterator(leaving); }

    /** Get an iterator pointing to the first edge incident on the vertex. */
    EdgeIterator edgesBegin() { return EdgeIterator(leaving); }

    /** Get an iterator pointing to the position beyond the last edge incident on the vertex. */
    EdgeConstIterator edgesEnd() const { return EdgeConstIterator(leaving, false); }

    /** Get an iterator pointing to the position beyond the last edge incident on the vertex. */
    EdgeIterator edgesEnd() { return EdgeIterator(leaving, false); }

    /** Get an iterator pointing to the first face incident on the vertex. */
    FaceConstIterator facesBegin() const { return FaceConstIterator(leaving); }

    /** Get an iterator pointing to the first face incident on the vertex. */
    FaceIterator facesBegin() { return FaceIterator(leaving); }

    /** Get an iterator pointing to the position beyond the last face incident on the vertex. */
    FaceConstIterator facesEnd() const { return FaceConstIterator(leaving, false); }

    /** Get an iterator pointing to the position beyond the last face incident on the vertex. */
    FaceIterator facesEnd() { return FaceIterator(leaving, false); }

    /** Check if the vertex lies on a mesh boundary. */
    bool isBoundaryVertex() const
    {
      // Cycle round the halfedges emanating from this vertex, looking for one on the boundary
      if (leaving)
      {
        Halfedge * e = leaving;
        do
        {
          if (e->isBoundaryEdge())
            return true;

          e = e->nextAroundOrigin();

        } while (e != leaving);

        return false;
      }
      else
        return true;
    }

    /** Get the index of the vertex, typically in the source file (or negative if unindexed). */
    intx getIndex() const { return index; }

    /** Set the index of the vertex, typically from the source file (or negative if unindexed). */
    void setIndex(intx index_) { index = index_; }

  private:
    template <typename V, typename H, typename F> friend class DcelMesh;

    /**
     * Add normal information from a new face at this vertex.
     *
     * @param n Unit (or weighted) normal of the new face.
     */
    void addFaceNormal(Vector3 const & n)
    {
      if (!has_precomputed_normal)
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() + n;
        normal_normalization_factor = sum_normals.norm();
        if (normal_normalization_factor < 1e-20)
          setNormal(Vector3::Zero());
        else
          setNormal(sum_normals / normal_normalization_factor);
      }
    }

    /** Estimate the result of addFaceNormal() without actually updating any data. */
    Vector3 estimateUpdatedNormal(Vector3 const & new_normal) const
    {
      if (has_precomputed_normal)
        return getNormal();  // it won't change
      else
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() + new_normal;
        return sum_normals.normalized();
      }
    }

    /** Get the index of the vertex in a GPU array. */
    uint32 getPackingIndex() const { return packing_index; }

    /** Set the index of the vertex in a GPU array. */
    void setPackingIndex(uint32 packing_index_) const { packing_index = packing_index_; }

    Halfedge * leaving;
    intx index;
    bool has_precomputed_normal;
    float normal_normalization_factor;
    mutable uint32 packing_index;

}; // class DcelVertex

} // namespace Graphics

namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <typename VT, typename HT, typename FT>
class IsPointN<Graphics::DcelVertex<VT, HT, FT>, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <typename VT, typename HT, typename FT>
class PointTraitsN<Graphics::DcelVertex<VT, HT, FT>, 3>
{
  public:
    static Vector3 const & getPosition(Graphics::DcelVertex<VT, HT, FT> const & t) { return t.getPosition(); }
};

} // namespace Algorithms

} // namespace Thea

#endif
