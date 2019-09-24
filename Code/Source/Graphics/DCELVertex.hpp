//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#ifndef __Thea_Graphics_DCELVertex_hpp__
#define __Thea_Graphics_DCELVertex_hpp__

#include "../Common.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../AttributedObject.hpp"
#include "DCELHalfedge.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DCELMesh;

/**
 * Vertex of DCELMesh.
 *
 * Adapted from: DCELVertex class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DCELVertex
: public PositionAttribute<Vector3>,
  public NormalAttribute<Vector3>,
  public AttributedObject<VertexAttribute>
{
  public:
    typedef DCELMesh    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Mesh;      ///< Parent mesh class.
    typedef DCELHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute>  Halfedge;  ///< Halfedge of the mesh.
    typedef Halfedge                                                         Edge;      /**< Typedef for interoperability with
                                                                                             GeneralMesh. */
    typedef DCELFace    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Face;      ///< Face of the mesh.

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
    DCELVertex()
    : NormalBaseType(Vector3::Zero()), leaving(NULL), index(-1), has_precomputed_normal(false), normal_normalization_factor(0)
    {}

    /** Sets the vertex to have a location. */
    explicit DCELVertex(Vector3 const & p)
    : PositionBaseType(p), NormalBaseType(Vector3::Zero()), leaving(NULL), index(-1), has_precomputed_normal(false),
      normal_normalization_factor(0)
    {}

    /**
     * Sets the vertex to have a location and a precomputed normal. Adding faces will <b>not</b> change the vertex normal if
     * this constructor is used.
     */
    DCELVertex(Vector3 const & p, Vector3 const & n)
    : PositionBaseType(p), NormalBaseType(n), leaving(NULL), index(-1), has_precomputed_normal(true),
      normal_normalization_factor(0)
    {}

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge const * getEdgeTo(DCELVertex const * v) const  { return const_cast<DCELVertex *>(this)->getEdgeTo(v); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge * getEdgeTo(DCELVertex const * v)
    {
      Halfedge const * rval = NULL;

      if (leaving)
      {
        if (leaving->twin()->getOrigin() == v)
          rval = leaving;
        else
        {
          Halfedge const * test = leaving->twin()->next();
          while (rval == NULL && test != leaving)
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
    template <typename V, typename H, typename F> friend class DCELMesh;

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
    void setPackingIndex(uint32 packing_index_) { packing_index = packing_index_; }

    Halfedge * leaving;
    intx index;
    bool has_precomputed_normal;
    float normal_normalization_factor;
    uint32 packing_index;

}; // class DCELVertex

} // namespace Graphics

namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <typename VT, typename HT, typename FT>
class IsPointN<Graphics::DCELVertex<VT, HT, FT>, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <typename VT, typename HT, typename FT>
class PointTraitsN<Graphics::DCELVertex<VT, HT, FT>, 3>
{
  public:
    static Vector3 const & getPosition(Graphics::DCELVertex<VT, HT, FT> const & t) { return t.getPosition(); }
};

} // namespace Algorithms

} // namespace Thea

#endif
