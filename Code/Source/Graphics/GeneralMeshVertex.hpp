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

#ifndef __Thea_Graphics_GeneralMeshVertex_hpp__
#define __Thea_Graphics_GeneralMeshVertex_hpp__

#include "../Common.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../AttributedObject.hpp"
#include "../List.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {

namespace Graphics {

// Forward declarations
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshEdge;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshFace;

/** Vertex of GeneralMesh. */
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class /* THEA_API */ GeneralMeshVertex
: public PositionAttribute<Vector3>,
  public NormalAttribute<Vector3>,
  public AttributedObject<VertexAttributeT>
{
  public:
    typedef GeneralMeshEdge<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;  ///< Edge of the mesh.
    typedef GeneralMeshFace<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;  ///< Face of the mesh.

  private:
    typedef PositionAttribute<Vector3>  PositionBaseType;
    typedef NormalAttribute<Vector3>    NormalBaseType;

    typedef TheaList< Edge *, AllocatorT<Edge *> > EdgeList;
    typedef TheaList< Face *, AllocatorT<Face *> > FaceList;

  public:
    typedef typename EdgeList::iterator        EdgeIterator;       ///< Iterator over edges.
    typedef typename EdgeList::const_iterator  EdgeConstIterator;  ///< Const iterator over edges.
    typedef typename FaceList::iterator        FaceIterator;       ///< Iterator over faces.
    typedef typename FaceList::const_iterator  FaceConstIterator;  ///< Const iterator over faces.

    /** Default constructor. */
    GeneralMeshVertex()
    : NormalBaseType(Vector3::zero()), has_precomputed_normal(false), normal_normalization_factor(0), marked(false) {}

    /** Sets the vertex to have a location. */
    explicit GeneralMeshVertex(Vector3 const & p)
    : PositionBaseType(p), NormalBaseType(Vector3::zero()), has_precomputed_normal(false), normal_normalization_factor(0),
      marked(false)
    {}

    /**
     * Sets the vertex to have a location and a precomputed normal. Adding faces will <b>not</b> change the vertex normal if
     * this constructor is used.
     */
    GeneralMeshVertex(Vector3 const & p, Vector3 const & n)
    : PositionBaseType(p), NormalBaseType(n), has_precomputed_normal(true), normal_normalization_factor(0), marked(false)
    {}

    /**
     * Get the number of edges incident on the vertex. Equivalent to degree().
     *
     * @see degree();
     */
    int numEdges() const { return (int)edges.size(); }

    /**
     * Get the degree of the vertex, i.e. number of edges incident on it. Equivalent to numEdges().
     *
     * @see numEdges();
     */
    int degree() const { return (int)edges.size(); }

    /** Get the number of faces incident on the vertex. */
    int numFaces() const { return (int)faces.size(); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Edge const * getEdgeTo(GeneralMeshVertex const * v) const
    { return const_cast<GeneralMeshVertex *>(this)->getEdgeTo(v); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Edge * getEdgeTo(GeneralMeshVertex const * v)
    {
      if (v == this) return NULL;

      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
      {
        Edge * e = *ei;
        if (e->hasEndpoint(v)) return e;
      }

      return NULL;
    }

    /** Check if the vertex is adjacent to a given edge. */
    bool hasEdgeTo(GeneralMeshVertex const * v) const { return getEdgeTo(v) != NULL; }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentEdge(Edge const * edge) const
    {
      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if (*ei == edge)
          return true;

      return false;
    }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentFace(Face const * face) const
    {
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (*fi == face)
          return true;

      return false;
    }

    /** Get an iterator pointing to the first edge. */
    EdgeConstIterator edgesBegin() const { return edges.begin(); }

    /** Get an iterator pointing to the first edge. */
    EdgeIterator edgesBegin() { return edges.begin(); }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeConstIterator edgesEnd() const { return edges.end(); }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeIterator edgesEnd() { return edges.end(); }

    /** Get an iterator pointing to the first face. */
    FaceConstIterator facesBegin() const { return faces.begin(); }

    /** Get an iterator pointing to the first face. */
    FaceIterator facesBegin() { return faces.begin(); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceConstIterator facesEnd() const { return faces.end(); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceIterator facesEnd() { return faces.end(); }

    /** Check if the vertex lies on a mesh boundary. */
    bool isBoundary() const
    {
      if (edges.empty()) return true;

      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if ((*ei)->isBoundary()) return true;

      return false;
    }

    /** Check if the vertex has a precomputed normal. */
    bool hasPrecomputedNormal() const { return has_precomputed_normal; }

    /**
     * Update the vertex normal by recomputing it from face data. Useful when the faces have been modified externally. Destroys
     * any prior precomputed normal.
     */
    void updateNormal()
    {
      if (!faces.empty())
      {
        Vector3 sum_normals = Vector3::zero();
        for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
          sum_normals += (*fi)->getNormal();  // weight by face area?

        normal_normalization_factor = sum_normals.length();
        setNormal(normal_normalization_factor < 1e-20f ? Vector3::zero() : sum_normals / normal_normalization_factor);
      }
      else
      {
        setNormal(Vector3::zero());
        normal_normalization_factor = 0;
      }

      has_precomputed_normal = false;
    }

  private:
    template <typename V, typename E, typename F, template <typename T> class A> friend class GeneralMesh;

    /** Add a reference to an edge incident at this vertex. */
    void addEdge(Edge * edge) { edges.push_back(edge); }

    /** Remove a reference to an edge incident at this vertex. */
    EdgeIterator removeEdge(EdgeIterator loc) { return edges.erase(loc); }

    /** Remove a reference to an edge incident at this vertex. */
    void removeEdge(Edge * edge)
    {
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); )
      {
        if (*ei == edge)
          ei = edges.erase(ei);
        else
          ++ei;
      }
    }

    /** Replace all references to an edge with references to another edge. */
    void replaceEdge(Edge * old_edge, Edge * new_edge)
    {
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if (*ei == old_edge)
          *ei = new_edge;
    }

    /** Add a reference to a face incident at this vertex. */
    void addFace(Face * face, bool update_normal = true)
    {
      faces.push_back(face);
      if (update_normal && !has_precomputed_normal) addFaceNormal(face->getNormal());
    }

    /** Remove all references to a face incident at this vertex. */
    void removeFace(Face * face)
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); )
      {
        if (*fi == face)
        {
          fi = faces.erase(fi);
          if (!has_precomputed_normal) removeFaceNormal(face->getNormal());

          // Keep going, just in case the face somehow got added twice
        }
        else
          ++fi;
      }
    }

    /** Replace all references to a face with references to another face. */
    void replaceFace(Face * old_face, Face * new_face)
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (*fi == old_face)
          *fi = new_face;
    }

    /**
     * Add normal information from a new face at this vertex, unless the vertex has a precomputed normal.
     *
     * @param n Unit (or weighted) normal of the new face.
     */
    void addFaceNormal(Vector3 const & n)
    {
      if (!has_precomputed_normal)
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() + n;
        normal_normalization_factor = sum_normals.length();
        setNormal(normal_normalization_factor < 1e-20f ? Vector3::zero() : sum_normals / normal_normalization_factor);
      }
    }

    /**
     * Remove normal information from a new face at this vertex, unless the vertex has a precomputed normal.
     *
     * @param n Unit (or weighted) normal of the face to be removed.
     */
    void removeFaceNormal(Vector3 const & n)
    {
      if (!has_precomputed_normal)
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() - n;
        normal_normalization_factor = sum_normals.length();
        setNormal(normal_normalization_factor < 1e-20f ? Vector3::zero() : sum_normals / normal_normalization_factor);
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
        return sum_normals.unit();
      }
    }

    /** Mark the vertex. */
    void mark()
    {
      marked = true;
    }

    /** Unmark the vertex. */
    void unmark()
    {
      marked = false;
    }

    /** Check if the vertex is marked. */
    bool isMarked() const { return marked; }

    /** Get the index of the vertex. */
    uint32 getIndex() const { return index; }

    /** Set the index of the vertex. */
    void setIndex(uint32 index_) { index = index_; }

    /** Make an exact copy of the vertex. */
    void copyTo(GeneralMeshVertex & dst,
                TheaUnorderedMap<GeneralMeshVertex const *, GeneralMeshVertex *> const & vertex_map,
                TheaUnorderedMap<Edge const *, Edge *> const & edge_map,
                TheaUnorderedMap<Face const *, Face *> const & face_map) const
    {
      dst.setPosition(this->getPosition());
      dst.setNormal(this->getNormal());
      dst.setAttr(this->attr());  // assume attributes can be simply copied

      dst.edges.resize(edges.size());
      EdgeConstIterator ei = edges.begin();
      EdgeIterator dei = dst.edges.begin();
      for ( ; ei != edges.end(); ++ei, ++dei)
        *dei = edge_map.find(*ei)->second;  // assume it always exists

      dst.faces.resize(faces.size());
      FaceConstIterator fi = faces.begin();
      FaceIterator dfi = dst.faces.begin();
      for ( ; fi != faces.end(); ++fi, ++dfi)
        *dfi = face_map.find(*fi)->second;  // assume it always exists

      dst.has_precomputed_normal = has_precomputed_normal;
      dst.normal_normalization_factor = normal_normalization_factor;
      dst.index = index;
      dst.marked = marked;
    }

    EdgeList edges;
    FaceList faces;
    bool has_precomputed_normal;
    float normal_normalization_factor;
    uint32 index;
    bool marked;

}; // class GeneralMeshVertex

} // namespace Graphics

namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <typename VT, typename ET, typename FT, template <typename T> class AT>
class IsPointN<Graphics::GeneralMeshVertex<VT, ET, FT, AT>, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <typename VT, typename ET, typename FT, template <typename T> class AT>
class PointTraitsN<Graphics::GeneralMeshVertex<VT, ET, FT, AT>, 3>
{
  public:
    static Vector3 const & getPosition(Graphics::GeneralMeshVertex<VT, ET, FT, AT> const & t) { return t.getPosition(); }
};

} // namespace Algorithms

} // namespace Thea

#endif
