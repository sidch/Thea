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

#ifndef __Thea_Graphics_GeneralMeshFace_hpp__
#define __Thea_Graphics_GeneralMeshFace_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "../List.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshVertex;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshEdge;

/** Face of GeneralMesh. */
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class /* THEA_API */ GeneralMeshFace : public NormalAttribute<Vector3>, public AttributedObject<FaceAttributeT>
{
  public:
    typedef GeneralMeshVertex<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;  ///< Vertex of the mesh.
    typedef GeneralMeshEdge  <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;    ///< Edge of the mesh.

  private:
    typedef NormalAttribute<Vector3> NormalBaseType;

    typedef TheaList< Vertex *, AllocatorT<Vertex *> >  VertexList;
    typedef TheaList< Edge *,   AllocatorT<Edge *>   >  EdgeList;

  public:
    typedef typename VertexList::iterator                VertexIterator;              ///< Iterator over vertices.
    typedef typename VertexList::const_iterator          VertexConstIterator;         ///< Const iterator over vertices.
    typedef typename VertexList::reverse_iterator        VertexReverseIterator;       ///< Reverse iterator over vertices.
    typedef typename VertexList::const_reverse_iterator  VertexConstReverseIterator;  ///< Const reverse iterator over vertices.
    typedef typename EdgeList::iterator                  EdgeIterator;                ///< Iterator over edges.
    typedef typename EdgeList::const_iterator            EdgeConstIterator;           ///< Const iterator over edges.
    typedef typename EdgeList::reverse_iterator          EdgeReverseIterator;         ///< Reverse iterator over edges.
    typedef typename EdgeList::const_reverse_iterator    EdgeConstReverseIterator;    ///< Const reverse iterator over edges.

    /** Construct with the given normal. */
    GeneralMeshFace(Vector3 const & normal = Vector3::zero()) : NormalBaseType(normal), marked(false) {}

    /** Check if the face has a given vertex. */
    bool hasVertex(Vertex const * vertex) const
    {
      for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        if (*vi == vertex)
          return true;

      return false;
    }

    /** Check if the face has a given edge. */
    bool hasEdge(Edge const * edge) const
    {
      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if (*ei == edge)
          return true;

      return false;
    }

    /** Get the predecessor of a vertex around the face. Assumes the iterator points to a valid vertex of the face. */
    Vertex const * getPredecessor(VertexConstIterator vertex) const
    {
      return const_cast<GeneralMeshFace *>(this)->getPredecessor(vertex);
    }

    /** Get the predecessor of a vertex around the face. Assumes the iterator points to a valid vertex of the face. */
    Vertex * getPredecessor(VertexConstIterator vertex)
    {
      if (vertices.empty())
        return NULL;
      else if (vertex == vertices.begin())
        return vertices.back();
      else
      {
        --vertex;
        return *vertex;
      }
    }

    /** Get the predecessor of a vertex around the face. Returns null if the vertex does not belong to the face. */
    Vertex const * getPredecessor(Vertex const * vertex) const
    {
      return const_cast<GeneralMeshFace *>(this)->getPredecessor(vertex);
    }

    /** Get the predecessor of a vertex around the face. Returns null if the vertex does not belong to the face. */
    Vertex * getPredecessor(Vertex const * vertex)
    {
      for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        if (*vi == vertex)
          return getPredecessor(vi);

      return NULL;
    }

    /** Get the successor of a vertex around the face. Assumes the iterator points to a valid vertex of the face. */
    Vertex const * getSuccessor(VertexConstIterator vertex) const
    {
      return const_cast<GeneralMeshFace *>(this)->getSuccessor(vertex);
    }

    /** Get the successor of a vertex around the face. Assumes the iterator points to a valid vertex of the face. */
    Vertex * getSuccessor(VertexConstIterator vertex)
    {
      if (vertices.empty())
        return NULL;
      else
      {
        ++vertex;
        return vertex == vertices.end() ? vertices.front() : *vertex;
      }
    }

    /** Get the successor of a vertex around the face. Returns null if the vertex does not belong to the face. */
    Vertex const * getSuccessor(Vertex const * vertex) const
    {
      return const_cast<GeneralMeshFace *>(this)->getSuccessor(vertex);
    }

    /** Get the successor of a vertex around the face. Returns null if the vertex does not belong to the face. */
    Vertex * getSuccessor(Vertex const * vertex)
    {
      for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        if (*vi == vertex)
          return getSuccessor(vi);

      return NULL;
    }

    /** Get the predecessor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge const * getPredecessor(EdgeConstIterator edge) const
    {
      return const_cast<GeneralMeshFace *>(this)->getPredecessor(edge);
    }

    /** Get the predecessor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge * getPredecessor(EdgeConstIterator edge)
    {
      if (edges.empty())
        return NULL;
      else if (edge == edges.begin())
        return edges.back();
      else
      {
        --edge;
        return *edge;
      }
    }

    /** Get the predecessor of an edge around the face. Returns null if the edge does not belong to the face. */
    Edge const * getPredecessor(Edge const * edge) const { return const_cast<GeneralMeshFace *>(this)->getPredecessor(edge); }

    /** Get the predecessor of an edge around the face. Returns null if the edge does not belong to the face. */
    Edge * getPredecessor(Edge const * edge)
    {
      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if (*ei == edge)
          return getPredecessor(ei);

      return NULL;
    }

    /** Get the successor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge const * getSuccessor(EdgeConstIterator edge) const { return const_cast<GeneralMeshFace *>(this)->getSuccessor(edge); }

    /** Get the successor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge * getSuccessor(EdgeConstIterator edge)
    {
      if (edges.empty())
        return NULL;
      else
      {
        ++edge;
        return edge == edges.end() ? edges.front() : *edge;
      }
    }

    /** Get the successor of an edge around the face. Returns null if the edge does not belong to the face. */
    Edge const * getSuccessor(Edge const * edge) const { return const_cast<GeneralMeshFace *>(this)->getSuccessor(edge); }

    /** Get the successor of an edge around the face. Returns null if the edge does not belong to the face. */
    Edge * getSuccessor(Edge const * edge)
    {
      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if (*ei == edge)
          return getSuccessor(ei);

      return NULL;
    }

    /** Get an iterator pointing to the first vertex. */
    VertexConstIterator verticesBegin() const { return vertices.begin(); }

    /** Get an iterator pointing to the first vertex. */
    VertexIterator verticesBegin() { return vertices.begin(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexConstIterator verticesEnd() const { return vertices.end(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexIterator verticesEnd() { return vertices.end(); }

    /** Get a reverse iterator pointing to the last vertex. */
    VertexConstReverseIterator verticesReverseBegin() const { return vertices.rbegin(); }

    /** Get a reverse iterator pointing to the last vertex. */
    VertexReverseIterator verticesReverseBegin() { return vertices.rbegin(); }

    /** Get a reverse iterator pointing to the position before the first vertex. */
    VertexConstReverseIterator verticesReverseEnd() const { return vertices.rend(); }

    /** Get a reverse iterator pointing to the position before the first vertex. */
    VertexReverseIterator verticesReverseEnd() { return vertices.rend(); }

    /** Get an iterator pointing to the first edge. */
    EdgeConstIterator edgesBegin() const { return edges.begin(); }

    /** Get an iterator pointing to the first edge. */
    EdgeIterator edgesBegin() { return edges.begin(); }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeConstIterator edgesEnd() const { return edges.end(); }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeIterator edgesEnd() { return edges.end(); }

    /** Get a reverse iterator pointing to the last edge. */
    EdgeConstReverseIterator edgesReverseBegin() const { return edges.rbegin(); }

    /** Get a reverse iterator pointing to the last edge. */
    EdgeReverseIterator edgesReverseBegin() { return edges.rbegin(); }

    /** Get a reverse iterator pointing to the position before the first edge. */
    EdgeConstReverseIterator edgesReverseEnd() const { return edges.rend(); }

    /** Get a reverse iterator pointing to the position before the first edge. */
    EdgeReverseIterator edgesReverseEnd() { return edges.rend(); }

    /** Get the number of vertices of the face. */
    int numVertices() const
    {
      debugAssertM(vertices.size() == edges.size(), "GeneralMeshFace: Numbers of edges != number of vertices");
      return (int)vertices.size();
    }

    /** Get the number of edges bordering the face. */
    int numEdges() const
    {
      debugAssertM(vertices.size() == edges.size(), "GeneralMeshFace: Numbers of edges != number of vertices");
      return (int)edges.size();
    }

    /** Check if the face is a triangle. */
    bool isTriangle() const { return  vertices.size() == 3; }

    /** Check if the face is a quad. */
    bool isQuad() const { return vertices.size() == 4; }

    /** Reverse the order in which vertices and edges wind around the face. The face normal is <b>not</b> modified. */
    void reverseWinding()
    {
      vertices.reverse();
      edges.reverse();
    }

    /** Update the face normal by recomputing it from vertex data. */
    void updateNormal()
    {
      // Assume the face is planar.
      VertexConstIterator vi2 = vertices.begin();
      VertexConstIterator vi0 = vi2++;
      VertexConstIterator vi1 = vi2++;

      if (vertices.size() > 3)
      {
        // vi1 might be a concave corner -- we need to add up the cross products at all vertices
        Vector3 sum_cross = Vector3::zero();
        for ( ; vi0 != vertices.end(); ++vi0, ++vi1, ++vi2)
        {
          if (vi1 == vertices.end()) vi1 = vertices.begin();
          if (vi2 == vertices.end()) vi2 = vertices.begin();

          Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
          Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
          sum_cross += e2.cross(e1);
        }

        setNormal(sum_cross.unit());
      }
      else
      {
        Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
        Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
        setNormal(e2.cross(e1).unit());  // counter-clockwise
      }
    }

    /** Compute the centroid of the face. */
    Vector3 centroid() const
    {
      Vector3 c(0, 0, 0);
      if (!vertices.empty())
      {
        for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
          c += (*vi)->getPosition();

        c /= vertices.size();
      }

      return c;
    }

    /**
     * Test if the face contains a point (which is assumed to lie on the plane of the face -- for efficiency the function does
     * <b>not</b> explicitly verify that this holds).
     */
    bool contains(Vector3 const & p) const
    {
      if (vertices.empty()) return false;

      // Generate a ray for the even-odd test, from p to the midpoint of the first halfedge. Ignore degenerate situations for
      // now.
      VertexConstIterator vi    =  verticesBegin();
      VertexConstIterator last  =  vi++;
      Vector3 u = 0.5 * ((*last)->getPosition() + (*vi)->getPosition()) - p;

      long count = 1;  // first halfedge is obviously intersected, since we generated the ray through its midpoint
      for ( ; last != verticesEnd(); ++vi)
      {
        if (vi == verticesEnd()) vi = verticesBegin();

        Vector3 v0 = (*last)->getPosition() - p;
        Vector3 v1 = (*vi)->getPosition()   - p;

        // If winding order is: vector to first vertex, ray, vector to second vertex, then intersects
        Vector3 c0 = v0.cross(u);
        Vector3 c1 = u.cross(v1);
        if (c0.dot(c1) > 0)  // intersects, now check forward or reverse
        {
          // Forward if the vector to the point nearest to p on the line containing the edge makes an acute angle with u.
          //
          // The point p' on line v + t * e closest to point p is v + t0 * e, where t0 = e.dot(p - v) / e.dot(e)
          // (see www.geometrictools.com/Documentation/DistancePointLine.pdf).
          //
          // We translate p to the origin for simpler computations.
          Vector3 edge = v1 - v0;
          Real t0 = -edge.dot(v0) / edge.dot(edge);
          Vector3 u0 = v0 + t0 * edge;

          if (u0.dot(u) > 0)
            count++;
        }

        last = vi;
      }

      return (count % 2 == 1);
    }

  private:
    template <typename V, typename E, typename F, template <typename T> class A> friend class GeneralMesh;
    template <typename V, typename E, typename F, template <typename T> class A> friend class GeneralMeshEdge;

    /** Add a reference to a vertex of this face. */
    void addVertex(Vertex * vertex) { vertices.push_back(vertex); }

    /** Remove a reference to a vertex. */
    VertexIterator removeVertex(VertexIterator loc) { return vertices.erase(loc); }

    /** Remove all references to a vertex. */
    void removeVertex(Vertex * vertex)
    {
      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); )
      {
        if (*vi == vertex)
          vi = vertices.erase(vi);
        else
          ++vi;
      }
    }

    /** Replace all references to an vertex with references to another vertex. */
    void replaceVertex(Vertex * old_vertex, Vertex * new_vertex)
    {
      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        if (*vi == old_vertex)
          *vi = new_vertex;
    }

    /** Add a reference to an edge of this face. */
    void addEdge(Edge * edge) { edges.push_back(edge); }

    /** Remove a reference to an edge. */
    EdgeIterator removeEdge(EdgeIterator loc) { return edges.erase(loc); }

    /** Remove all references to an edge. */
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

    /** Mark the face. */
    void mark()
    {
      marked = true;
    }

    /** Unmark the face. */
    void unmark()
    {
      marked = false;
    }

    /** Check if the face is marked. */
    bool isMarked() const { return marked; }

    /** Make an exact copy of the face. */
    void copyTo(GeneralMeshFace & dst,
                TheaUnorderedMap<Vertex const *, Vertex *> const & vertex_map,
                TheaUnorderedMap<Edge const *, Edge *> const & edge_map,
                TheaUnorderedMap<GeneralMeshFace const *, GeneralMeshFace *> const & face_map) const
    {
      dst.setNormal(this->getNormal());
      dst.setAttr(this->attr());  // assume attributes can be simply copied

      dst.vertices.resize(vertices.size());
      VertexConstIterator vi = vertices.begin();
      VertexIterator dvi = dst.vertices.begin();
      for ( ; vi != vertices.end(); ++vi, ++dvi)
        *dvi = vertex_map.find(*vi)->second;  // assume it always exists

      dst.edges.resize(edges.size());
      EdgeConstIterator ei = edges.begin();
      EdgeIterator dei = dst.edges.begin();
      for ( ; ei != edges.end(); ++ei, ++dei)
        *dei = edge_map.find(*ei)->second;  // assume it always exists

      dst.marked = marked;
    }

    /** Reset the face to an empty state. */
    void clear()
    {
      vertices.clear();
      edges.clear();
      marked = false;
    }

    VertexList vertices;
    EdgeList edges;
    bool marked;

}; // class GeneralMeshFace

} // namespace Graphics
} // namespace Thea

#endif
