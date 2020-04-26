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
class GeneralMesh;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshVertex;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshEdge;

/** Face of GeneralMesh. */
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class /* THEA_API */ GeneralMeshFace : public NormalAttribute<Vector3>, public AttributedObject<FaceAttributeT>
{
  public:
    typedef GeneralMesh      <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Mesh;    ///< Parent mesh class.
    typedef GeneralMeshVertex<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;  ///< Vertex of the mesh.
    typedef GeneralMeshEdge  <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;    ///< Edge of the mesh.

  private:
    typedef NormalAttribute<Vector3> NormalBaseType;

    typedef List< Vertex *, AllocatorT<Vertex *> >  VertexList;
    typedef List< Edge *,   AllocatorT<Edge *>   >  EdgeList;

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
    GeneralMeshFace(Vector3 const & normal = Vector3::Zero()) : NormalBaseType(normal), index(-1), marked(false) {}

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
        return nullptr;
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

      return nullptr;
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
        return nullptr;
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

      return nullptr;
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
        return nullptr;
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

      return nullptr;
    }

    /** Get the successor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge const * getSuccessor(EdgeConstIterator edge) const { return const_cast<GeneralMeshFace *>(this)->getSuccessor(edge); }

    /** Get the successor of an edge around the face. Assumes the iterator points to a valid edge of the face. */
    Edge * getSuccessor(EdgeConstIterator edge)
    {
      if (edges.empty())
        return nullptr;
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

      return nullptr;
    }

    /** Get an iterator pointing to the first vertex of the face. */
    VertexConstIterator verticesBegin() const { return vertices.begin(); }

    /** Get an iterator pointing to the first vertex of the face. */
    VertexIterator verticesBegin() { return vertices.begin(); }

    /** Get an iterator pointing to one position beyond the last vertex of the face. */
    VertexConstIterator verticesEnd() const { return vertices.end(); }

    /** Get an iterator pointing to one position beyond the last vertex of the face. */
    VertexIterator verticesEnd() { return vertices.end(); }

    /** Get a reverse iterator pointing to the last vertex of the face. */
    VertexConstReverseIterator verticesReverseBegin() const { return vertices.rbegin(); }

    /** Get a reverse iterator pointing to the last vertex of the face. */
    VertexReverseIterator verticesReverseBegin() { return vertices.rbegin(); }

    /** Get a reverse iterator pointing to the position before the first vertex of the face. */
    VertexConstReverseIterator verticesReverseEnd() const { return vertices.rend(); }

    /** Get a reverse iterator pointing to the position before the first vertex of the face. */
    VertexReverseIterator verticesReverseEnd() { return vertices.rend(); }

    /** Get an iterator pointing to the first edge of the face. */
    EdgeConstIterator edgesBegin() const { return edges.begin(); }

    /** Get an iterator pointing to the first edge of the face. */
    EdgeIterator edgesBegin() { return edges.begin(); }

    /** Get an iterator pointing to one position beyond the last edge of the face. */
    EdgeConstIterator edgesEnd() const { return edges.end(); }

    /** Get an iterator pointing to one position beyond the last edge of the face. */
    EdgeIterator edgesEnd() { return edges.end(); }

    /** Get a reverse iterator pointing to the last edge of the face. */
    EdgeConstReverseIterator edgesReverseBegin() const { return edges.rbegin(); }

    /** Get a reverse iterator pointing to the last edge of the face. */
    EdgeReverseIterator edgesReverseBegin() { return edges.rbegin(); }

    /** Get a reverse iterator pointing to the position before the first edge of the face. */
    EdgeConstReverseIterator edgesReverseEnd() const { return edges.rend(); }

    /** Get a reverse iterator pointing to the position before the first edge of the face. */
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

    /** Get the index of the face, typically in the source file (or negative if unindexed). */
    intx getIndex() const { return index; }

    /** Set the index of the face, typically from the source file (or negative if unindexed). */
    void setIndex(intx index_) { index = index_; }

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
        Vector3 sum_cross = Vector3::Zero();
        for ( ; vi0 != vertices.end(); ++vi0, ++vi1, ++vi2)
        {
          if (vi1 == vertices.end()) vi1 = vertices.begin();
          if (vi2 == vertices.end()) vi2 = vertices.begin();

          Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
          Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
          sum_cross += e2.cross(e1);
        }

        setNormal(sum_cross.normalized());
      }
      else
      {
        Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
        Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
        setNormal(e2.cross(e1).normalized());  // counter-clockwise
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

      intx count = 1;  // first halfedge is obviously intersected, since we generated the ray through its midpoint
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
                UnorderedMap<Vertex const *, Vertex *> const & vertex_map,
                UnorderedMap<Edge const *, Edge *> const & edge_map,
                UnorderedMap<GeneralMeshFace const *, GeneralMeshFace *> const & face_map) const
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

      dst.index = index;
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
    intx index;
    bool marked;

}; // class GeneralMeshFace

} // namespace Graphics
} // namespace Thea

#endif
