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

#ifndef __Thea_Graphics_DCELFace_hpp__
#define __Thea_Graphics_DCELFace_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "DCELHalfedge.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DCELMesh;

/**
 * Face of DCELMesh.
 *
 * Adapted from: DCELFace class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DCELFace : public NormalAttribute<Vector3>, public AttributedObject<FaceAttribute>
{
  public:
    typedef DCELMesh    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Mesh;      ///< Parent mesh class.
    typedef DCELHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute>  Halfedge;  ///< Halfedge of the mesh.
    typedef Halfedge                                                         Edge;      /**< Typedef for interoperability with
                                                                                             GeneralMesh. */
    typedef DCELVertex  <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Vertex;    ///< Vertex of the mesh.

  private:
    typedef NormalAttribute<Vector3> NormalBaseType;

    /** "Dereference" a halfedge and return the edge itself. */
    struct VertexDeref { Vertex * operator()(Halfedge const * e) const { return const_cast<Vertex *>(e->getOrigin()); } };

    /** Move to the halfedge that links to the next non-null face associated with this vertex. */
    struct VertexIncrement { void operator()(Halfedge const ** e) const { (*e) = (*e)->next(); } };

    /** "Dereference" a halfedge and return the edge itself. */
    struct EdgeDeref { Edge * operator()(Halfedge const * e) const { return const_cast<Edge *>(e); } };

    /** Move to the halfedge that links to the next non-null face associated with this vertex. */
    struct EdgeIncrement { void operator()(Halfedge const ** e) const { (*e) = (*e)->next(); } };

  public:
    /** An iterator over the vertices of the face. */
    typedef DCELInternal::FwdIterator<Halfedge, Vertex, VertexDeref, VertexIncrement> VertexIterator;

    /** A const iterator over the vertices of the face. */
    typedef DCELInternal::FwdIterator<Halfedge const, Vertex const, VertexDeref, VertexIncrement> VertexConstIterator;

    /** An iterator over the edges of the face. */
    typedef DCELInternal::FwdIterator<Halfedge, Edge, EdgeDeref, EdgeIncrement> EdgeIterator;

    /** A const iterator over the edges of the face. */
    typedef DCELInternal::FwdIterator<Halfedge const, Edge const, EdgeDeref, EdgeIncrement> EdgeConstIterator;

    /** Default constructor. */
    DCELFace() : halfedge(NULL), index(-1), num_edges(0) {}

    /** Get a canonical halfedge on the boundary of the face. */
    Halfedge const * getHalfedge() const
    {
      debugAssertM(halfedge, "DCELFace: Face has a null boundary");
      return halfedge;
    }

    /** Get a canonical halfedge on the boundary of the face. */
    Halfedge * getHalfedge()
    {
      debugAssertM(halfedge, "DCELFace: Face has a null boundary");
      return halfedge;
    }

    /** Get an iterator pointing to the first vertex of the face. */
    VertexConstIterator verticesBegin() const { return VertexConstIterator(halfedge); }

    /** Get an iterator pointing to the first vertex of the face. */
    VertexIterator verticesBegin() { return VertexIterator(halfedge); }

    /** Get an iterator pointing to the position beyond the last vertex of the face. */
    VertexConstIterator verticesEnd() const { return VertexConstIterator(halfedge, false); }

    /** Get an iterator pointing to the position beyond the last vertex of the face. */
    VertexIterator verticesEnd() { return VertexIterator(halfedge, false); }

    /** Get an iterator pointing to the first edge of the face. */
    EdgeConstIterator edgesBegin() const { return EdgeConstIterator(halfedge); }

    /** Get an iterator pointing to the first edge of the face. */
    EdgeIterator edgesBegin() { return EdgeIterator(halfedge); }

    /** Get an iterator pointing to the position beyond the last edge of the face. */
    EdgeConstIterator edgesEnd() const { return EdgeConstIterator(halfedge, false); }

    /** Get an iterator pointing to the position beyond the last edge of the face. */
    EdgeIterator edgesEnd() { return EdgeIterator(halfedge, false); }

    /** Get the number of edges bordering this face. */
    int numEdges() const { return num_edges; }

    /** Get the number of vertices of this face (identical to numEdges()). */
    int numVertices() const { return num_edges; }

    /** Check if the face is a triangle. */
    bool isTriangle() const { return num_edges == 3; }

    /** Check if the face is a quad. */
    bool isQuad() const { return num_edges == 4; }

    /** Get the index of the face, typically in the source file (or negative if unindexed). */
    intx getIndex() const { return index; }

    /** Set the index of the face, typically from the source file (or negative if unindexed). */
    void setIndex(intx index_) { index = index_; }

    /** Update the face normal by recomputing it from vertex data. */
    void updateNormal()
    {
      // Assume the normal is consistent across the face
      Vector3 e1 = halfedge->getOrigin()->getPosition()      - halfedge->getEnd()->getPosition();
      Vector3 e2 = halfedge->next()->getEnd()->getPosition() - halfedge->getEnd()->getPosition();
      setNormal(e2.cross(e1).normalized());  // counter-clockwise
    }

    /** Compute the centroid of the face. */
    Vector3 centroid() const
    {
      if (!halfedge)
        return Vector3::Zero();

      Vector3 c = halfedge->getOrigin()->getPosition();
      intx nv = 1;
      Halfedge const * e = halfedge->next();
      while (e != halfedge)
      {
        c += e->getOrigin()->getPosition();
        nv++;
      }

      return c / nv;
    }

    /**
     * Test if the face contains a point (which is assumed to lie on the plane of the face -- for efficiency the function does
     * <b>not</b> explicitly verify that this holds).
     */
    bool contains(Vector3 const & p) const
    {
      if (!halfedge) return false;

      // Generate a ray for the even-odd test, from p to the midpoint of the first halfedge. Ignore degenerate situations for
      // now.
      Vector3 u = 0.5 * (halfedge->getOrigin()->getPosition() + halfedge->getEnd()->getPosition()) - p;

      int count = 1;  // first halfedge is obviously intersected, since we generated the ray through its midpoint
      Halfedge const * e = halfedge->next();
      while (e != halfedge)
      {
        Vector3 v0 = e->getOrigin()->getPosition() - p;
        Vector3 v1 = e->getEnd()->getPosition()    - p;

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

        e = e->next();
      }

      return (count % 2 == 1);
    }

  private:
    template <typename _VertexAttribute, typename _HalfedgeAttribute, typename _FaceAttribute> friend class DCELMesh;

    Halfedge * halfedge;
    intx index;
    int num_edges;

}; // class DCELFace

} // namespace Graphics
} // namespace Thea

#endif
