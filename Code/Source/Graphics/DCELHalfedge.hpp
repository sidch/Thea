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

#ifndef __Thea_Graphics_DCELHalfedge_hpp__
#define __Thea_Graphics_DCELHalfedge_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DCELFace;
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute> class DCELVertex;

/**
 * Halfedge of DCELMesh.
 *
 * Adapted from: DCELHalfedge class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DCELHalfedge : public AttributedObject<HalfedgeAttribute>
{
  public:
    typedef DCELFace  <VertexAttribute, HalfedgeAttribute, FaceAttribute> Face;    ///< Face of the mesh.
    typedef DCELVertex<VertexAttribute, HalfedgeAttribute, FaceAttribute> Vertex;  ///< Vertex of the mesh.

    /** Default constructor. */
    DCELHalfedge(long index_ = -1) : index(index_), twin_he(NULL), next_he(NULL), face(NULL), origin(NULL), bits(0) {}

    /** Get the vertex from which this halfedge originates. */
    Vertex const * getOrigin() const
    {
      debugAssertM(origin, "DCELHalfedge: Halfedge has no origin");
      return origin;
    }

    /** Get the vertex from which this halfedge originates. */
    Vertex * getOrigin()
    {
      debugAssertM(origin, "DCELHalfedge: Halfedge has no origin");
      return origin;
    }

    /** Get the vertex at which this halfedge ends. */
    Vertex const * getEnd() const { return twin()->getOrigin(); }

    /** Get the vertex at which this halfedge ends. */
    Vertex * getEnd() { return twin()->getOrigin(); }

    /** Get the next halfedge around the face. */
    DCELHalfedge const * next() const
    {
      debugAssertM(next_he, "DCELHalfedge: Halfedge has no successor");
      return next_he;
    }

    /** Get the next halfedge around the face. */
    DCELHalfedge * next()
    {
      debugAssertM(next_he, "DCELHalfedge: Halfedge has no successor");
      return next_he;
    }

    /** Get the halfedge between the same two vertices, but in the opposite direction. */
    DCELHalfedge const * twin() const
    {
      debugAssertM(twin_he, "DCELHalfedge: Halfedge has no twin_he");
      return twin_he;
    }

    /** Get the next halfedge around the face. */
    DCELHalfedge * twin()
    {
      debugAssertM(twin_he, "DCELHalfedge: Halfedge has no twin_he");
      return twin_he;
    }

    /**
     * Get the next halfedge around the originating vertex.
     */
    DCELHalfedge const * nextAroundOrigin() const { return twin()->next(); }

    /**
     * Get the next halfedge around the originating vertex.
     */
    DCELHalfedge * nextAroundOrigin() { return twin()->next(); }

    /**
     * Get the previous halfedge around the originating vertex.
     */
    DCELHalfedge const * prevAroundOrigin() const { return const_cast<DCELHalfedge *>(this)->prevAroundOrigin(); }

    /**
     * Get the previous halfedge around the originating vertex.
     */
    DCELHalfedge * prevAroundOrigin()
    {
      DCELHalfedge * rval = this;
      DCELHalfedge * next_around_vertex = twin()->next();

      while (next_around_vertex != this)
      {
        rval = next_around_vertex;
        next_around_vertex = next_around_vertex->twin()->next();
      }

      return rval;
    }

    /** Get the face adjoining this halfedge (or null if this is a border halfedge). */
    Face const * getFace() const { return face; }

    /** Get the face adjoining this halfedge (or null if this is a border halfedge). */
    Face * getFace() { return face; }

    /** Check if this is a boundary halfedge, i.e. if its face pointer is null. */
    bool isBoundary() const { return !face; }

    /** Check if this is a boundary edge, i.e. if either this halfedge or its twin has a null face pointer. */
    bool isBoundaryEdge() const { return !face || !twin()->face; }

    /** Check if one or more marker bits are set. */
    bool areBitsSet(unsigned char mask) const { return ((bits & mask) == mask); };

    /** Set one or more marker bits on or off. */
    void setBits(unsigned char mask, bool value) { bits = (value ? (bits | mask) : (bits & ~mask)); };

    /** Set all marker bits to off. */
    void clearAllBits() { bits = 0; };

  private:
    template <typename _VertexAttribute, typename _HalfedgeAttribute, typename _FaceAttribute> friend class DCELMesh;

    long index;
    DCELHalfedge * twin_he;
    DCELHalfedge * next_he;
    Face * face;
    Vertex * origin;
    unsigned char bits;

#if 0
    /** Check if one or more internal bits are set. */
    bool areInternalBitsSet(unsigned char mask) const { return ((internal_bits & mask) == mask); };

    /** Set one or more internal bits on or off. */
    void setInternalBits(unsigned char mask, bool value) const
    { internal_bits = (value ? (internal_bits | mask) : (internal_bits & ~mask)); };

    /** Set all internal bits to off. */
    void clearAllInternalBits() const { internal_bits = 0; };

    mutable unsigned char internal_bits;  // only for use by DCELMesh
#endif

}; // class DCELHalfedge

} // namespace Graphics
} // namespace Thea

#endif
