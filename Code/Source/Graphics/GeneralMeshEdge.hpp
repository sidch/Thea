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

#ifndef __Thea_Graphics_GeneralMeshEdge_hpp__
#define __Thea_Graphics_GeneralMeshEdge_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "../List.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshVertex;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshFace;

/** Edge of GeneralMesh. */
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class /* THEA_API */ GeneralMeshEdge : public AttributedObject<EdgeAttributeT>
{
  public:
    typedef GeneralMeshVertex<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;  ///< Vertex of the mesh.
    typedef GeneralMeshFace  <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;    ///< Face of the mesh.

  private:
    typedef TheaList< Face *, AllocatorT<Face *> > FaceList;

  public:
    typedef typename FaceList::iterator        FaceIterator;       ///< Iterator over faces.
    typedef typename FaceList::const_iterator  FaceConstIterator;  ///< Const iterator over faces.

    /** Construct from two endpoints. */
    GeneralMeshEdge(Vertex * v0 = NULL, Vertex * v1 = NULL) : marked(false), bits(0), internal_bits(0)
    {
      endpoints[0] = v0;
      endpoints[1] = v1;
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a i = 1 the second. */
    Vertex const * getEndpoint(int i) const
    {
      debugAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
      return endpoints[i];
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a i = 1 the second. */
    Vertex * getEndpoint(int i)
    {
      debugAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
      return endpoints[i];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex const * getOtherEndpoint(Vertex const * endpoint) const
    {
      debugAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
      return endpoints[0] == endpoint ? endpoints[1] : endpoints[0];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex * getOtherEndpoint(Vertex const * endpoint)
    {
      debugAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
      return endpoints[0] == endpoint ? endpoints[1] : endpoints[0];
    }

    /** Get the index (0 or 1) of an endpoint given a pointer to it, or a negative value if the neither endpoint matches. */
    int getEndpointIndex(Vertex const * endpoint) const
    {
      return endpoints[0] == endpoint ? 0 : (endpoints[1] == endpoint ? 1 : -1);
    }

    /** Check if the edge has a given vertex as an endpoint. */
    bool hasEndpoint(Vertex const * v) const { return endpoints[0] == v || endpoints[1] == v; }

    /** Check if the edge is adjacent to a given face. */
    bool hasIncidentFace(Face const * face) const
    {
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (*fi == face)
          return true;

      return false;
    }

    /** Check if two edges share the same endpoints. */
    bool isCoincidentTo(GeneralMeshEdge const & other) const
    {
      return (endpoints[0] == other.endpoints[0] && endpoints[1] == other.endpoints[1])
          || (endpoints[0] == other.endpoints[1] && endpoints[1] == other.endpoints[0]);
    }

    /** Check if this edge shares an endpoint with another. */
    bool isConnectedTo(GeneralMeshEdge const & other) const
    {
      return (endpoints[0] == other.endpoints[0] || endpoints[1] == other.endpoints[1]
           || endpoints[0] == other.endpoints[1] || endpoints[1] == other.endpoints[0]);
    }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error, returns null.
     */
    GeneralMeshEdge const * nextAroundEndpoint(int i) const
    { return const_cast<GeneralMeshEdge *>(this)->nextAroundEndpoint(i); }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error, returns null.
     */
    GeneralMeshEdge * nextAroundEndpoint(int i)
    {
      debugAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");

      if (numFaces() > 2)  // non-manifold
        return NULL;

      // Find which incident face has this endpoint as the origin of the edge when stepping round the face. The required edge
      // is then the predecessor of this edge around the face.
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      {
        Face * face = *fi;
        GeneralMeshEdge * prev = face->getPredecessor(this);
        if (prev->hasEndpoint(endpoints[i]))  // found it!
          return prev;
      }

      return NULL;
    }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error, returns null.
     */
    GeneralMeshEdge const * nextAroundEndpoint(Vertex const * endpoint) const
    { return const_cast<GeneralMeshEdge *>(this)->nextAroundEndpoint(endpoint); }

    /**
     * Get the next edge when stepping counter-clockwise (when viewed from the "outside" of the mesh) around a specified
     * endpoint, assuming the neighborhood is manifold. This also assumes that face vertices/edges have consistent
     * counter-clockwise winding order when viewed from the outside. On error, returns null.
     */
    GeneralMeshEdge * nextAroundEndpoint(Vertex const * endpoint)
    {
      debugAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
      return nextAroundEndpoint(getEndpointIndex(endpoint));
    }

    /** Get the number of faces incident on the edge. */
    long numFaces() const { return (long)faces.size(); }

    /** Get an iterator pointing to the first face. */
    FaceConstIterator facesBegin() const { return faces.begin(); }

    /** Get an iterator pointing to the first face. */
    FaceIterator facesBegin() { return faces.begin(); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceConstIterator facesEnd() const { return faces.end(); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceIterator facesEnd() { return faces.end(); }

    /** Check if this is a boundary edge, i.e. if it is adjacent to at most one face. */
    bool isBoundary() const { return numFaces() <= 1; }

    /** Check if one or more marker bits are set. */
    bool areBitsSet(unsigned char mask) const { return ((bits & mask) == mask); };

    /** Set one or more marker bits on or off. */
    void setBits(unsigned char mask, bool value) { bits = (value ? (bits | mask) : (bits & ~mask)); };

    /** Set all marker bits to off. */
    void clearAllBits() { bits = 0; };

  private:
    template <typename V, typename E, typename F, template <typename T> class A> friend class GeneralMesh;

    /** Set an endpoint of the edge. */
    void setEndpoint(int i, Vertex * vertex)
    {
      debugAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
      endpoints[i] = vertex;
    }

    /** Replace all references to a vertex with references to another vertex. */
    void replaceVertex(Vertex * old_vertex, Vertex * new_vertex)
    {
      if (endpoints[0] == old_vertex) endpoints[0] = new_vertex;
      if (endpoints[1] == old_vertex) endpoints[1] = new_vertex;
    }

    /** Add a reference to a face incident at this vertex. */
    void addFace(Face * face) { faces.push_back(face); }

    /** Remove all references to a face incident on this edge. */
    void removeFace(Face * face)
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); )
      {
        if (*fi == face)
          fi = faces.erase(fi);
        else
          ++fi;
      }
    }

    /** Remove a face incident on this edge. */
    FaceIterator removeFace(FaceIterator loc) { return faces.erase(loc); }

    /** Replace all references to a face with references to another face. */
    void replaceFace(Face * old_face, Face * new_face)
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (*fi == old_face)
          *fi = new_face;
    }

    /** Is the edge a self-loop (both endpoints same)? */
    bool isSelfLoop() const { return endpoints[0] == endpoints[1]; }

    /** Mark the edge. */
    void mark()
    {
      marked = true;
    }

    /** Unmark the edge. */
    void unmark()
    {
      marked = false;
    }

    /** Check if the edge is marked. */
    bool isMarked() const { return marked; }

    /** Check if one or more internal bits are set. */
    bool areInternalBitsSet(unsigned char mask) const { return ((internal_bits & mask) == mask); };

    /** Set one or more internal bits on or off. */
    void setInternalBits(unsigned char mask, bool value) const
    { internal_bits = (value ? (internal_bits | mask) : (internal_bits & ~mask)); };

    /** Set all internal bits to off. */
    void clearAllInternalBits() const { internal_bits = 0; };

    /** Make an exact copy of the edge. */
    void copyTo(GeneralMeshEdge & dst,
                TheaUnorderedMap<Vertex const *, Vertex *> const & vertex_map,
                TheaUnorderedMap<GeneralMeshEdge const *, GeneralMeshEdge *> const & edge_map,
                TheaUnorderedMap<Face const *, Face *> const & face_map) const
    {
      dst.setAttr(this->attr());  // assume attributes can be simply copied

      dst.endpoints[0] = vertex_map.find(endpoints[0])->second;  // assume it always exists
      dst.endpoints[1] = vertex_map.find(endpoints[1])->second;  // assume it always exists

      dst.faces.resize(faces.size());
      FaceConstIterator fi = faces.begin();
      FaceIterator dfi = dst.faces.begin();
      for ( ; fi != faces.end(); ++fi, ++dfi)
        *dfi = face_map.find(*fi)->second;  // assume it always exists

      dst.marked = marked;
      dst.bits = bits;
      dst.internal_bits = internal_bits;
    }

    Vertex * endpoints[2];
    FaceList faces;
    bool marked;
    unsigned char bits;
    mutable unsigned char internal_bits;  // only for use by GeneralMesh

}; // class GeneralMeshEdge

} // namespace Graphics
} // namespace Thea

#endif
