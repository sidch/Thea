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

#ifndef __Thea_Graphics_GeneralMeshEdge_hpp__
#define __Thea_Graphics_GeneralMeshEdge_hpp__

#include "../Common.hpp"
#include "../AttributedObject.hpp"
#include "../List.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMesh;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshVertex;

template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMeshFace;

/** Edge of GeneralMesh. */
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class /* THEA_API */ GeneralMeshEdge : public AttributedObject<EdgeAttributeT>
{
  public:
    typedef GeneralMesh      <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Mesh;    ///< Parent mesh class.
    typedef GeneralMeshVertex<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;  ///< Vertex of the mesh.
    typedef GeneralMeshFace  <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;    ///< Face of the mesh.

  private:
    typedef List< Face *, AllocatorT<Face *> > FaceCollection;  ///< Collection of faces.

  public:
    typedef typename FaceCollection::iterator        FaceIterator;       ///< Iterator over faces.
    typedef typename FaceCollection::const_iterator  FaceConstIterator;  ///< Const iterator over faces.

    /** Construct from two endpoints. */
    GeneralMeshEdge(Vertex * v0 = nullptr, Vertex * v1 = nullptr) : marked(false), bits(0), internal_bits(0)
    {
      endpoints[0] = v0;
      endpoints[1] = v1;
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a i = 1 the second. */
    Vertex const * getEndpoint(int i) const
    {
      theaAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
      return endpoints[i];
    }

    /** Get an endpoint of the edge. \a i = 0 returns the first endpoint and \a i = 1 the second. */
    Vertex * getEndpoint(int i)
    {
      theaAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
      return endpoints[i];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex const * getOtherEndpoint(Vertex const * endpoint) const
    {
      theaAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
      return endpoints[0] == endpoint ? endpoints[1] : endpoints[0];
    }

    /**
     * Given one endpoint of the edge, get the other one. This function assumes the supplied vertex is indeed a valid endpoint
     * of the edge, and (in release mode) does not check that this is so.
     */
    Vertex * getOtherEndpoint(Vertex const * endpoint)
    {
      theaAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
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
     * Get the index (0 or 1) of the endpoint that this edge shares with another if they are connected, or a negative value if
     * they are not. If the two edges have the same endpoints, the index of an arbitrary one is returned.
     */
    int getCommonEndpoint(GeneralMeshEdge const & other) const
    {
      if (endpoints[0] == other.endpoints[0] || endpoints[0] == other.endpoints[1]) { return 0; }
      if (endpoints[1] == other.endpoints[0] || endpoints[1] == other.endpoints[1]) { return 1; }

      return -1;
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
      theaAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");

      if (numFaces() > 2)  // non-manifold
        return nullptr;

      // Find which incident face has this endpoint as the origin of the edge when stepping round the face. The required edge
      // is then the predecessor of this edge around the face.
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      {
        Face * face = *fi;
        GeneralMeshEdge * prev = face->getPredecessor(this);
        if (prev->hasEndpoint(endpoints[i]))  // found it!
          return prev;
      }

      return nullptr;
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
      theaAssertM(hasEndpoint(endpoint), "GeneralMeshEdge: Vertex is not an endpoint of the edge");
      return nextAroundEndpoint(getEndpointIndex(endpoint));
    }

    /** Get the number of faces incident on the edge. */
    intx numFaces() const { return (intx)faces.size(); }

    /** Get an iterator pointing to the first face incident on the edge. */
    FaceConstIterator facesBegin() const { return faces.begin(); }

    /** Get an iterator pointing to the first face incident on the edge. */
    FaceIterator facesBegin() { return faces.begin(); }

    /** Get an iterator pointing to one position beyond the last face incident on the edge. */
    FaceConstIterator facesEnd() const { return faces.end(); }

    /** Get an iterator pointing to one position beyond the last face incident on the edge. */
    FaceIterator facesEnd() { return faces.end(); }

    /**
     * Check if this is a boundary edge. A boundary edge, by default, is an edge adjacent to <b><i>at most</i></b> one face.
     * Note that this definition considers isolated edges, not adjacent to any faces, to be boundary edges. To exclude such
     * edges, set \a isolated_is_boundary to be false.
     */
    bool isBoundaryEdge(bool isolated_is_boundary = true) const
    { return isolated_is_boundary ? numFaces() <= 1 : numFaces() == 1; }

    /**
     * Check if this is a boundary edge of a patch defined by a subset of mesh faces. Such an edge will be adjacent to
     * <b><i>exactly</i></b> one face of the patch. FaceCollectionT should be a collection with a <tt>find()</tt> member function
     * compatible with <tt>std::set<Face const *>::find()</tt>.
     */
    template <typename FaceCollectionT>
    bool isBoundaryEdge(FaceCollectionT const & patch) const
    {
      intx count = 0;
      for (auto f : faces)
        if (patch.find(f) != patch.end())
          count++;

      return count == 1;
    }

    /** Get the length of the edge (involves a square root since the value is not cached). */
    Real length() const { return (endpoints[0]->getPosition() - endpoints[1]->getPosition()).norm(); }

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
      theaAssertM(i == 0 || i == 1, "GeneralMeshEdge: Invalid endpoint index");
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
                UnorderedMap<Vertex const *, Vertex *> const & vertex_map,
                UnorderedMap<GeneralMeshEdge const *, GeneralMeshEdge *> const & edge_map,
                UnorderedMap<Face const *, Face *> const & face_map) const
    {
      dst.setAttr(this->attr());  // assume attributes can be simply copied

      auto e0 = vertex_map.find(endpoints[0]);
      auto e1 = vertex_map.find(endpoints[1]);
      theaAssertM(e0 != vertex_map.end() && e1 != vertex_map.end(), "GeneralMeshEdge: Edge endpoint not mapped to target mesh");

      dst.endpoints[0] = e0->second;
      dst.endpoints[1] = e1->second;

      dst.faces.clear();
      for (auto fi = faces.begin() ; fi != faces.end(); ++fi)
      {
        auto loc = face_map.find(*fi);
        if (loc != face_map.end()) { dst.faces.push_back(loc->second); }
      }

      dst.marked = marked;
      dst.bits = bits;
      dst.internal_bits = internal_bits;
    }

    Vertex * endpoints[2];
    FaceCollection faces;
    bool marked;
    unsigned char bits;
    mutable unsigned char internal_bits;  // only for use by GeneralMesh

}; // class GeneralMeshEdge

} // namespace Graphics
} // namespace Thea

#endif
