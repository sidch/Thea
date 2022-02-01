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

#ifndef __Thea_Graphics_GeneralMeshVertex_hpp__
#define __Thea_Graphics_GeneralMeshVertex_hpp__

#include "../Common.hpp"
#include "../Algorithms/NormalTraitsN.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../AttributedObject.hpp"
#include "../List.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {

namespace Graphics {

// Forward declarations
template <typename VertexAttributeT, typename EdgeAttributeT, typename FaceAttributeT, template <typename T> class AllocatorT>
class GeneralMesh;

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
    typedef GeneralMesh    <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Mesh;  ///< Parent mesh class.
    typedef GeneralMeshEdge<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;  ///< Edge of the mesh.
    typedef GeneralMeshFace<VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;  ///< Face of the mesh.

  private:
    typedef PositionAttribute<Vector3>  PositionBaseType;
    typedef NormalAttribute<Vector3>    NormalBaseType;

    typedef List< Edge *, AllocatorT<Edge *> > EdgeList;
    typedef List< Face *, AllocatorT<Face *> > FaceList;

  public:
    typedef typename EdgeList::iterator        EdgeIterator;       ///< Iterator over edges.
    typedef typename EdgeList::const_iterator  EdgeConstIterator;  ///< Const iterator over edges.
    typedef typename FaceList::iterator        FaceIterator;       ///< Iterator over faces.
    typedef typename FaceList::const_iterator  FaceConstIterator;  ///< Const iterator over faces.

    /** Default constructor. */
    GeneralMeshVertex()
    : NormalBaseType(Vector3::Zero()), index(-1), has_precomputed_normal(false), normal_normalization_factor(0), marked(false)
    {}

    /** Sets the vertex to have a location. */
    explicit GeneralMeshVertex(Vector3 const & p)
    : PositionBaseType(p), NormalBaseType(Vector3::Zero()), index(-1), has_precomputed_normal(false),
      normal_normalization_factor(0), marked(false)
    {}

    /**
     * Sets the vertex to have a location and a precomputed normal. Adding faces will <b>not</b> change the vertex normal if
     * this constructor is used.
     */
    GeneralMeshVertex(Vector3 const & p, Vector3 const & n)
    : PositionBaseType(p), NormalBaseType(n), index(-1), has_precomputed_normal(true), normal_normalization_factor(0),
      marked(false)
    {}

    /**
     * Get the number of edges incident on the vertex. Equivalent to degree().
     *
     * @see degree()
     */
    intx numEdges() const { return (intx)edges.size(); }

    /**
     * Get the degree of the vertex, i.e. number of edges incident on it. Equivalent to numEdges().
     *
     * @see numEdges()
     */
    intx degree() const { return (intx)edges.size(); }

    /** Get the number of faces incident on the vertex. */
    intx numFaces() const { return (intx)faces.size(); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Edge const * getEdgeTo(GeneralMeshVertex const * v) const
    { return const_cast<GeneralMeshVertex *>(this)->getEdgeTo(v); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Edge * getEdgeTo(GeneralMeshVertex const * v)
    {
      if (v == this) return nullptr;

      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
      {
        Edge * e = *ei;
        if (e->hasEndpoint(v)) return e;
      }

      return nullptr;
    }

    /** Check if the vertex is adjacent to a given edge. */
    bool hasEdgeTo(GeneralMeshVertex const * v) const { return getEdgeTo(v) != nullptr; }

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

    /** Get an iterator pointing to the first edge incident on the vertex. */
    EdgeConstIterator edgesBegin() const { return edges.begin(); }

    /** Get an iterator pointing to the first edge incident on the vertex. */
    EdgeIterator edgesBegin() { return edges.begin(); }

    /** Get an iterator pointing to one position beyond the last edge incident on the vertex. */
    EdgeConstIterator edgesEnd() const { return edges.end(); }

    /** Get an iterator pointing to one position beyond the last edge incident on the vertex. */
    EdgeIterator edgesEnd() { return edges.end(); }

    /** Get an iterator pointing to the first face incident on the vertex. */
    FaceConstIterator facesBegin() const { return faces.begin(); }

    /** Get an iterator pointing to the first face incident on the vertex. */
    FaceIterator facesBegin() { return faces.begin(); }

    /** Get an iterator pointing to one position beyond the last face incident on the vertex. */
    FaceConstIterator facesEnd() const { return faces.end(); }

    /** Get an iterator pointing to one position beyond the last face incident on the vertex. */
    FaceIterator facesEnd() { return faces.end(); }

    /** Check if the vertex lies on a mesh boundary. */
    bool isBoundaryVertex() const
    {
      if (edges.empty()) return true;

      for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        if ((*ei)->isBoundaryEdge()) return true;

      return false;
    }

    /** Get the index of the vertex, typically in the source file (or negative if unindexed). */
    intx getIndex() const { return index; }

    /** Set the index of the vertex, typically from the source file (or negative if unindexed). */
    void setIndex(intx index_) { index = index_; }

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
        Vector3 sum_normals = Vector3::Zero();
        for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
          sum_normals += (*fi)->getNormal();  // weight by face area?

        normal_normalization_factor = sum_normals.norm();
        if (normal_normalization_factor < 1e-20f)
          setNormal(Vector3::Zero());
        else
          setNormal(sum_normals / normal_normalization_factor);
      }
      else
      {
        setNormal(Vector3::Zero());
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
        normal_normalization_factor = sum_normals.norm();
        if (normal_normalization_factor < 1e-20f)
          setNormal(Vector3::Zero());
        else
          setNormal(sum_normals / normal_normalization_factor);
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
        normal_normalization_factor = sum_normals.norm();
        if (normal_normalization_factor < 1e-20f)
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

    /** Get the index of the vertex in a GPU array. */
    uint32 getPackingIndex() const { return packing_index; }

    /** Set the index of the vertex in a GPU array. */
    void setPackingIndex(uint32 packing_index_) const { packing_index = packing_index_; }

    /** Make an exact copy of the vertex. */
    void copyTo(GeneralMeshVertex & dst,
                UnorderedMap<GeneralMeshVertex const *, GeneralMeshVertex *> const & vertex_map,
                UnorderedMap<Edge const *, Edge *> const & edge_map,
                UnorderedMap<Face const *, Face *> const & face_map) const
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

      dst.index = index;
      dst.has_precomputed_normal = has_precomputed_normal;
      dst.normal_normalization_factor = normal_normalization_factor;
      dst.packing_index = packing_index;
      dst.marked = marked;
    }

    EdgeList edges;
    FaceList faces;
    intx index;
    bool has_precomputed_normal;
    float normal_normalization_factor;
    mutable uint32 packing_index;
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

// Specify that a mesh vertex has a normal vector. */
template <typename VT, typename ET, typename FT, template <typename T> class AT>
class HasNormalN<Graphics::GeneralMeshVertex<VT, ET, FT, AT>, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its normal vector. */
template <typename VT, typename ET, typename FT, template <typename T> class AT>
class NormalTraitsN<Graphics::GeneralMeshVertex<VT, ET, FT, AT>, 3>
{
  public:
    static Vector3 const & getNormal(Graphics::GeneralMeshVertex<VT, ET, FT, AT> const & t) { return t.getNormal(); }
};

} // namespace Algorithms

} // namespace Thea

#endif
