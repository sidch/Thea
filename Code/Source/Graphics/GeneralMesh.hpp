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

#ifndef __Thea_Graphics_GeneralMesh_hpp__
#define __Thea_Graphics_GeneralMesh_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Colors.hpp"
#include "../List.hpp"
#include "../MatrixWrapper.hpp"
#include "../NamedObject.hpp"
#include "../Polygon3.hpp"
#include "../UnorderedMap.hpp"
#include "IMesh.hpp"
#include "DefaultMeshCodecs.hpp"
#include "GeneralMeshFace.hpp"
#include "GeneralMeshVertex.hpp"
#include "GeneralMeshEdge.hpp"
#include "GraphicsAttributes.hpp"
#include "IncrementalGeneralMeshBuilder.hpp"
#include "EdgeWelder.hpp"
#include <type_traits>

namespace Thea {
namespace Graphics {

/**
 * A class for storing meshes with arbitrary topologies. Optionally allows GPU-buffered rendering (see important note below).
 *
 * @note When using GPU-buffered rendering, any methods of this class which change the mesh will automatically re-initialize the
 * buffers. <b>However</b>, if <i>external</i> methods change the mesh, such as methods of the GeneralMeshVertex,
 * GeneralMeshEdge and GeneralMeshFace classes, then the user <i>must</i> manually indicate that the mesh needs to be
 * resynchronized with the GPU. The invalidateGpuBuffers() function should be used for this.
 *
 * @todo Add support for GPU-buffered texture coordinates with 1, 3 or 4 dimensions.
 * @todo Instantiate different types of GPU buffers for different types of colors/texture coordinates.
 */
template < typename VertexAttributeT               =  Graphics::NullAttribute,
           typename EdgeAttributeT                 =  Graphics::NullAttribute,
           typename FaceAttributeT                 =  Graphics::NullAttribute,
           template <typename T> class AllocatorT  =  std::allocator >
class /* THEA_API */ GeneralMesh : public NamedObject, public virtual IMesh
{
  public:
    THEA_DECL_SMART_POINTERS(GeneralMesh)

    /** Mesh type tag. */
    struct GENERAL_MESH_TAG {};

    /**< Vertex of the mesh. */
    typedef GeneralMeshVertex <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;

    /**< Edge of the mesh. */
    typedef GeneralMeshEdge   <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;

    /**< Face of the mesh. */
    typedef GeneralMeshFace   <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;

  private:
    typedef List<Vertex, AllocatorT<Vertex> >  VertexList;
    typedef List<Edge,   AllocatorT<Edge>   >  EdgeList;
    typedef List<Face,   AllocatorT<Face>   >  FaceList;

  public:
    typedef typename VertexList::iterator        VertexIterator;       ///< Iterator over vertices.
    typedef typename VertexList::const_iterator  VertexConstIterator;  ///< Const iterator over vertices.
    typedef typename EdgeList::iterator          EdgeIterator;         ///< Iterator over edges.
    typedef typename EdgeList::const_iterator    EdgeConstIterator;    ///< Const iterator over edges.
    typedef typename FaceList::iterator          FaceIterator;         ///< Iterator over faces.
    typedef typename FaceList::const_iterator    FaceConstIterator;    ///< Const iterator over faces.

    // Generic typedefs, each mesh class must define these for builder and codec compatibility
    typedef Vertex        *  VertexHandle;       ///< Handle to a mesh vertex.
    typedef Vertex const  *  VertexConstHandle;  ///< Handle to an immutable mesh vertex.
    typedef Face          *  FaceHandle;         ///< Handle to a mesh face.
    typedef Face   const  *  FaceConstHandle;    ///< Handle to an immutable mesh face.

    /** Identifiers for the various buffers (enum class). */
    struct BufferId
    {
      /** Supported values. */
      enum Value
      {
        ALL              =  0xFFFF,  ///< The set of all GPU buffers.
        VERTEX_POSITION  =  0x0001,  ///< Buffer containing vertex positions.
        VERTEX_NORMAL    =  0x0002,  ///< Buffer containing vertex normals.
        VERTEX_COLOR     =  0x0004,  ///< Buffer containing vertex colors.
        VERTEX_TEXCOORD  =  0x0008,  ///< Buffer containing vertex texture coordinates.
        FACE_NORMAL      =  0x0010,  ///< Buffer containing face normals (currently not used).
        TOPOLOGY         =  0x0020   ///< Buffer(s) containing face indices.
      };

      THEA_ENUM_CLASS_BODY(BufferId)

    }; // struct BufferId

    /** Constructor. */
    GeneralMesh(std::string const & name = "AnonymousMesh")
    : NamedObject(name),
      max_vertex_index(-1),
      max_face_index(-1),
      changed_packed(BufferId::ALL),
      has_large_polys(false),
      buffered_rendering(false),
      buffered_wireframe(false),
      changed_buffers(BufferId::ALL),
      num_tri_indices(0),
      num_quad_indices(0),
      buf_pool(nullptr),
      vertex_positions_buf(nullptr),
      vertex_normals_buf(nullptr),
      vertex_colors_buf(nullptr),
      vertex_texcoords_buf(nullptr),
      tris_buf(nullptr),
      quads_buf(nullptr),
      edges_buf(nullptr),
      vertex_matrix(nullptr, 3, 0),
      tri_matrix(nullptr, 3, 0),
      quad_matrix(nullptr, 4, 0),
      vertex_wrapper(&vertex_matrix),
      tri_wrapper(&tri_matrix),
      quad_wrapper(&quad_matrix)
    {}

    /**
     * Copy constructor. Creates a deep copy of the mesh (including copies of the attributes). <b>Currently not implemented.</b>
     */
    GeneralMesh(GeneralMesh const & src) : NamedObject(src)
    {
      throw Error("GeneralMesh: Copy constructor not currently implemented");
    }

    /**
     * Make an exact copy of the mesh, optionally returning mapping from source to destination vertices/edges/faces. Previous
     * data in the maps is <b>not</b> cleared.
     */
    void copyTo(GeneralMesh & dst,
                UnorderedMap<Vertex const *, Vertex *> * vertex_map = nullptr,
                UnorderedMap<Edge const *, Edge *> * edge_map = nullptr,
                UnorderedMap<Face const *, Face *> * face_map = nullptr) const
    {
      typedef UnorderedMap<Vertex const *, Vertex *> VertexMap;
      typedef UnorderedMap<Edge   const *, Edge   *> EdgeMap;
      typedef UnorderedMap<Face   const *, Face   *> FaceMap;

      VertexMap  tmp_vertex_map;
      EdgeMap    tmp_edge_map;
      FaceMap    tmp_face_map;
      if (!vertex_map) vertex_map  =  &tmp_vertex_map;
      if (!edge_map  ) edge_map    =  &tmp_edge_map;
      if (!face_map  ) face_map    =  &tmp_face_map;

      dst.clear();

      dst.vertices.resize(vertices.size());
      dst.edges.resize(edges.size());
      dst.faces.resize(faces.size());

      // Initialize vertex mapping from source to destination
      {
        VertexConstIterator vi = vertices.begin();
        VertexIterator dvi = dst.vertices.begin();
        for ( ; vi != vertices.end(); ++vi, ++dvi)
          (*vertex_map)[&(*vi)] = &(*dvi);
      }

      // Initialize edge mapping from source to destination
      {
        EdgeConstIterator ei = edges.begin();
        EdgeIterator dei = dst.edges.begin();
        for ( ; ei != edges.end(); ++ei, ++dei)
          (*edge_map)[&(*ei)] = &(*dei);
      }

      // Initialize face mapping from source to destination
      {
        FaceConstIterator fi = faces.begin();
        FaceIterator dfi = dst.faces.begin();
        for ( ; fi != faces.end(); ++fi, ++dfi)
          (*face_map)[&(*fi)] = &(*dfi);
      }

      // Copy vertices
      {
        VertexConstIterator vi = vertices.begin();
        VertexIterator dvi = dst.vertices.begin();
        for ( ; vi != vertices.end(); ++vi, ++dvi)
          vi->copyTo(*dvi, *vertex_map, *edge_map, *face_map);
      }

      // Copy edges
      {
        EdgeConstIterator ei = edges.begin();
        EdgeIterator dei = dst.edges.begin();
        for ( ; ei != edges.end(); ++ei, ++dei)
          ei->copyTo(*dei, *vertex_map, *edge_map, *face_map);
      }

      // Copy faces
      {
        FaceConstIterator fi = faces.begin();
        FaceIterator dfi = dst.faces.begin();
        for ( ; fi != faces.end(); ++fi, ++dfi)
          fi->copyTo(*dfi, *vertex_map, *edge_map, *face_map);
      }

      dst.max_vertex_index = max_vertex_index;
      dst.max_face_index = max_face_index;
      dst.bounds = bounds;
      dst.buffered_rendering = buffered_rendering;
      dst.buffered_wireframe = buffered_wireframe;
      dst.has_large_polys = has_large_polys;
    }

    // Abstract mesh interface
    IDenseMatrix<Real> const * THEA_ICALL getVertexMatrix() const
    {
      // Assume Vector3 is tightly packed and has no padding
      packVertexPositions();
      Vector3 const * buf = (packed_vertex_positions.empty() ? nullptr : &packed_vertex_positions[0]);
      new (&vertex_matrix) VertexMatrix(reinterpret_cast<Real *>(const_cast<Vector3 *>(buf)), 3, numVertices());
      return &vertex_wrapper;
    }

    IDenseMatrix<uint32> const * THEA_ICALL getTriangleMatrix() const
    {
      packTopology();
      uint32 const * buf = (packed_tris.empty() ? nullptr : &packed_tris[0]);
      new (&tri_matrix) TriangleMatrix(const_cast<uint32 *>(buf), 3, num_tri_indices / 3);
      return &tri_wrapper;
    }

    IDenseMatrix<uint32> const * THEA_ICALL getQuadMatrix() const
    {
      packTopology();
      uint32 const * buf = (packed_quads.empty() ? nullptr : &packed_quads[0]);
      new (&quad_matrix) QuadMatrix(const_cast<uint32 *>(buf), 4, num_quad_indices / 4);
      return &quad_wrapper;
    }

    /** Get an iterator pointing to the first vertex. */
    VertexConstIterator verticesBegin() const { return vertices.begin(); }

    /** Get an iterator pointing to the first vertex. */
    VertexIterator verticesBegin() { return vertices.begin(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexConstIterator verticesEnd() const { return vertices.end(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexIterator verticesEnd() { return vertices.end(); }

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

    /** Deletes all data in the mesh and resets automatic element indexing. */
    void clear()
    {
      vertices.clear();
      edges.clear();
      faces.clear();
      bounds = AxisAlignedBox3();
      max_vertex_index = -1;
      max_face_index = -1;

      packed_vertex_positions.clear();
      packed_vertex_normals.clear();
      packed_vertex_colors.clear();
      packed_vertex_texcoords.clear();
      packed_tris.clear();
      packed_quads.clear();
      packed_edges.clear();

      has_large_polys = false;

      invalidateGpuBuffers();
    }

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const { return vertices.empty() && faces.empty() && edges.empty(); }

    /** Get the number of vertices. */
    intx numVertices() const { return (intx)vertices.size(); };

    /** Get the number of edges. */
    intx numEdges() const { return (intx)edges.size(); };

    /** Get the number of faces. */
    intx numFaces() const { return (intx)faces.size(); };

    /** Compute the number of triangles in the mesh. */
    intx numTriangles() const
    {
      intx rval = 0;
      for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isTriangle())
          rval++;

      return rval;
    }

    /** Compute the number of quads in the mesh. */
    intx numQuads() const
    {
      intx rval = 0;
      for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isQuad())
          rval++;

      return rval;
    }

    /** Recompute and cache the bounding box for the mesh. Make sure this has been called before calling getBounds(). */
    void updateBounds()
    {
      bounds = AxisAlignedBox3();
      for (auto vi = verticesBegin(); vi != verticesEnd(); ++vi)
        bounds.merge(vi->getPosition());
    }

    /**
     * Get the cached bounding box of the mesh. Will be out-of-date unless updateBounds() has been called after all
     * modifications.
     */
    AxisAlignedBox3 const & getBounds() const { return bounds; }

    /** Do the mesh vertices have attached colors? */
    bool hasVertexColors() const { return HasColor<Vertex>::value; }

    /** Do the mesh vertices have attached texture coordinates? */
    bool hasVertexTexCoords() const { return HasTexCoord<Vertex>::value; }

    /**
     * Add a vertex to the mesh, with optional precomputed normal, color, texture coordinates and index. If the index is
     * negative, a new, unique index is generated for the vertex. Automatically calls invalidateGpuBuffers() to schedule a
     * resync with the GPU.
     *
     * @return A pointer to the newly created vertex on success, null on failure.
     */
    Vertex * addVertex(Vector3 const & point, intx index = -1, Vector3 const * normal = nullptr,
                       ColorRgba const * color = nullptr, Vector2 const * texcoord = nullptr)
    {
      if (normal)
        vertices.push_back(Vertex(point, *normal));
      else
        vertices.push_back(Vertex(point));

      Vertex * vertex = &(*vertices.rbegin());
      if (color)     setVertexColor<Vertex>(vertex, *color);
      if (texcoord)  setVertexTexCoord<Vertex>(vertex, *texcoord);

      if (index < 0)
        index = (++max_vertex_index);
      else if (index > max_vertex_index)
        max_vertex_index = index;

      vertex->setIndex(index);

      invalidateGpuBuffers();
      return vertex;
    }

    /**
     * Add a face to the mesh, specified by the sequence of vertices obtained by dereferencing [vbegin, vend) and an optional
     * index. VertexInputIterator must dereference to a pointer to a Vertex. If the index is negative, a new, unique index is
     * generated for the face. Unless the mesh is already in an inconsistent state, failure to add the face will not affect the
     * mesh.
     *
     * Automatically calls invalidateGpuBuffers() to schedule a resync with the GPU.
     *
     * @return A pointer to the newly created face, or null on error.
     */
    template <typename VertexInputIterator>
    Face * addFace(VertexInputIterator vbegin, VertexInputIterator vend, intx index = -1)
    {
      // Create the (initially empty) face
      faces.push_back(Face());
      Face * face = &(*faces.rbegin());

      // Initialize the face
      face = initFace(face, vbegin, vend);  // invalidates GPU buffers
      if (face)
      {
        if (index < 0)
          index = (++max_face_index);
        else if (index > max_face_index)
          max_face_index = index;

        face->setIndex(index);
      }
      else
        faces.pop_back();

      return face;
    }

    /**
     * Remove a face of the mesh. This does NOT remove any vertices or edges. Use removeIsolatedEdges() and
     * removeIsolatedVertices() after calling this function one or more times. Iterators to the face list remain valid unless
     * the iterator pointed to the removed face.
     *
     * This is a relatively slow operation since the face needs to be looked up in the face list
     * (linear in number of faces). For speed, use removeFace(FaceIterator).
     *
     * @return True if the face was found and removed, else false.
     */
    bool removeFace(Face * face)
    {
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if (&(*fi) == face)
          return removeFace(fi);  // invalidates GPU buffers

      return false;
    }

    /**
     * Remove a face of the mesh. This does NOT remove any vertices or edges. Use removeIsolatedEdges() and
     * removeIsolatedVertices() after calling this function one or more times. Iterators to the face list remain valid unless
     * the iterator pointed to the removed face.
     *
     * Use this version in preference to removeFace(Face const *) where possible.
     *
     * @return True if the face was found and removed, else false.
     */
    bool removeFace(FaceIterator face)
    {
      Face * fp = &(*face);
      unlinkFace(fp);
      faces.erase(face);

      invalidateGpuBuffers();
      return true;
    }

    /**
     * Check if the mesh represents a topological manifold.
     *
     * @param require_closed If true, the mesh cannot have boundary edges.
     * @param require_connected If true, the mesh must have a single connected component.
     *
     * @return True if the mesh is manifold, else false.
     */
    bool isManifold(bool require_closed = true, bool require_connected = false) const
    {
      alwaysAssertM(!require_connected, "GeneralMesh: Testing for connected manifolds not implemented");

      // Check that each edge has 2 or (if open surfaces are ok) 1 incident face(s).
      for (auto ei = edges.begin(); ei != edges.end(); ++ei)
      {
        switch (ei->numFaces())
        {
          case 1:
            if (require_closed) return false;
            break;

          case 2: break;
          default: return false;
        }
      }

      // Now check that each vertex has a manifold neighborhood
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
        if (!isManifoldVertex(&(*vi)))
          return false;

      return true;
    }

    /**
     * Weld boundary edges that are approximately coincident.
     *
     * @todo Test this properly.
     */
    void sealSeams(Real weld_radius)
    {
      EdgeWelder welder(weld_radius);

      for (auto ei = edges.begin(); ei != edges.end(); ++ei)
      {
        if (!ei->isBoundaryEdge())
          continue;

        if (ei->isSelfLoop())
          continue;

        Edge * edge = &(*ei);
        Edge * existing = (Edge *)welder.getUndirectedEdge(edge->getEndpoint(0)->getPosition(),
                                                           edge->getEndpoint(1)->getPosition());
        bool can_seal = (existing != nullptr);
        if (can_seal)
        {
          // We can seal only if the two edges don't share a face
          for (auto fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
            if (existing->hasIncidentFace(*fi))
            {
              can_seal = false;
              break;
            }
        }

        if (can_seal)
          replaceEdge(edge, existing);
        else
          welder.addEdge(edge, edge->getEndpoint(0)->getPosition(), edge->getEndpoint(1)->getPosition());
      }

      removeIsolatedEdges();
      invalidateGpuBuffers();
    }

    /**
     * Replace one vertex with another. The old vertex becomes isolated after this operation, but is <b>not</b> automatically
     * removed from the mesh. You might want to update the normal of \a new_vertex after this operation.
     */
    void replaceVertex(Vertex * old_vertex, Vertex * new_vertex)
    {
      if (new_vertex == old_vertex)  // nothing to do
        return;

      // Check if the vertices are connected by an edge
      Edge * connecting_edge = new_vertex->getEdgeTo(old_vertex);
      if (connecting_edge)
      {
        collapseEdge(connecting_edge, connecting_edge->getEndpointIndex(new_vertex));
      }
      else
      {
        for (auto ei = old_vertex->edgesBegin(); ei != old_vertex->edgesEnd(); ++ei)
        {
          (*ei)->replaceVertex(old_vertex, new_vertex);

          if (!new_vertex->hasIncidentEdge(*ei))
            new_vertex->addEdge(*ei);
        }

        for (auto fi = old_vertex->facesBegin(); fi != old_vertex->facesEnd(); ++fi)
        {
          (*fi)->replaceVertex(old_vertex, new_vertex);

          if (!new_vertex->hasIncidentFace(*fi))
            new_vertex->addFace(*fi);
        }

        old_vertex->edges.clear();
        old_vertex->faces.clear();
      }

      invalidateGpuBuffers();
    }

    /**
     * Collapse an edge, retaining just one of its endpoints. The edge and the other endpoint become isolated and may be removed
     * with removeIsolatedEdges() and removeIsolatedVertices().
     *
     * @param edge The edge to collapse.
     * @param endpoint_to_preserve The index (0 or 1) of the endpoint of the edge that will replace the latter upon collapse.
     */
    void collapseEdge(Edge * edge, int endpoint_to_preserve = 0)
    {
      if (!edge)
        return;

      if (edge->isSelfLoop())
      {
        THEA_WARNING << getName() << ": Edge is self-loop";
        return;
      }

      debugAssertM(endpoint_to_preserve == 0 || endpoint_to_preserve == 1,
                   getNameStr() + ": Endpoint to preserve during edge collapse must be indexed by 0 or 1");

      Vertex * vertex_to_preserve  =  edge->getEndpoint(endpoint_to_preserve);
      Vertex * vertex_to_remove    =  edge->getEndpoint(1 - endpoint_to_preserve);

      // Remove all references to this edge, and the vertex to remove, from all adjacent faces
      bool preserve_edge = false;
      for (auto fi = edge->facesBegin(); fi != edge->facesEnd(); )
      {
        Face * face = *fi;
        auto fei = face->edgesBegin();
        auto fvi = face->verticesBegin();

        bool preserve_edge_ref_to_face = false;
        while (fei != face->edgesEnd())
        {
          if (*fei == edge)
          {
            if (*fvi == vertex_to_remove)
            {
              fei = face->removeEdge(fei);
              fvi = face->removeVertex(fvi);
            }
            else
            {
              // We expect the next vertex to be the one to remove. We'll have to advance the vertex iterator, delete the
              // vertex, and then advance the edge iterator as well
              auto next = fvi; ++next;
              if (next == face->verticesEnd())
                next = face->verticesBegin();

              if (*next == vertex_to_remove)
              {
                fei = face->removeEdge(fei);
                fvi = face->removeVertex(next);

                // Now we must skip over the next edge as well, else edge and vertex iterators will be out of sync
                if (fei != face->edgesEnd())
                {
                  if (*fei == edge)  // somehow this has happened
                  {
                    THEA_WARNING << getName() << ": Face has repeated edge";

                    preserve_edge_ref_to_face = true;
                  }

                  ++fei;
                }
              }
              else
              {
                THEA_WARNING << getName() << ": Edge does not correspond to pair of consecutive vertices";

                preserve_edge_ref_to_face = true;

                ++fei;
                ++fvi;
              }
            }
          }
          else
          {
            ++fei;
            ++fvi;
          }
        }

        if (!preserve_edge_ref_to_face)
          fi = edge->removeFace(fi);
        else
        {
          preserve_edge = true;
          ++fi;
        }
      }

      // Update references of elements adjacent to the vertex that's going to be removed
      for (auto ei = vertex_to_remove->edgesBegin(); ei != vertex_to_remove->edgesEnd(); ++ei)
      {
        // We do the replacement even for the edge that is getting collapsed, since if a face repeats this edge it is going to
        // be present in the output (as a self-loop)

        (*ei)->replaceVertex(vertex_to_remove, vertex_to_preserve);

        if (!(*ei)->hasEndpoint(vertex_to_preserve))
          vertex_to_preserve->addEdge(*ei);
      }

      for (auto fi = vertex_to_remove->facesBegin(); fi != vertex_to_remove->facesEnd(); ++fi)
      {
        (*fi)->replaceVertex(vertex_to_remove, vertex_to_preserve);

        if (!vertex_to_preserve->hasIncidentFace(*fi))
          vertex_to_preserve->addFace(*fi);
      }

      // Remove references to this edge from the vertex that will be preserved (the other one will have its references deleted
      // wholesale later)
      if (!preserve_edge)
        vertex_to_preserve->removeEdge(edge);

      vertex_to_remove->edges.clear();
      vertex_to_remove->faces.clear();

      invalidateGpuBuffers();
    }

    /**
     * Split an edge by routing it through an existing vertex. The edge is updated by setting the existing vertex as its second
     * endpoint, and a new edge is created from the existing vertex to the initial second endpoint.
     *
     * @return The new edge created by the operation, from the existing vertex to the initial second endpoint of the existing
     *   edge.
     */
    Edge * splitEdge(Edge * edge, Vertex * vertex)
    {
      if (!edge)
        return nullptr;

      if (edge->isSelfLoop())
      {
        THEA_DEBUG << getName() << ": Can't split self-loop edge";
        return nullptr;
      }

      if (edge->hasEndpoint(vertex))
      {
        THEA_DEBUG << getName() << ": Can't split edge at existing endpoint";
        return nullptr;
      }

      for (auto efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        if ((*efi)->hasVertex(vertex))
        {
          THEA_DEBUG << getName() << ": Can't split edge at vertex on the same face";
          return nullptr;
        }

      Vertex * old_e1 = edge->getEndpoint(1);
      edge->setEndpoint(1, vertex);
      edges.push_back(Edge(vertex, old_e1));
      Edge * new_edge = &edges.back();

      vertex->addEdge(edge);
      vertex->addEdge(new_edge);

      old_e1->removeEdge(edge);
      old_e1->addEdge(new_edge);

      for (auto fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
      {
        new_edge->addFace(*fi);

        if (!vertex->hasIncidentFace(*fi))
          vertex->addFace(*fi);
      }

      // Now everything's ok except except the references in the faces
      for (auto fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
      {
        Face * face = *fi;

        auto ei = face->edgesBegin();
        auto vi = face->verticesBegin();

        while (ei != face->edgesEnd())
        {
          if (*ei == edge)
          {
            // Insert before or after?
            if (*vi == old_e1)  // sequence is [v1, v0], insert before
            {
              face->edges.insert(ei, new_edge);
            }
            else  // sequence is [v0, v1], insert after
            {
              auto next_ei = ei; ++next_ei;
              face->edges.insert(next_ei, new_edge);
            }

            // Vertex goes in the middle, i.e. always after the current vertex
            auto next_vi = vi; ++next_vi;
            face->vertices.insert(next_vi, vertex);

            break;
          }
          else
          {
            ++ei;
            ++vi;
          }
        }
      }

      invalidateGpuBuffers();
      return new_edge;
    }

    /** Split an edge by creating a new vertex along its length at a given position. */
    Vertex * splitEdge(Edge * edge, Vector3 const & p)
    {
      Real s = (p - edge->getEndpoint(1)->getPosition()).norm();
      Real t = (p - edge->getEndpoint(0)->getPosition()).norm();
      Real sum = s + t;
      s /= sum;
      t /= sum;
      Vector3 n = s * edge->getEndpoint(0)->getNormal() + t * edge->getEndpoint(1)->getNormal();
      Vertex * new_vx = addVertex(p, &n);

      if (!splitEdge(edge, new_vx))  // should generally never happen
      {
        vertices.pop_back();  // remove the vertex we just added
        return nullptr;
      }

      invalidateGpuBuffers();
      return new_vx;
    }

    /** Remove all isolated vertices and edges, and empty faces. */
    void removeDanglers()
    {
      // This sequence removes isolated edges and vertices twice, which is currently necessary
      removeDegenerateEdges();
      removeDegenerateFaces();
    }

    /** Remove degenerate edges (isolated or self-loops) from the mesh. */
    void removeDegenerateEdges()
    {
      // Remove self-loops
      for (auto ei = edgesBegin(); ei != edgesEnd(); ++ei)
      {
        if (!ei->isSelfLoop())
          continue;

        Edge * edge = &(*ei);
        for (auto efi = edge->facesBegin(); efi != edge->facesEnd(); )
        {
          Face * face = *efi;

          auto fei = face->edgesBegin();
          auto fvi = face->verticesBegin();
          while (fei != face->edgesEnd())
          {
            if (*fei == edge)
            {
              alwaysAssertM(*fvi == edge->getEndpoint(0), std::string(getName()) + ": Edge and vertex sequences out of sync");

              fei = face->removeEdge(fei);
              fvi = face->removeVertex(fvi);
            }
            else
            {
              ++fei;
              ++fvi;
            }
          }

          efi = edge->removeFace(efi);
        }
      }

      // Remove isolated edges, including those that became isolated above
      removeIsolatedEdges();
      invalidateGpuBuffers();
    }

    /** Remove all empty faces (fewer than 3 edges) from the mesh. */
    void removeDegenerateFaces()
    {
      for (auto fi = faces.begin(); fi != facesEnd(); )
      {
        if (fi->numVertices() < 3)
        {
          Face * face = &(*fi);

          for (auto vi = face->verticesBegin(); vi != face->verticesEnd(); ++vi)
            (*vi)->removeFace(face);

          for (auto ei = face->edgesBegin(); ei != face->edgesEnd(); ++ei)
            (*ei)->removeFace(face);

          fi = faces.erase(fi);
        }
        else
          ++fi;
      }

      removeIsolatedEdges();  // also removes isolated vertices
      invalidateGpuBuffers();
    }

    /** Remove all isolated edges (no incident faces) from the mesh. */
    void removeIsolatedEdges()
    {
      for (auto ei = edgesBegin(); ei != edgesEnd(); )
      {
        if (ei->faces.empty())
        {
          ei->getEndpoint(0)->removeEdge(&(*ei));
          ei->getEndpoint(1)->removeEdge(&(*ei));
          ei = edges.erase(ei);
        }
        else
          ++ei;
      }

      removeIsolatedVertices();  // some vertices might have become isolated because of edge removal
      invalidateGpuBuffers(BufferId::TOPOLOGY);
    }

    /** Remove all isolated vertices (no incident edges or faces) from the mesh. */
    void removeIsolatedVertices()
    {
      for (auto vi = verticesBegin(); vi != verticesEnd(); )
      {
        if (vi->faces.empty() && vi->edges.empty())
          vi = vertices.erase(vi);
        else
          ++vi;
      }

      invalidateGpuBuffers();
    }

    /** Remove all vertices marked for deletion by other operations. */
    void removeMarkedVertices()
    {
      for (auto vi = verticesBegin(); vi != verticesEnd(); )
      {
        if (vi->isMarked())
          vi = vertices.erase(vi);
        else
          ++vi;
      }

      invalidateGpuBuffers();
    }

    /** Remove all edges marked for deletion by other operations. */
    void removeMarkedEdges()
    {
      for (auto ei = edgesBegin(); ei != edgesEnd(); )
      {
        if (ei->isMarked())
          ei = edges.erase(ei);
        else
          ++ei;
      }

      invalidateGpuBuffers(BufferId::TOPOLOGY);
    }

    /** Remove all faces marked for deletion by other operations. */
    void removeMarkedFaces()
    {
      for (auto fi = facesBegin(); fi != facesEnd(); )
      {
        if (fi->isMarked())
          fi = faces.erase(fi);
        else
          ++fi;
      }

      invalidateGpuBuffers(BufferId::TOPOLOGY);
    }

    /**
     * Triangulate all faces with more than 3 vertices.
     *
     * @param epsilon A tolerance threshold to decide if a triangle is degenerate or not. A negative value selects a default
     *   setting.
     *
     * @return The number of triangulated faces, or a negative number on error. (The number of generated triangles can be
     *   obtained by comparing the number of mesh faces before and after the operation.)
     */
    intx triangulate(Real epsilon = -1)
    {
      intx orig_num_faces = numFaces();
      intx num_visited_faces = 0;
      intx num_triangulated_faces = 0;

      // New faces will be added to the end of the face list, so we can just keep track of when we've processed the original
      // number of faces
      for (auto fi = facesBegin(); num_visited_faces < orig_num_faces; ++fi, ++num_visited_faces)
        if (fi->numVertices() > 3)
        {
          intx nt = triangulate(&(*fi), epsilon);
          if (nt < 0)
            return nt;

          num_triangulated_faces++;
        }

      return num_triangulated_faces;
    }

    /**
     * Triangulate a face if it has more than 3 vertices.
     *
     * @param face The face to triangulate.
     * @param epsilon A tolerance threshold to decide if a triangle is degenerate or not. A negative value selects a default
     *   setting.
     *
     * @return The number of triangles resulting from the operation (1 if the face is already a triangle, negative on error).
     */
    intx triangulate(Face * face, Real epsilon = -1)
    {
      if (!face)
        return 0;

      if (face->numVertices() <= 3)
        return 1;

      intx ntris = 0;
      if (face->numVertices() == 4)
      {
        Vertex * face_vertices[4];
        {
          size_t i = 0;
          for (auto fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi, ++i)
            face_vertices[i] = *fvi;
        }

        intx tri_indices[6];
        ntris = Polygon3::triangulateQuad(face_vertices[0]->getPosition(),
                                          face_vertices[1]->getPosition(),
                                          face_vertices[2]->getPosition(),
                                          face_vertices[3]->getPosition(),
                                          tri_indices[0], tri_indices[1], tri_indices[2],
                                          tri_indices[3], tri_indices[4], tri_indices[5], epsilon);
        if (ntris >= 1)
        {
          if (!replaceFaceWithTriangulation(face, &face_vertices[0], ntris, &tri_indices[0]))
            return -1;
        }
      }
      else
      {
        Array<Vertex *> face_vertices((size_t)face->numVertices());
        Polygon3 poly;
        {
          size_t i = 0;
          for (auto fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi, ++i)
          {
            face_vertices[i] = *fvi;
            poly.addVertex((*fvi)->getPosition());
          }
        }

        Array<intx> tri_indices;
        ntris = poly.triangulate(tri_indices, epsilon);
        if (ntris >= 1)
        {
          if (!replaceFaceWithTriangulation(face, &face_vertices[0], ntris, &tri_indices[0]))  // invalidates GPU buffers
            return -1;
        }
      }

      // If the triangulation generated no triangles, the face is degenerate but we must replace it with SOME triangle
      if (ntris <= 0)
      {
        Vertex * face_vertices[3];
        auto fvi = face->verticesBegin();
        face_vertices[0] = *(fvi++);
        face_vertices[1] = *(fvi++);
        face_vertices[2] = *(fvi);

        intx tri_indices[3] = { 0, 1, 2 };
        ntris = 1;

        if (!replaceFaceWithTriangulation(face, &face_vertices[0], ntris, &tri_indices[0]))
          return -1;
      }

      return ntris;
    }

    /**
     * Check if GPU-buffered rendering is on or off. If it is on, you <b>must manually call</b> invalidateGpuBuffers() every
     * time the mesh changes, to make sure the GPU buffers are update when the mesh is next rendered.
     *
     * @see setGpuBufferedRendering()
     */
    bool renderingIsGpuBuffered() const { return buffered_rendering; }

    /**
     * Turn GPU-buffered rendering on/off. If you enable this function, you <b>must manually call</b> invalidateGpuBuffers()
     * every time the mesh changes, to make sure the GPU buffers are synchronized when the mesh is next rendered.
     *
     * @see renderingIsGpuBuffered()
     */
    void setGpuBufferedRendering(bool value)
    {
      if (value == buffered_rendering)
        return;

      buffered_rendering = value;
      invalidateGpuBuffers();
    }

    /** Invalidate part or all of the current GPU data for the mesh. */
    void invalidateGpuBuffers(int changed_buffers_ = BufferId::ALL)
    {
      changed_buffers |= changed_buffers_;
      changed_packed  |= changed_buffers_;
    }

    /**
     * Enable/disable drawing the edges of the mesh in GPU-buffered mode. Enabling this function will <b>not</b> draw any edges
     * unless you turn on the appropriate RenderOptions flag. The edges will be uploaded to the graphics system on the next call
     * to uploadToGraphicsSystem().
     *
     * Wireframe mode is initially disabled to save video memory. This flag is ignored in non-buffered (immediate) mode.
     *
     * @see wireframeIsGpuBuffered()
     */
    void setGpuBufferedWireframe(bool value)
    {
      if (value == buffered_wireframe)
        return;

      buffered_wireframe = value;
      invalidateGpuBuffers();
    }

    /**
     * Check if wireframe drawing is enabled in GPU-buffered mode. This is initially disabled to save video memory. This flag is
     * ignored in non-buffered (immediate) mode.
     *
     * @see setGpuBufferedWireframe()
     */
    bool wireframeIsGpuBuffered() const { return buffered_wireframe; }

    /**
     * Pack all mesh data to arrays that can be directly transferred to the GPU. Packed arrays that are already synchronized
     * with the mesh are not re-packed.
     */
    void packArrays() const
    {
      packVertexPositions();
      packVertexNormals();
      packVertexColors<Vertex>();
      packVertexTexCoords<Vertex>();
      packTopology();
    }

    int8 THEA_ICALL draw(IRenderSystem * render_system, IRenderOptions const * options = nullptr) const
    {
      if (!render_system) { THEA_ERROR << getName() << ": Can't display mesh on a null rendersystem"; return false; }
      if (!options) options = RenderOptions::defaults();

      if (buffered_rendering)
        return drawBuffered(*render_system, *options);
      else
        return drawImmediate(*render_system, *options);
    }

  private:
    /**
     * Initialize a pre-constructed face, which will be assigned the sequence of vertices obtained by dereferencing
     * [vbegin, vend). VertexInputIterator must dereference to a pointer to a Vertex. Unless the mesh is already in an
     * inconsistent state, failure to add the face will not affect the mesh.
     *
     * Automatically calls invalidateGpuBuffers() to schedule a resync with the GPU.
     *
     * @param face The previously constructed face, assumed to be <b>uninitialized</b>.
     * @param vbegin Points to the beginning of the vertex sequence.
     * @param vend Points to (one past) the end of the vertex sequence.
     *
     * @return A pointer to the face, or null on error.
     */
    template <typename VertexInputIterator>
    Face * initFace(Face * face, VertexInputIterator vbegin, VertexInputIterator vend)
    {
      debugAssertM(face, getNameStr() + ": Null face cannot be initialized");

      // Check for errors and compute normal
      size_t num_verts = 0;
      Vector3 v[3];
      for (auto vi = vbegin; vi != vend; ++vi, ++num_verts)
      {
        debugAssertM(*vi, getNameStr() + ": Null vertex pointer specified for new face");
        if (num_verts < 3) v[num_verts] = (*vi)->getPosition();
      }

      if (num_verts < 3)
      {
        THEA_WARNING << getName() << ": Skipping face -- too few vertices (" << num_verts << ')';
        return nullptr;
      }

      face->clear();

      // Add the loop of vertices to the face
      VertexInputIterator next = vbegin;
      for (auto vi = next++ ; vi != vend; ++vi, ++next)
      {
        if (next == vend) next = vbegin;

        face->addVertex(*vi);
        (*vi)->addFace(face, false);  // we'll update the normals later

        Edge * edge = (*vi)->getEdgeTo(*next);
        if (!edge)
        {
          edges.push_back(Edge(*vi, *next));
          edge = &(*edges.rbegin());

          (*vi)->addEdge(edge);
          (*next)->addEdge(edge);
        }

        edge->addFace(face);
        face->addEdge(edge);
      }

      // Update the face and vertex normals;
      face->updateNormal();
      for (auto fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi)
        (*fvi)->addFaceNormal(face->getNormal());  // weight by face area?

      invalidateGpuBuffers();
      return face;
    }

    /** Set vertex color. */
    template < typename VertexT, typename std::enable_if< HasColor<VertexT>::value, int >::type = 0 >
    void setVertexColor(VertexT * vertex, ColorRgba const & color)
    {
      vertex->attr().setColor(color);
      invalidateGpuBuffers(BufferId::VERTEX_COLOR);
    }

    /** Set vertex color (no-op, called if vertex does not have color attribute). */
    template < typename VertexT, typename std::enable_if< !HasColor<VertexT>::value, int >::type = 0 >
    void setVertexColor(VertexT * vertex, ColorRgba const & color)
    {}

    /** Set vertex texture coordinates. */
    template < typename VertexT, typename std::enable_if< HasTexCoord<VertexT>::value, int >::type = 0 >
    void setVertexTexCoord(VertexT * vertex, Vector2 const & texcoord)
    {
      vertex->attr().setTexCoord(texcoord);
      invalidateGpuBuffers(BufferId::VERTEX_TEXCOORD);
    }

    /** Set vertex texture coordinates (no-op, called if vertex does not have texture coordinate attribute). */
    template < typename VertexT, typename std::enable_if< !HasTexCoord<VertexT>::value, int >::type = 0 >
    void setVertexTexCoord(VertexT * vertex, Vector2 const & texcoord)
    {}

    /**
     * "Step around" the vertex. Given a vertex V, an incident face F (which may be null), and an edge E incident on both V and
     * F, the function finds another face F' incident on V and E, and another edge E' incident on V and F'.
     *
     * @param vertex The vertex to step around.
     * @param face The face F (may be null).
     * @param edge The edge E.
     * @param next_face Used to return the face F'.
     * @param next_edge Used to return the edge E'.
     *
     * @return True if a successor was successfully found, else false.
     */
    bool getSuccessor(Vertex const * vertex, Face const * face, Edge const * edge,
                      Face const ** next_face, Edge const ** next_edge) const
    {
      *next_edge = nullptr;
      *next_face = nullptr;

      // Find the next face around the vertex
      for (auto efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        if (*efi != face)
        {
          *next_face = *efi;
          break;
        }

      // Have a we hit a boundary?
      if (!(*next_face))
        return false;

      // Find the next edge around the vertex
      for (auto fei = (*next_face)->edgesBegin(); fei != (*next_face)->edgesEnd(); ++fei)
        if (*fei != edge && (*fei)->hasEndpoint(vertex))
        {
          *next_edge = *fei;
          break;
        }

      return (bool)(*next_edge);
    }

    /** Check that the region around a vertex represents a manifold surface. */
    bool isManifoldVertex(Vertex const * vertex) const
    {
      // Check if the vertex is isolated
      if (vertex->numEdges() < 2)
        return false;

      // The incident edges must form a continuous sequence around the vertex

      // Find an edge to start from. If the vertex is a boundary vertex, we *must* start from a boundary edge
      Edge const * first_edge = *vertex->edgesBegin();
      for (auto vei = vertex->edgesBegin(); vei != vertex->edgesEnd(); ++vei)
      {
        (*vei)->clearAllInternalBits();

        if ((*vei)->isBoundaryEdge())
          first_edge = *vei;
      }

      // Check if the edge is isolated
      if (first_edge->numFaces() <= 0)
        return false;

      // Start from the edge found above and an incident face, and step around the vertex
      Edge const * edge = first_edge;
      Face const * face = edge->isBoundaryEdge() ? nullptr : *edge->facesBegin();  // if starting from a boundary edge, initial
                                                                                   // "face" is empty space before first edge
      intx num_visited_edges = 0;
      while (true)
      {
        if (++num_visited_edges >= vertex->numEdges())
          break;

        // Step around the vertex
        Edge const * next_edge = nullptr;
        Face const * next_face = nullptr;
        if (!getSuccessor(vertex, face, edge, &next_face, &next_edge))
          break;

        // If we've already visited the successor edge, we've looped round too early
        if (next_edge->areInternalBitsSet(0xFF))
          return false;

        // Mark the current edge as visited
        edge->setInternalBits(0xFF, true);

        // Update the pointers
        edge = next_edge;
        face = next_face;
      }

      // All incident edges should have been visited at this stage
      return (num_visited_edges >= vertex->numEdges());
    }

    /**
     * Replace all references to \a old_edge with references to \a new_edge.
     *
     * @return True if the operation succeeded, else false.
     */
    bool replaceEdge(Edge * old_edge, Edge * new_edge)
    {
      if (old_edge == new_edge)
        return false;

      Vertex * old_v[2] = { old_edge->getEndpoint(0), old_edge->getEndpoint(1) };
      Vertex * new_v[2] = { new_edge->getEndpoint(0), new_edge->getEndpoint(1) };

      if (old_v[0] == old_v[1])  // can't handle degenerate cases
      {
        THEA_DEBUG << getName() << ": Can't replace a self-loop";
        return false;
      }

      if (new_v[0] == new_v[1])  // can't handle degenerate cases
      {
        THEA_DEBUG << getName() << ": Can't replace with a self-loop";
        return false;
      }

      Edge const * e00 = old_v[0]->getEdgeTo(new_v[0]);
      Edge const * e01 = old_v[0]->getEdgeTo(new_v[1]);
      Edge const * e10 = old_v[1]->getEdgeTo(new_v[0]);
      Edge const * e11 = old_v[1]->getEdgeTo(new_v[1]);
      if ((e00 && e00 != old_edge && e00 != new_edge)
       || (e01 && e01 != old_edge && e01 != new_edge)
       || (e10 && e10 != old_edge && e10 != new_edge)
       || (e11 && e11 != old_edge && e11 != new_edge))
      {
        THEA_DEBUG << getName() << ": Can't replace an edge with another to which it has a connecting edge";
        return false;
      }

      for (auto efi = old_edge->facesBegin(); efi != old_edge->facesEnd(); ++efi)
        if (new_edge->hasIncidentFace(*efi))
        {
          THEA_DEBUG << getName() << ": Can't replace an edge with another on the same face";
          return false;
        }

      bool swap_endpts = false;
      if (old_v[0] == new_v[0] || old_v[1] == new_v[1])
        swap_endpts = false;
      else if (old_v[0] == new_v[1] || old_v[1] == new_v[0])
        swap_endpts = true;
      else
      {
        Real err0 = (old_v[0]->getPosition() - new_v[0]->getPosition()).squaredNorm();
        Real err1 = (old_v[0]->getPosition() - new_v[1]->getPosition()).squaredNorm();
        swap_endpts = (err1 < err0);
      }

      if (swap_endpts)
      {
        Vertex * tmp = new_v[0];
        new_v[0] = new_v[1];
        new_v[1] = tmp;
      }

      // Replace the endpoints
      for (int i = 0; i < 2; ++i)
        if (old_v[i] != new_v[i])
        {
          replaceVertex(old_v[i], new_v[i]);
        }

      // Now there is a double edge, which we will remove
      new_v[0]->removeEdge(old_edge);
      new_v[1]->removeEdge(old_edge);
      for (auto fi = old_edge->facesBegin(); fi != old_edge->facesEnd(); ++fi)
      {
        (*fi)->replaceEdge(old_edge, new_edge);
        if (!new_edge->hasIncidentFace(*fi))
          new_edge->addFace(*fi);
      }

      old_edge->faces.clear();

      invalidateGpuBuffers();
      return true;
    }

    /** Unmark all vertices. */
    void unmarkAllVertices()
    {
      for (auto vi = verticesBegin(); vi != verticesEnd(); ++vi)
        vi->unmark();
    }

    /** Unmark all edges. */
    void unmarkAllEdges()
    {
      for (auto ei = edgesBegin(); ei != edgesEnd(); ++ei)
        ei->unmark();
    }

    /** Unmark all faces. */
    void unmarkAllFaces()
    {
      for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
        fi->unmark();
    }

    /** Remove all references to a face from its adjacent elements. */
    void unlinkFace(Face * face)
    {
      debugAssertM(face, getNameStr() + ": Can't unlink null face");

      for (auto fvi = face->vertices.begin(); fvi != face->vertices.end(); ++fvi)
        (*fvi)->removeFace(face);

      for (auto fei = face->edges.begin(); fei != face->edges.end(); ++fei)
        (*fei)->removeFace(face);

      invalidateGpuBuffers();
    }

    /** Replace a higher-degree face with multiple triangular faces. */
    bool replaceFaceWithTriangulation(Face * face, Vertex ** face_vertices, intx num_tris, intx * tri_indices)
    {
      debugAssertM(face, getNameStr() + ": Can't replace null face with triangulation");

      Vertex * tri_vertices[3];
      for (intx i = 0; i < num_tris; ++i)
      {
        for (int j = 0; j < 3; ++j)
          tri_vertices[j] = face_vertices[tri_indices[3 * i + j]];

        if (i == 0)
        {
          // Repurpose the existing face as the first triangle. Avoids having to remove the face from the mesh (linear-time).
          unlinkFace(face);
          if (!initFace(face, &tri_vertices[0], &tri_vertices[0] + 3))
            return false;
        }
        else
        {
          // Construct a new triangular face
          if (!addFace(&tri_vertices[0], &tri_vertices[0] + 3))
            return false;
        }
      }

      invalidateGpuBuffers();
      return true;
    }

    /** Draw the mesh in immediate rendering mode. */
    int8 drawImmediate(IRenderSystem & render_system, IRenderOptions const & options) const;

    /**
     * Utility function to draw a face. Must be enclosed in the appropriate
     * IRenderSystem::beginPrimitive()/IRenderSystem::endPrimitive() block.
     */
    void drawFace(Face const & face, IRenderSystem & render_system, IRenderOptions const & options) const
    {
      if (!options.useVertexNormals() && options.sendNormals())
        face.drawNormal(render_system, options);

      if (!options.useVertexData())
        face.attr().draw(render_system, options);

      for (auto vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      {
        Vertex const & vertex = **vi;

        if (options.useVertexNormals() && options.sendNormals())
          vertex.drawNormal(render_system, options);

        if (options.useVertexData())  // vertex attributes (per-vertex color, texcoord etc)
          vertex.attr().draw(render_system, options);

        vertex.drawPosition(render_system, options);  // finally send the position of the vertex
      }
    }

    /** Check if a packed array is synchronized with the mesh or not. */
    bool packedArrayIsValid(BufferId buffer) const { return (changed_packed & (int)buffer) == 0; }

    /** Mark a specific packed array as being synchronized with the mesh. */
    void setPackedArrayIsValid(BufferId buffer) const { changed_packed &= (~(int)buffer); }

    /** Clear the set of changed packed arrays. */
    void setAllPackedArraysAreValid() const { changed_packed = 0; }

    /** Check if a GPU buffer is synchronized with the mesh or not. */
    bool gpuBufferIsValid(BufferId buffer) const { return (changed_buffers & (int)buffer) == 0; }

    /** Clear the set of changed buffers. */
    void setAllGpuBuffersAreValid() { setAllPackedArraysAreValid(); changed_buffers = 0; }

    /** Upload GPU resources to the graphics system. */
    bool uploadToGraphicsSystem(IRenderSystem & render_system);

    /** Draw the mesh in GPU-buffered rendering mode. */
    int8 drawBuffered(IRenderSystem & render_system, IRenderOptions const & options) const;

    /** Pack vertex positions densely in an array. */
    void packVertexPositions() const
    {
      if (packedArrayIsValid(BufferId::VERTEX_POSITION)) return;

      packed_vertex_positions.resize(vertices.size());
      size_t i = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_positions[i] = vi->getPosition();

      setPackedArrayIsValid(BufferId::VERTEX_POSITION);
    }

    /** Pack vertex positions densely in an array. */
    void packVertexNormals() const
    {
      if (packedArrayIsValid(BufferId::VERTEX_NORMAL)) return;

      packed_vertex_normals.resize(vertices.size());
      size_t i = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_normals[i] = vi->getNormal();

      setPackedArrayIsValid(BufferId::VERTEX_NORMAL);
    }

    /** Pack vertex colors densely in an array. */
    template < typename VertexT, typename std::enable_if< HasColor<VertexT>::value, int >::type = 0 >
    void packVertexColors() const
    {
      if (packedArrayIsValid(BufferId::VERTEX_COLOR)) return;

      packed_vertex_colors.resize(vertices.size());
      size_t i = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_colors[i] = ColorRgba(vi->attr().getColor());

      setPackedArrayIsValid(BufferId::VERTEX_COLOR);
    }

    /** Clear the array of packed vertex colors (called when vertices don't have attached colors). */
    template < typename VertexT, typename std::enable_if< !HasColor<VertexT>::value, int >::type = 0 >
    void packVertexColors() const
    {
      if (!packedArrayIsValid(BufferId::VERTEX_COLOR))
      {
        packed_vertex_colors.clear();
        setPackedArrayIsValid(BufferId::VERTEX_COLOR);
      }
    }

    /** Pack vertex texture coordinates densely in an array. */
    template < typename VertexT, typename std::enable_if< HasTexCoord<VertexT>::value, int >::type = 0 >
    void packVertexTexCoords() const
    {
      if (packedArrayIsValid(BufferId::VERTEX_TEXCOORD)) return;

      packed_vertex_texcoords.resize(vertices.size());
      size_t i = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_texcoords[i] = vi->attr().getTexCoord();

      setPackedArrayIsValid(BufferId::VERTEX_TEXCOORD);
    }

    /** Clear the array of packed vertex texture coordinates (called when vertices don't have attached texture coordinates). */
    template < typename VertexT, typename std::enable_if< !HasTexCoord<VertexT>::value, int >::type = 0 >
    void packVertexTexCoords() const
    {
      if (!packedArrayIsValid(BufferId::VERTEX_TEXCOORD))
      {
        packed_vertex_texcoords.clear();
        setPackedArrayIsValid(BufferId::VERTEX_TEXCOORD);
      }
    }

    /** Pack face and edge indices densely in an array. */
    void packTopology() const
    {
      if (packedArrayIsValid(BufferId::TOPOLOGY)) return;

      packed_tris.clear();
      packed_quads.clear();

      uint32 packing_index = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
        vi->setPackingIndex(packing_index++);

      has_large_polys = false;
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
      {
        if (fi->isTriangle())
        {
          auto vi = fi->verticesBegin();
          packed_tris.push_back((*(vi++))->getPackingIndex());
          packed_tris.push_back((*(vi++))->getPackingIndex());
          packed_tris.push_back((*vi)->getPackingIndex());
        }
        else if (fi->isQuad())
        {
          auto vi = fi->verticesBegin();
          packed_quads.push_back((*(vi++))->getPackingIndex());
          packed_quads.push_back((*(vi++))->getPackingIndex());
          packed_quads.push_back((*(vi++))->getPackingIndex());
          packed_quads.push_back((*vi)->getPackingIndex());
        }
        else
        {
          has_large_polys = true;
          THEA_DEBUG << getName() << ": Mesh has polygons with 5 or more vertices. These will not be drawn with GPU buffers.";
        }
      }

      if (buffered_wireframe)
      {
        packed_edges.resize(2 * edges.size());
        size_t i = 0;
        for (auto ei = edges.begin(); ei != edges.end(); ++ei, i += 2)
        {
          packed_edges[i    ] = ei->getEndpoint(0)->getPackingIndex();
          packed_edges[i + 1] = ei->getEndpoint(1)->getPackingIndex();
        }
      }
      else
        packed_edges.clear();

      num_tri_indices   =  (intx)packed_tris.size();
      num_quad_indices  =  (intx)packed_quads.size();

      setPackedArrayIsValid(BufferId::TOPOLOGY);
    }

    typedef Array<Vector3>    PositionArray;  ///< Array of vertex positions.
    typedef Array<Vector3>    NormalArray;    ///< Array of normals.
    typedef Array<ColorRgba>  ColorArray;     ///< Array of colors.
    typedef Array<Vector2>    TexCoordArray;  ///< Array of texture coordinates.
    typedef Array<uint32>     IndexArray;     ///< Array of indices.

    FaceList    faces;        ///< Set of mesh faces.
    VertexList  vertices;     ///< Set of mesh vertices.
    EdgeList    edges;        ///< Set of mesh edges.

    intx max_vertex_index;    ///< The largest index of a vertex in the mesh.
    intx max_face_index;      ///< The largest index of a face in the mesh.

    AxisAlignedBox3 bounds;   ///< Mesh bounding box.

    mutable int changed_packed;    ///< A bitwise OR of the flags of the packed arrays that need to be updated.
    mutable bool has_large_polys;  ///< Does the mesh have polygons with more than 4 vertices?
    bool buffered_rendering;       ///< Should the mesh be rendered using GPU buffers?
    bool buffered_wireframe;       ///< Can edges be drawn in GPU-buffered mode?
    int changed_buffers;           ///< A bitwise OR of the flags of the buffers that need to be updated.

    mutable PositionArray  packed_vertex_positions;  ///< Array containing packed set of vertex positions.
    mutable NormalArray    packed_vertex_normals;    ///< Array containing packed set of vertex normals.
    mutable ColorArray     packed_vertex_colors;     ///< Array containing packed set of vertex colors.
    mutable TexCoordArray  packed_vertex_texcoords;  ///< Array containing packed set of vertex texture coordinates.
    mutable IndexArray     packed_tris;              ///< Array containing packed set of triangle indices.
    mutable IndexArray     packed_quads;             ///< Array containing packed set of quad indices.
    mutable IndexArray     packed_edges;             ///< Array containing packed set of edge indices.
    mutable intx           num_tri_indices;          ///< Number of triangle indices in the mesh.
    mutable intx           num_quad_indices;         ///< Number of quad indices in the mesh.

    IBufferPool * buf_pool;          ///< GPU buffer pool.
    IBuffer * vertex_positions_buf;  ///< GPU buffer for vertex positions.
    IBuffer * vertex_normals_buf;    ///< GPU buffer for vertex normals.
    IBuffer * vertex_colors_buf;     ///< GPU buffer for vertex colors.
    IBuffer * vertex_texcoords_buf;  ///< GPU buffer for texture coordinates.
    IBuffer * tris_buf;              ///< GPU buffer for triangle indices.
    IBuffer * quads_buf;             ///< GPU buffer for quad indices.
    IBuffer * edges_buf;             ///< GPU buffer for edges.

    mutable Array<Vertex *> face_vertices;  ///< Internal cache of vertex pointers for a face.

    // Map the packed vertex and index buffers as matrices for the IMesh API
    typedef MatrixMap<3, Eigen::Dynamic, Real,   MatrixLayout::COLUMN_MAJOR>  VertexMatrix;    ///< Wraps vertices as a matrix.
    typedef MatrixMap<3, Eigen::Dynamic, uint32, MatrixLayout::COLUMN_MAJOR>  TriangleMatrix;  ///< Wraps triangles as a matrix.
    typedef MatrixMap<4, Eigen::Dynamic, uint32, MatrixLayout::COLUMN_MAJOR>  QuadMatrix;      ///< Wraps quads as a matrix.

    mutable VertexMatrix    vertex_matrix;  ///< Vertex data as a dense 3xN column-major matrix.
    mutable TriangleMatrix  tri_matrix;     ///< Triangle indices as a dense 3xN column-major matrix.
    mutable QuadMatrix      quad_matrix;    ///< Quad indices as a dense 4xN column-major matrix.

    mutable MatrixWrapper<VertexMatrix>    vertex_wrapper;
    mutable MatrixWrapper<TriangleMatrix>  tri_wrapper;
    mutable MatrixWrapper<QuadMatrix>      quad_wrapper;

}; // class GeneralMesh

template <typename V, typename E, typename F, template <typename T> class A>
bool
GeneralMesh<V, E, F, A>::uploadToGraphicsSystem(IRenderSystem & render_system)
{
  if (!buffered_rendering || changed_buffers == 0) return true;

  if (!gpuBufferIsValid(BufferId::TOPOLOGY))
    invalidateGpuBuffers(BufferId::ALL);

  if (changed_buffers == BufferId::ALL)
  {
    if (buf_pool) buf_pool->reset();

    vertex_positions_buf  =  nullptr;
    vertex_normals_buf    =  nullptr;
    vertex_colors_buf     =  nullptr;
    vertex_texcoords_buf  =  nullptr;
    tris_buf              =  nullptr;
    quads_buf             =  nullptr;
    edges_buf             =  nullptr;

    if (vertices.empty() || (faces.empty() && edges.empty()))
    {
      if (buf_pool)
      {
        render_system.destroyBufferPool(buf_pool);
        buf_pool = nullptr;
      }

      setAllGpuBuffersAreValid();
      return true;
    }

    packArrays();

    static int const PADDING = 32;
    intx vertex_position_bytes = !packed_vertex_positions.empty() ? 3 * 4 * (intx)packed_vertex_positions.size() + PADDING : 0;
    intx vertex_normal_bytes   = !packed_vertex_normals.empty()   ? 3 * 4 * (intx)packed_vertex_normals.size()   + PADDING : 0;
    intx vertex_color_bytes    = !packed_vertex_colors.empty()    ? 4 * 4 * (intx)packed_vertex_colors.size()    + PADDING : 0;
    intx vertex_texcoord_bytes = !packed_vertex_texcoords.empty() ? 2 * 4 * (intx)packed_vertex_texcoords.size() + PADDING : 0;

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
    intx num_bytes = vertex_position_bytes + vertex_normal_bytes + vertex_color_bytes + vertex_texcoord_bytes + PADDING;
#else
    intx tri_bytes   =  !packed_tris.empty()   ?  4 * (intx)packed_tris.size()   +  PADDING : 0;  // uint32
    intx quad_bytes  =  !packed_quads.empty()  ?  4 * (intx)packed_quads.size()  +  PADDING : 0;  // uint32
    intx edge_bytes  =  !packed_edges.empty()  ?  4 * (intx)packed_edges.size()  +  PADDING : 0;  // uint32

    intx num_bytes = vertex_position_bytes
                   + vertex_normal_bytes
                   + vertex_color_bytes
                   + vertex_texcoord_bytes
                   + tri_bytes
                   + quad_bytes
                   + edge_bytes
                   + PADDING;
#endif

    // THEA_CONSOLE << "num_bytes = " << num_bytes;
    // THEA_CONSOLE << "packed_vertex_positions.size() = " << packed_vertex_positions.size();
    // THEA_CONSOLE << "packed_vertex_normals.size() = " << packed_vertex_normals.size();
    // THEA_CONSOLE << "packed_vertex_colors.size() = " << packed_vertex_colors.size();
    // THEA_CONSOLE << "packed_vertex_texcoords.size() = " << packed_vertex_texcoords.size();
    // THEA_CONSOLE << "packed_tris.size() = " << packed_tris.size();
    // THEA_CONSOLE << "packed_quads.size() = " << packed_quads.size();
    // THEA_CONSOLE << "packed_edges.size() = " << packed_edges.size();
    // THEA_CONSOLE << "num_tri_indices = " << num_tri_indices;
    // THEA_CONSOLE << "num_quad_indices = " << num_quad_indices;

    if (buf_pool)
    {
      if (buf_pool->getCapacity() <= num_bytes || buf_pool->getCapacity() > (intx)(1.5 * num_bytes))
      {
        render_system.destroyBufferPool(buf_pool);

        std::string pool_name = getNameStr() + " buffer pool";
        buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
        if (!buf_pool) return false;
      }
      // Else no need to reset buf_pool, we've done it above
    }
    else
    {
      std::string pool_name = getNameStr() + " buffer pool";
      buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
      if (!buf_pool) return false;
    }

    if (!packed_vertex_positions.empty())
    {
      if (!(vertex_positions_buf = buf_pool->createBuffer(vertex_position_bytes))) return false;
      if (!vertex_positions_buf->updateAttributes(0, (int64)packed_vertex_positions.size(), 3, NumericType::REAL,
                                                  &packed_vertex_positions[0])) return false;
    }

    if (!packed_vertex_normals.empty())
    {
      if (!(vertex_normals_buf = buf_pool->createBuffer(vertex_normal_bytes))) return false;
      if (!vertex_normals_buf->updateAttributes(0, (int64)packed_vertex_normals.size(), 3, NumericType::REAL,
                                                &packed_vertex_normals[0])) return false;
    }

    if (!packed_vertex_colors.empty())
    {
      if (!(vertex_colors_buf = buf_pool->createBuffer(vertex_color_bytes))) return false;
      if (!vertex_colors_buf->updateAttributes(0, (int64)packed_vertex_colors.size(), 4, NumericType::REAL,
                                               &packed_vertex_colors[0])) return false;
    }

    if (!packed_vertex_texcoords.empty())
    {
      if (!(vertex_texcoords_buf = buf_pool->createBuffer(vertex_texcoord_bytes))) return false;
      if (!vertex_texcoords_buf->updateAttributes(0, (int64)packed_vertex_texcoords.size(), 2, NumericType::REAL,
                                                  &packed_vertex_texcoords[0])) return false;
    }

#ifndef THEA_GENERAL_MESH_NO_INDEX_ARRAY
    if (!packed_tris.empty())
    {
      if (!(tris_buf = buf_pool->createBuffer(tri_bytes))) return false;
      if (!tris_buf->updateIndices(0, (int64)packed_tris.size(), NumericType::UINT32, &packed_tris[0])) return false;
    }

    if (!packed_quads.empty())
    {
      if (!(quads_buf = buf_pool->createBuffer(quad_bytes))) return false;
      if (!quads_buf->updateIndices(0, (int64)packed_quads.size(), NumericType::UINT32, &packed_quads[0])) return false;
    }

    if (!packed_edges.empty())
    {
      if (!(edges_buf = buf_pool->createBuffer(edge_bytes))) return false;
      if (!edges_buf->updateIndices(0, (int64)packed_edges.size(), NumericType::UINT32, &packed_edges[0])) return false;
    }
#endif
  }
  else
  {
    packArrays();

    if (!gpuBufferIsValid(BufferId::VERTEX_POSITION) && !vertices.empty()
     && !vertex_positions_buf->updateAttributes(0, (int64)packed_vertex_positions.size(), 3, NumericType::REAL,
                                                &packed_vertex_positions[0])) return false;

    if (!gpuBufferIsValid(BufferId::VERTEX_NORMAL) && !vertices.empty()
     && !vertex_normals_buf->updateAttributes(0, (int64)packed_vertex_normals.size(), 3, NumericType::REAL,
                                              &packed_vertex_normals[0])) return false;

    if (!gpuBufferIsValid(BufferId::VERTEX_COLOR) && hasVertexColors()
     && !vertex_colors_buf->updateAttributes(0, (int64)packed_vertex_colors.size(), 4, NumericType::REAL,
                                             &packed_vertex_colors[0])) return false;

    if (!gpuBufferIsValid(BufferId::VERTEX_TEXCOORD) && hasVertexTexCoords()
     && !vertex_texcoords_buf->updateAttributes(0, (int64)packed_vertex_texcoords.size(), 2, NumericType::REAL,
                                                &packed_vertex_texcoords[0])) return false;
  }

  setAllGpuBuffersAreValid();

  return true;
}

template <typename V, typename E, typename F, template <typename T> class A>
int8
GeneralMesh<V, E, F, A>::drawBuffered(IRenderSystem & render_system, IRenderOptions const & options) const
{
  if (options.drawEdges() && !buffered_wireframe)
  { THEA_ERROR << getName() << ": Can't draw mesh edges with GPU-buffered wireframe disabled"; return false; }

  if (!const_cast<GeneralMesh *>(this)->uploadToGraphicsSystem(render_system)) return false;

  if (!vertex_positions_buf) return true;
  if (!options.drawFaces() && !options.drawEdges()) return true;
  if (!options.drawFaces() && edges.size() <= 0) return true;
  if (!options.drawEdges() && faces.size() <= 0) return true;

  render_system.beginIndexedPrimitives();

    render_system.setVertexBuffer  (vertex_positions_buf);
    render_system.setNormalBuffer  (options.sendNormals()      && vertex_normals_buf    ?  vertex_normals_buf    :  nullptr);
    render_system.setColorBuffer   (options.sendColors()       && vertex_colors_buf     ?  vertex_colors_buf     :  nullptr);
    render_system.setTexCoordBuffer(0, options.sendTexCoords() && vertex_texcoords_buf  ?  vertex_texcoords_buf  :  nullptr);

    if (options.drawFaces())
    {
      if (options.drawEdges())
      {
        render_system.pushShapeFlags();
        render_system.setPolygonOffset(true, 2);
      }

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
        if (num_tri_indices > 0)
          render_system.sendIndices(IRenderSystem::Primitive::TRIANGLES, num_tri_indices, &packed_tris[0]);

        if (num_quad_indices > 0)
          render_system.sendIndices(IRenderSystem::Primitive::QUADS, num_quad_indices, &packed_quads[0]);
#else
        if (num_tri_indices > 0)
        {
          render_system.setIndexBuffer(tris_buf);
          render_system.sendIndicesFromBuffer(IRenderSystem::Primitive::TRIANGLES, 0, num_tri_indices);
        }

        if (num_quad_indices > 0)
        {
          render_system.setIndexBuffer(quads_buf);
          render_system.sendIndicesFromBuffer(IRenderSystem::Primitive::QUADS, 0, num_quad_indices);
        }
#endif

        // Finish off with all larger polygons
        if (has_large_polys)
        {
          for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
            if (fi->numEdges() > 4)
            {
              render_system.beginPrimitive(IRenderSystem::Primitive::POLYGON);
                drawFace(*fi, render_system, options);
              render_system.endPrimitive();
            }
        }

      if (options.drawEdges())
        render_system.popShapeFlags();
    }

    if (options.drawEdges())
    {
      render_system.pushShader();
      render_system.pushColorFlags();

        render_system.setShader(nullptr);
        render_system.setColorBuffer(nullptr);
        render_system.setTexCoordBuffer(0, nullptr);
        render_system.setNormalBuffer(nullptr);
        render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
        if (!edges.empty())
          render_system.sendIndices(IRenderSystem::Primitive::LINES, 2 * (int64)edges.size(), &packed_edges[0]);
#else
        if (!edges.empty())
        {
          render_system.setIndexBuffer(edges_buf);
          render_system.sendIndicesFromBuffer(IRenderSystem::Primitive::LINES, 0, 2 * (int64)edges.size());
        }
#endif

      render_system.popColorFlags();
      render_system.popShader();
    }

  render_system.endIndexedPrimitives();

  if (char const * err = render_system.getLastError())
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

template <typename V, typename E, typename F, template <typename T> class A>
int8
GeneralMesh<V, E, F, A>::drawImmediate(IRenderSystem & render_system, IRenderOptions const & options) const
{
  // Three separate passes over the faces is probably faster (TODO: profile) than using Primitive::POLYGON for each face

  if (options.drawFaces())
  {
    if (options.drawEdges())
    {
      render_system.pushShapeFlags();
      render_system.setPolygonOffset(true, 1);
    }

    // First try to render as much stuff using triangles as possible
    render_system.beginPrimitive(IRenderSystem::Primitive::TRIANGLES);
      for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isTriangle()) drawFace(*fi, render_system, options);
    render_system.endPrimitive();

    // Now render all quads
    render_system.beginPrimitive(IRenderSystem::Primitive::QUADS);
      for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isQuad()) drawFace(*fi, render_system, options);
    render_system.endPrimitive();

    // Finish off with all larger polygons
    for (auto fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->numEdges() > 4)
      {
        render_system.beginPrimitive(IRenderSystem::Primitive::POLYGON);
          drawFace(*fi, render_system, options);
        render_system.endPrimitive();
      }

    if (options.drawEdges())
      render_system.popShapeFlags();
  }

  if (options.drawEdges())
  {
    render_system.pushShader();
    render_system.pushColorFlags();

      render_system.setShader(nullptr);
      render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)

      render_system.beginPrimitive(IRenderSystem::Primitive::LINES);
        for (auto ei = edgesBegin(); ei != edgesEnd(); ++ei)
        {
          ei->attr().draw(render_system, options);
          ei->getEndpoint(0)->drawPosition(render_system, options);
          ei->getEndpoint(1)->drawPosition(render_system, options);
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }

  if (char const * err = render_system.getLastError())
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

} // namespace Graphics
} // namespace Thea

#endif
