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

#ifndef __Thea_Graphics_GeneralMesh_hpp__
#define __Thea_Graphics_GeneralMesh_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Colors.hpp"
#include "../List.hpp"
#include "../NamedObject.hpp"
#include "../Polygon3.hpp"
#include "../UnorderedMap.hpp"
#include "GeneralMeshFace.hpp"
#include "GeneralMeshVertex.hpp"
#include "GeneralMeshEdge.hpp"
#include "GraphicsAttributes.hpp"
#include "IncrementalGeneralMeshBuilder.hpp"
#include "DefaultMeshCodecs.hpp"
#include "DrawableObject.hpp"
#include "EdgeWelder.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Graphics {

/**
 * A class for storing meshes with arbitrary topologies. Optionally allows GPU-buffered rendering, which requires the user to
 * manually indicate when mesh contents have changed and need to be resynchronized with the GPU.
 *
 * @todo Automatically invalidate appropriate GPU buffers on all modifications.
 * @todo Add support for GPU-buffered texture coordinates with 1, 3 or 4 dimensions.
 * @todo Instantiate different types of GPU buffers for different types of colors/texture coordinates.
 */
template < typename VertexAttributeT               =  Graphics::NullAttribute,
           typename EdgeAttributeT                 =  Graphics::NullAttribute,
           typename FaceAttributeT                 =  Graphics::NullAttribute,
           template <typename T> class AllocatorT  =  std::allocator >
class /* THEA_API */ GeneralMesh : public virtual NamedObject, public DrawableObject
{
  public:
    THEA_DEF_POINTER_TYPES(GeneralMesh, shared_ptr, weak_ptr)

    /** Mesh type tag. */
    struct GENERAL_MESH_TAG {};

    /**< Vertex of the mesh. */
    typedef GeneralMeshVertex <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Vertex;

    /**< Edge of the mesh. */
    typedef GeneralMeshEdge   <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Edge;

    /**< Face of the mesh. */
    typedef GeneralMeshFace   <VertexAttributeT, EdgeAttributeT, FaceAttributeT, AllocatorT>  Face;

  private:
    typedef TheaList<Vertex, AllocatorT<Vertex> >  VertexList;
    typedef TheaList<Edge,   AllocatorT<Edge>   >  EdgeList;
    typedef TheaList<Face,   AllocatorT<Face>   >  FaceList;

  public:
    typedef typename VertexList::iterator        VertexIterator;       ///< Iterator over vertices.
    typedef typename VertexList::const_iterator  VertexConstIterator;  ///< Const iterator over vertices.
    typedef typename EdgeList::iterator          EdgeIterator;         ///< Iterator over edges.
    typedef typename EdgeList::const_iterator    EdgeConstIterator;    ///< Const iterator over edges.
    typedef typename FaceList::iterator          FaceIterator;         ///< Iterator over faces.
    typedef typename FaceList::const_iterator    FaceConstIterator;    ///< Const iterator over faces.

    /** Identifiers for the various buffers (enum class). */
    struct BufferID
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

      THEA_ENUM_CLASS_BODY(BufferID)

    }; // struct BufferID

    /** Constructor. */
    GeneralMesh(std::string const & name = "AnonymousMesh")
    : NamedObject(name),
      buffered_rendering(false),
      buffered_wireframe(false),
      changed_buffers(BufferID::ALL),
      has_large_polys(false),
      num_tri_indices(0),
      num_quad_indices(0),
      var_area(NULL),
      vertex_positions_var(NULL),
      vertex_normals_var(NULL),
      vertex_colors_var(NULL),
      vertex_texcoords_var(NULL),
      tris_var(NULL),
      quads_var(NULL),
      edges_var(NULL)
    {}

    /**
     * Copy constructor. Creates a deep copy of the mesh (including copies of the attributes). <b>Currently not implemented.</b>
     */
    GeneralMesh(GeneralMesh const & src) : NamedObject(src)
    {
      throw Error("GeneralMesh: Copy constructor not currently implemented");
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

    /** Deletes all data in the mesh. */
    void clear()
    {
      vertices.clear();
      edges.clear();
      faces.clear();
      bounds = AxisAlignedBox3();

      packed_vertex_positions.clear();
      packed_vertex_normals.clear();
      packed_vertex_colors.clear();
      packed_vertex_texcoords.clear();
      packed_tris.clear();
      packed_quads.clear();
      packed_edges.clear();

      invalidateGPUBuffers();
      has_large_polys = false;
    }

    /**
     * Make an exact copy of the mesh, optionally returning mapping from source to destination vertices/edges/faces. Previous
     * data in the maps is <b>not</b> cleared.
     */
    void copyTo(GeneralMesh & dst,
                TheaUnorderedMap<Vertex const *, Vertex *> * vertex_map = NULL,
                TheaUnorderedMap<Edge const *, Edge *> * edge_map = NULL,
                TheaUnorderedMap<Face const *, Face *> * face_map = NULL) const
    {
      typedef TheaUnorderedMap<Vertex const *, Vertex *> VertexMap;
      typedef TheaUnorderedMap<Edge   const *, Edge   *> EdgeMap;
      typedef TheaUnorderedMap<Face   const *, Face   *> FaceMap;

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

      dst.bounds = bounds;
      dst.buffered_rendering = buffered_rendering;
      dst.buffered_wireframe = buffered_wireframe;
      dst.has_large_polys = has_large_polys;
    }

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const { return vertices.empty() && faces.empty() && edges.empty(); }

    /** Get the number of vertices. */
    long numVertices() const { return (long)vertices.size(); };

    /** Get the number of edges. */
    long numEdges() const { return (long)edges.size(); };

    /** Get the number of faces. */
    long numFaces() const { return (long)faces.size(); };

    /** Compute the number of triangles in the mesh. */
    long numTriangles() const
    {
      long rval = 0;
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isTriangle())
          rval++;

      return rval;
    }

    /** Compute the number of quads in the mesh. */
    long numQuads() const
    {
      long rval = 0;
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isQuad())
          rval++;

      return rval;
    }

    /** Do the mesh vertices have attached colors? */
    bool hasVertexColors() const { return HasColor<Vertex>::value; }

    /** Do the mesh vertices have attached texture coordinates? */
    bool hasVertexTexCoords() const { return HasTexCoord<Vertex>::value; }

    /**
     * Add a vertex to the mesh, with optional precomputed normal, color and texture coordinates. Automatically calls
     * invalidateGPUBuffers() to schedule a resync with the GPU.
     *
     * @return A pointer to the newly created vertex on success, null on failure.
     */
    Vertex * addVertex(Vector3 const & point, Vector3 const * normal = NULL, ColorRGBA const * color = NULL,
                       Vector2 const * texcoord = NULL)
    {
      if (normal)
        vertices.push_back(Vertex(point, *normal));
      else
        vertices.push_back(Vertex(point));

      invalidateGPUBuffers();

      Vertex * vertex = &(*vertices.rbegin());
      if (color)     setVertexColor<Vertex>(vertex, *color);
      if (texcoord)  setVertexTexCoord<Vertex>(vertex, *texcoord);

      return vertex;
    }

    /**
     * Add a face to the mesh, specified by the sequence of vertices obtained by dereferencing [vbegin, vend).
     * VertexInputIterator must dereference to a pointer to a Vertex. Unless the mesh is already in an inconsistent state,
     * failure to add the face will not affect the mesh.
     *
     * Automatically calls invalidateGPUBuffers() to schedule a resync with the GPU.
     *
     * @return A pointer to the newly created face, or null on error.
     */
    template <typename VertexInputIterator>
    Face * addFace(VertexInputIterator vbegin, VertexInputIterator vend)
    {
      // Create the (initially empty) face
      faces.push_back(Face());
      Face * face = &(*faces.rbegin());

      // Initialize the face
      face = initFace(face, vbegin, vend);
      if (!face)
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
      for (FaceIterator fi = faces.begin(); fi != faces.end(); ++fi)
        if (&(*fi) == face)
          return removeFace(fi);

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
      for (EdgeConstIterator ei = edges.begin(); ei != edges.end(); ++ei)
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
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
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

      for (EdgeIterator ei = edges.begin(); ei != edges.end(); ++ei)
      {
        if (!ei->isBoundary())
          continue;

        if (ei->isSelfLoop())
          continue;

        Edge * edge = &(*ei);
        Edge * existing = (Edge *)welder.getUndirectedEdge(edge->getEndpoint(0)->getPosition(),
                                                           edge->getEndpoint(1)->getPosition());
        bool can_seal = (existing != NULL);
        if (can_seal)
        {
          // We can seal only if the two edges don't share a face
          for (typename Edge::FaceConstIterator fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
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
        for (typename Vertex::EdgeIterator ei = old_vertex->edgesBegin(); ei != old_vertex->edgesEnd(); ++ei)
        {
          (*ei)->replaceVertex(old_vertex, new_vertex);

          if (!new_vertex->hasIncidentEdge(*ei))
            new_vertex->addEdge(*ei);
        }

        for (typename Vertex::FaceIterator fi = old_vertex->facesBegin(); fi != old_vertex->facesEnd(); ++fi)
        {
          (*fi)->replaceVertex(old_vertex, new_vertex);

          if (!new_vertex->hasIncidentFace(*fi))
            new_vertex->addFace(*fi);
        }

        old_vertex->edges.clear();
        old_vertex->faces.clear();
      }
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
      for (typename Edge::FaceIterator fi = edge->facesBegin(); fi != edge->facesEnd(); )
      {
        Face * face = *fi;
        typename Face::EdgeIterator fei = face->edgesBegin();
        typename Face::VertexIterator fvi = face->verticesBegin();

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
              typename Face::VertexIterator next = fvi; ++next;
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
      for (typename Vertex::EdgeIterator ei = vertex_to_remove->edgesBegin(); ei != vertex_to_remove->edgesEnd(); ++ei)
      {
        // We do the replacement even for the edge that is getting collapsed, since if a face repeats this edge it is going to
        // be present in the output (as a self-loop)

        (*ei)->replaceVertex(vertex_to_remove, vertex_to_preserve);

        if (!(*ei)->hasEndpoint(vertex_to_preserve))
          vertex_to_preserve->addEdge(*ei);
      }

      for (typename Vertex::FaceIterator fi = vertex_to_remove->facesBegin(); fi != vertex_to_remove->facesEnd(); ++fi)
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
        return NULL;

      if (edge->isSelfLoop())
      {
        THEA_DEBUG << getName() << ": Can't split self-loop edge";
        return NULL;
      }

      if (edge->hasEndpoint(vertex))
      {
        THEA_DEBUG << getName() << ": Can't split edge at existing endpoint";
        return NULL;
      }

      for (typename Edge::FaceConstIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        if ((*efi)->hasVertex(vertex))
        {
          THEA_DEBUG << getName() << ": Can't split edge at vertex on the same face";
          return NULL;
        }

      Vertex * old_e1 = edge->getEndpoint(1);
      edge->setEndpoint(1, vertex);
      edges.push_back(Edge(vertex, old_e1));
      Edge * new_edge = &edges.back();

      vertex->addEdge(edge);
      vertex->addEdge(new_edge);

      old_e1->removeEdge(edge);
      old_e1->addEdge(new_edge);

      for (typename Edge::FaceIterator fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
      {
        new_edge->addFace(*fi);

        if (!vertex->hasIncidentFace(*fi))
          vertex->addFace(*fi);
      }

      // Now everything's ok except except the references in the faces
      for (typename Edge::FaceIterator fi = edge->facesBegin(); fi != edge->facesEnd(); ++fi)
      {
        Face * face = *fi;

        typename Face::EdgeIterator ei = face->edgesBegin();
        typename Face::VertexIterator vi = face->verticesBegin();

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
              typename Face::EdgeIterator next_ei = ei; ++next_ei;
              face->edges.insert(next_ei, new_edge);
            }

            // Vertex goes in the middle, i.e. always after the current vertex
            typename Face::VertexIterator next_vi = vi; ++next_vi;
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

      return new_edge;
    }

    /** Split an edge by creating a new vertex along its length at a given position. */
    Vertex * splitEdge(Edge * edge, Vector3 const & p)
    {
      Real s = (p - edge->getEndpoint(1)->getPosition()).length();
      Real t = (p - edge->getEndpoint(0)->getPosition()).length();
      Real sum = s + t;
      s /= sum;
      t /= sum;
      Vector3 n = s * edge->getEndpoint(0)->getNormal() + t * edge->getEndpoint(1)->getNormal();
      Vertex * new_vx = addVertex(p, &n);

      if (!splitEdge(edge, new_vx))  // should generally never happen
      {
        vertices.pop_back();  // remove the vertex we just added
        return NULL;
      }

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
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
      {
        if (!ei->isSelfLoop())
          continue;

        Edge * edge = &(*ei);
        for (typename Edge::FaceIterator efi = edge->facesBegin(); efi != edge->facesEnd(); )
        {
          Face * face = *efi;

          typename Face::EdgeIterator fei = face->edgesBegin();
          typename Face::VertexIterator fvi = face->verticesBegin();
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
    }

    /** Remove all empty faces (fewer than 3 edges) from the mesh. */
    void removeDegenerateFaces()
    {
      for (FaceIterator fi = faces.begin(); fi != facesEnd(); )
      {
        if (fi->numVertices() < 3)
        {
          Face * face = &(*fi);

          for (typename Face::VertexIterator vi = face->verticesBegin(); vi != face->verticesEnd(); ++vi)
            (*vi)->removeFace(face);

          for (typename Face::EdgeIterator ei = face->edgesBegin(); ei != face->edgesEnd(); ++ei)
            (*ei)->removeFace(face);

          fi = faces.erase(fi);
        }
        else
          ++fi;
      }

      removeIsolatedEdges();  // also removes isolated vertices
    }

    /** Remove all isolated edges (no incident faces) from the mesh. */
    void removeIsolatedEdges()
    {
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); )
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
    }

    /** Remove all isolated vertices (no incident edges or faces) from the mesh. */
    void removeIsolatedVertices()
    {
      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); )
      {
        if (vi->faces.empty() && vi->edges.empty())
          vi = vertices.erase(vi);
        else
          ++vi;
      }
    }

    /** Remove all vertices marked for deletion by other operations. */
    void removeMarkedVertices()
    {
      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); )
      {
        if (vi->isMarked())
          vi = vertices.erase(vi);
        else
          ++vi;
      }
    }

    /** Remove all edges marked for deletion by other operations. */
    void removeMarkedEdges()
    {
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); )
      {
        if (ei->isMarked())
          ei = edges.erase(ei);
        else
          ++ei;
      }
    }

    /** Remove all faces marked for deletion by other operations. */
    void removeMarkedFaces()
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); )
      {
        if (fi->isMarked())
          fi = faces.erase(fi);
        else
          ++fi;
      }
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
    long triangulate(Real epsilon = -1)
    {
      long orig_num_faces = numFaces();
      long num_visited_faces = 0;
      long num_triangulated_faces = 0;

      // New faces will be added to the end of the face list, so we can just keep track of when we've processed the original
      // number of faces
      for (FaceIterator fi = facesBegin(); num_visited_faces < orig_num_faces; ++fi, ++num_visited_faces)
        if (fi->numVertices() > 3)
        {
          long nt = triangulate(&(*fi), epsilon);
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
    long triangulate(Face * face, Real epsilon = -1)
    {
      if (!face)
        return 0;

      if (face->numVertices() <= 3)
        return 1;

      long ntris = 0;
      if (face->numVertices() == 4)
      {
        Vertex * face_vertices[4];
        {
          array_size_t i = 0;
          for (typename Face::VertexConstIterator fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi, ++i)
            face_vertices[i] = *fvi;
        }

        long tri_indices[6];
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
        TheaArray<Vertex *> face_vertices((array_size_t)face->numVertices());
        Polygon3 poly;
        {
          array_size_t i = 0;
          for (typename Face::VertexConstIterator fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi, ++i)
          {
            face_vertices[i] = *fvi;
            poly.addVertex((*fvi)->getPosition());
          }
        }

        TheaArray<long> tri_indices;
        ntris = poly.triangulate(tri_indices, epsilon);
        if (ntris >= 1)
        {
          if (!replaceFaceWithTriangulation(face, &face_vertices[0], ntris, &tri_indices[0]))
            return -1;
        }
      }

      // If the triangulation generated no triangles, the face is degenerate but we must replace it with SOME triangle
      if (ntris <= 0)
      {
        Vertex * face_vertices[3];
        typename Face::VertexConstIterator fvi = face->verticesBegin();
        face_vertices[0] = *(fvi++);
        face_vertices[1] = *(fvi++);
        face_vertices[2] = *(fvi);

        long tri_indices[3] = { 0, 1, 2 };
        ntris = 1;

        if (!replaceFaceWithTriangulation(face, &face_vertices[0], ntris, &tri_indices[0]))
          return -1;
      }

      return ntris;
    }

    /**
     * Check if GPU-buffered rendering is on or off. If it is on, you <b>must manually call</b> invalidateGPUBuffers() every
     * time the mesh changes, to make sure the GPU buffers are update when the mesh is next rendered.
     *
     * @see setGPUBufferedRendering()
     */
    bool renderingIsGPUBuffered() const { return buffered_rendering; }

    /**
     * Turn GPU-buffered rendering on/off. If you enable this function, you <b>must manually call</b> invalidateGPUBuffers()
     * every time the mesh changes, to make sure the GPU buffers are synchronized when the mesh is next rendered.
     *
     * @see renderingIsGPUBuffered()
     */
    void setGPUBufferedRendering(bool value)
    {
      if (value == buffered_rendering)
        return;

      buffered_rendering = value;
      invalidateGPUBuffers();
    }

    /** Invalidate part or all of the current GPU data for the mesh. */
    void invalidateGPUBuffers(int changed_buffers_ = BufferID::ALL) { changed_buffers |= changed_buffers_; }

    /**
     * Enable/disable drawing the edges of the mesh in GPU-buffered mode. Enabling this function will <b>not</b> draw any edges
     * unless you turn on the appropriate RenderOptions flag. The edges will be uploaded to the graphics system on the next call
     * to uploadToGraphicsSystem().
     *
     * Wireframe mode is initially disabled to save video memory. This flag is ignored in non-buffered (immediate) mode.
     *
     * @see wireframeIsGPUBuffered()
     */
    void setGPUBufferedWireframe(bool value)
    {
      if (value == buffered_wireframe)
        return;

      buffered_wireframe = value;
      invalidateGPUBuffers();
    }

    /**
     * Check if wireframe drawing is enabled in GPU-buffered mode. This is initially disabled to save video memory. This flag is
     * ignored in non-buffered (immediate) mode.
     *
     * @see setGPUBufferedWireframe()
     */
    bool wireframeIsGPUBuffered() const { return buffered_wireframe; }

    void uploadToGraphicsSystem(RenderSystem & render_system);

    void draw(RenderSystem & render_system, RenderOptions const & options = RenderOptions::defaults()) const
    {
      if (buffered_rendering)
        drawBuffered(render_system, options);
      else
        drawImmediate(render_system, options);
    }

    void updateBounds()
    {
      bounds = AxisAlignedBox3();
      for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        bounds.merge(vi->getPosition());
    }

    AxisAlignedBox3 const & getBounds() const { return bounds; }

  private:
    /**
     * Initialize a pre-constructed face, which will be assigned the sequence of vertices obtained by dereferencing
     * [vbegin, vend). VertexInputIterator must dereference to a pointer to a Vertex. Unless the mesh is already in an
     * inconsistent state, failure to add the face will not affect the mesh.
     *
     * Automatically calls invalidateGPUBuffers() to schedule a resync with the GPU.
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
      array_size_t num_verts = 0;
      Vector3 v[3];
      for (VertexInputIterator vi = vbegin; vi != vend; ++vi, ++num_verts)
      {
        debugAssertM(*vi, getNameStr() + ": Null vertex pointer specified for new face");
        if (num_verts < 3) v[num_verts] = (*vi)->getPosition();
      }

      if (num_verts < 3)
      {
        THEA_WARNING << getName() << ": Skipping face -- too few vertices (" << num_verts << ')';
        return NULL;
      }

      face->clear();

      // Add the loop of vertices to the face
      VertexInputIterator next = vbegin;
      for (VertexInputIterator vi = next++ ; vi != vend; ++vi, ++next)
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
      for (typename Face::VertexIterator fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi)
        (*fvi)->addFaceNormal(face->getNormal());  // weight by face area?

      invalidateGPUBuffers();

      return face;
    }

    /** Set vertex color. */
    template <typename VertexT>
    void setVertexColor(VertexT * vertex, ColorRGBA const & color,
                        typename boost::enable_if< HasColor<VertexT> >::type * dummy = NULL)
    {
      vertex->attr().setColor(color);
    }

    /** Set vertex color (no-op, called if vertex does not have color attribute). */
    template <typename VertexT>
    void setVertexColor(VertexT * vertex, ColorRGBA const & color,
                        typename boost::disable_if< HasColor<VertexT> >::type * dummy = NULL)
    {}

    /** Set vertex texture coordinates. */
    template <typename VertexT>
    void setVertexTexCoord(VertexT * vertex, Vector2 const & texcoord,
                           typename boost::enable_if< HasTexCoord<VertexT> >::type * dummy = NULL)
    {
      vertex->attr().setTexCoord(texcoord);
    }

    /** Set vertex texture coordinates (no-op, called if vertex does not have texture coordinate attribute). */
    template <typename VertexT>
    void setVertexTexCoord(VertexT * vertex, Vector2 const & texcoord,
                           typename boost::disable_if< HasTexCoord<VertexT> >::type * dummy = NULL)
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
      *next_edge = NULL;
      *next_face = NULL;

      // Find the next face around the vertex
      for (typename Edge::FaceConstIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        if (*efi != face)
        {
          *next_face = *efi;
          break;
        }

      // Have a we hit a boundary?
      if (!(*next_face))
        return false;

      // Find the next edge around the vertex
      for (typename Face::EdgeConstIterator fei = (*next_face)->edgesBegin(); fei != (*next_face)->edgesEnd(); ++fei)
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
      for (typename Vertex::EdgeConstIterator vei = vertex->edgesBegin(); vei != vertex->edgesEnd(); ++vei)
      {
        (*vei)->clearAllInternalBits();

        if ((*vei)->isBoundary())
          first_edge = *vei;
      }

      // Check if the edge is isolated
      if (first_edge->numFaces() <= 0)
        return false;

      // Start from the edge found above and an incident face, and step around the vertex
      Edge const * edge = first_edge;
      Face const * face = edge->isBoundary() ? NULL : *edge->facesBegin();  // if we're starting from a boundary edge, the
                                                                            // initial "face" is the empty space before the
                                                                            // first edge
      long num_visited_edges = 0;
      while (true)
      {
        if (++num_visited_edges >= vertex->numEdges())
          break;

        // Step around the vertex
        Edge const * next_edge = NULL;
        Face const * next_face = NULL;
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

      for (typename Edge::FaceConstIterator efi = old_edge->facesBegin(); efi != old_edge->facesEnd(); ++efi)
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
        Real err0 = (old_v[0]->getPosition() - new_v[0]->getPosition()).squaredLength();
        Real err1 = (old_v[0]->getPosition() - new_v[1]->getPosition()).squaredLength();
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
      for (typename Edge::FaceIterator fi = old_edge->facesBegin(); fi != old_edge->facesEnd(); ++fi)
      {
        (*fi)->replaceEdge(old_edge, new_edge);
        if (!new_edge->hasIncidentFace(*fi))
          new_edge->addFace(*fi);
      }

      old_edge->faces.clear();

      return true;
    }

    /** Unmark all vertices. */
    void unmarkAllVertices()
    {
      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
        vi->unmark();
    }

    /** Unmark all edges. */
    void unmarkAllEdges()
    {
      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        ei->unmark();
    }

    /** Unmark all faces. */
    void unmarkAllFaces()
    {
      for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        fi->unmark();
    }

    /** Remove all references to a face from its adjacent elements. */
    void unlinkFace(Face * face)
    {
      debugAssertM(face, getNameStr() + ": Can't unlink null face");

      for (typename Face::VertexIterator fvi = face->vertices.begin(); fvi != face->vertices.end(); ++fvi)
        (*fvi)->removeFace(face);

      for (typename Face::EdgeIterator fei = face->edges.begin(); fei != face->edges.end(); ++fei)
        (*fei)->removeFace(face);
    }

    /** Replace a higher-degree face with multiple triangular faces. */
    bool replaceFaceWithTriangulation(Face * face, Vertex ** face_vertices, long num_tris, long * tri_indices)
    {
      debugAssertM(face, getNameStr() + ": Can't replace null face with triangulation");

      Vertex * tri_vertices[3];
      for (long i = 0; i < num_tris; ++i)
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

      return true;
    }

    /** Draw the mesh in immediate rendering mode. */
    void drawImmediate(RenderSystem & render_system, RenderOptions const & options) const;

    /**
     * Utility function to draw a face. Must be enclosed in the appropriate
     * RenderSystem::beginPrimitive()/RenderSystem::endPrimitive() block.
     */
    void drawFace(Face const & face, RenderSystem & render_system, RenderOptions const & options) const
    {
      if (!options.useVertexNormals() && options.sendNormals())
        face.drawNormal(render_system, options);

      if (!options.useVertexData())
        face.attr().draw(render_system, options);

      for (typename Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      {
        Vertex const & vertex = **vi;

        if (options.useVertexNormals() && options.sendNormals())
          vertex.drawNormal(render_system, options);

        if (options.useVertexData())  // vertex attributes (per-vertex color, texcoord etc)
          vertex.attr().draw(render_system, options);

        vertex.drawPosition(render_system, options);  // finally send the position of the vertex
      }
    }

    /** Check if a GPU buffer is synchronized with the mesh or not. */
    bool gpuBufferIsValid(BufferID buffer) const { return (changed_buffers & (int)buffer) == 0; }

    /** Clear the set of changed buffers. */
    void allGPUBuffersAreValid() { changed_buffers = 0; }

    /** Draw the mesh in GPU-buffered rendering mode. */
    void drawBuffered(RenderSystem & render_system, RenderOptions const & options) const;

    /** Pack vertex positions densely in an array. */
    void packVertexPositions()
    {
      packed_vertex_positions.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_positions[i] = vi->getPosition();
    }

    /** Pack vertex positions densely in an array. */
    void packVertexNormals()
    {
      packed_vertex_normals.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_normals[i] = vi->getNormal();
    }

    /** Pack vertex colors densely in an array. */
    template <typename VertexT>
    void packVertexColors(typename boost::enable_if< HasColor<VertexT> >::type * dummy = NULL)
    {
      packed_vertex_colors.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_colors[i] = ColorRGBA(vi->attr().getColor());
    }

    /** Clear the array of packed vertex colors (called when vertices don't have attached colors). */
    template <typename VertexT>
    void packVertexColors(typename boost::disable_if< HasColor<VertexT> >::type * dummy = NULL)
    {
      packed_vertex_colors.clear();
    }

    /** Pack vertex texture coordinates densely in an array. */
    template <typename VertexT>
    void packVertexTexCoords(typename boost::enable_if< HasTexCoord<VertexT> >::type * dummy = NULL)
    {
      packed_vertex_texcoords.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_texcoords[i] = vi->attr().getTexCoord();
    }

    /** Clear the array of packed vertex texture coordinates (called when vertices don't have attached texture coordinates). */
    template <typename VertexT>
    void packVertexTexCoords(typename boost::disable_if< HasTexCoord<VertexT> >::type * dummy = NULL)
    {
      packed_vertex_texcoords.clear();
    }

    /** Pack face and edge indices densely in an array. */
    void packTopology()
    {
      packed_tris.clear();
      packed_quads.clear();

      uint32 index = 0;
      for (VertexIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
        vi->setIndex(index++);

      has_large_polys = false;
      for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
      {
        if (fi->isTriangle())
        {
          typename Face::VertexConstIterator vi = fi->verticesBegin();
          packed_tris.push_back((*(vi++))->getIndex());
          packed_tris.push_back((*(vi++))->getIndex());
          packed_tris.push_back((*vi)->getIndex());
        }
        else if (fi->isQuad())
        {
          typename Face::VertexConstIterator vi = fi->verticesBegin();
          packed_quads.push_back((*(vi++))->getIndex());
          packed_quads.push_back((*(vi++))->getIndex());
          packed_quads.push_back((*(vi++))->getIndex());
          packed_quads.push_back((*vi)->getIndex());
        }
        else
        {
          has_large_polys = true;
          THEA_WARNING << getName() << ": Mesh has polygons with 5 or more vertices. These will not be drawn with GPU buffers.";
        }
      }

      if (buffered_wireframe)
      {
        packed_edges.resize(2 * edges.size());
        array_size_t i = 0;
        for (EdgeConstIterator ei = edges.begin(); ei != edges.end(); ++ei, i += 2)
        {
          packed_edges[i    ] = ei->getEndpoint(0)->getIndex();
          packed_edges[i + 1] = ei->getEndpoint(1)->getIndex();
        }
      }
      else
        packed_edges.clear();

      num_tri_indices   =  (long)packed_tris.size();
      num_quad_indices  =  (long)packed_quads.size();
    }

    typedef TheaArray<Vector3>  PositionArray;  ///< Array of vertex positions.
    typedef TheaArray<Vector3>  NormalArray;    ///< Array of normals.
    typedef TheaArray<ColorRGBA>   ColorArray;     ///< Array of colors.
    typedef TheaArray<Vector2>  TexCoordArray;  ///< Array of texture coordinates.
    typedef TheaArray<uint32>   IndexArray;     ///< Array of indices.

    FaceList         faces;     ///< Set of mesh faces.
    VertexList       vertices;  ///< Set of mesh vertices.
    EdgeList         edges;     ///< Set of mesh edges.

    AxisAlignedBox3  bounds;    ///< Mesh bounding box.

    bool buffered_rendering;  ///< Should the mesh be rendered using GPU buffers?
    bool buffered_wireframe;  ///< Can edges be drawn in GPU-buffered mode?
    int changed_buffers;      ///< A bitwise OR of the flags of the buffers that have changed.
    bool has_large_polys;     ///< Does the mesh have polygons with more than 4 vertices?

    PositionArray  packed_vertex_positions;  ///< Array containing packed set of vertex positions.
    NormalArray    packed_vertex_normals;    ///< Array containing packed set of vertex normals.
    ColorArray     packed_vertex_colors;     ///< Array containing packed set of vertex colors.
    TexCoordArray  packed_vertex_texcoords;  ///< Array containing packed set of vertex texture coordinates.
    IndexArray     packed_tris;              ///< Array containing packed set of triangle indices.
    IndexArray     packed_quads;             ///< Array containing packed set of quad indices.
    IndexArray     packed_edges;             ///< Array containing packed set of edge indices.
    long           num_tri_indices;          ///< Number of triangles in the mesh.
    long           num_quad_indices;         ///< Number of quads in the mesh.

    VARArea * var_area;          ///< GPU buffer area.
    VAR * vertex_positions_var;  ///< GPU buffer for vertex positions.
    VAR * vertex_normals_var;    ///< GPU buffer for vertex normals.
    VAR * vertex_colors_var;     ///< GPU buffer for vertex colors.
    VAR * vertex_texcoords_var;  ///< GPU buffer for texture coordinates.
    VAR * tris_var;              ///< GPU buffer for triangle indices.
    VAR * quads_var;             ///< GPU buffer for quad indices.
    VAR * edges_var;             ///< GPU buffer for edges.

    mutable TheaArray<Vertex *> face_vertices;  ///< Internal cache of vertex pointers for a face.

}; // class GeneralMesh

template <typename V, typename E, typename F, template <typename T> class A>
inline void
GeneralMesh<V, E, F, A>::uploadToGraphicsSystem(RenderSystem & render_system)
{
  if (!buffered_rendering || changed_buffers == 0) return;

  if (!gpuBufferIsValid(BufferID::TOPOLOGY))
    changed_buffers = BufferID::ALL;

  if (changed_buffers == BufferID::ALL)
  {
    if (var_area) var_area->reset();

    vertex_positions_var  =  NULL;
    vertex_normals_var    =  NULL;
    vertex_colors_var     =  NULL;
    vertex_texcoords_var  =  NULL;
    tris_var              =  NULL;
    quads_var             =  NULL;
    edges_var             =  NULL;

    if (vertices.empty() || (faces.empty() && edges.empty()))
    {
      if (var_area)
      {
        render_system.destroyVARArea(var_area);
        var_area = NULL;
      }

      allGPUBuffersAreValid();
      return;
    }

    static int const PADDING = 32;

    packVertexPositions();
    packVertexNormals();
    packVertexColors<Vertex>();
    packVertexTexCoords<Vertex>();

    long vertex_position_bytes = !packed_vertex_positions.empty() ? 3 * 4 * (long)packed_vertex_positions.size() + PADDING : 0;
    long vertex_normal_bytes   = !packed_vertex_normals.empty()   ? 3 * 4 * (long)packed_vertex_normals.size()   + PADDING : 0;
    long vertex_color_bytes    = !packed_vertex_colors.empty()    ? 4 * 4 * (long)packed_vertex_colors.size()    + PADDING : 0;
    long vertex_texcoord_bytes = !packed_vertex_texcoords.empty() ? 2 * 4 * (long)packed_vertex_texcoords.size() + PADDING : 0;

    packTopology();

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
    long num_bytes = vertex_position_bytes + vertex_normal_bytes + vertex_color_bytes + vertex_texcoord_bytes + PADDING;
#else
    long tri_bytes   =  !packed_tris.empty()   ?  4 * (long)packed_tris.size()   +  PADDING : 0;  // uint32
    long quad_bytes  =  !packed_quads.empty()  ?  4 * (long)packed_quads.size()  +  PADDING : 0;  // uint32
    long edge_bytes  =  !packed_edges.empty()  ?  4 * (long)packed_edges.size()  +  PADDING : 0;  // uint32

    long num_bytes = vertex_position_bytes
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

    if (var_area)
    {
      if (var_area->getCapacity() <= num_bytes || var_area->getCapacity() > (long)(1.5 * num_bytes))
      {
        render_system.destroyVARArea(var_area);

        std::string vararea_name = getNameStr() + " VAR area";
        var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, VARArea::Usage::WRITE_OCCASIONALLY, true);
        if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
      }
      // Else no need to reset var_area, we've done it above
    }
    else
    {
      std::string vararea_name = getNameStr() + " VAR area";
      var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, VARArea::Usage::WRITE_OCCASIONALLY, true);
      if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
    }

    if (!packed_vertex_positions.empty())
    {
      vertex_positions_var = var_area->createArray(vertex_position_bytes);
      if (!vertex_positions_var) throw Error(getNameStr() + ": Couldn't create vertices VAR");
      vertex_positions_var->updateVectors(0, (long)packed_vertex_positions.size(), &packed_vertex_positions[0]);
    }

    if (!packed_vertex_normals.empty())
    {
      vertex_normals_var = var_area->createArray(vertex_normal_bytes);
      if (!vertex_normals_var) throw Error(getNameStr() + ": Couldn't create normals VAR");
      vertex_normals_var->updateVectors(0, (long)packed_vertex_normals.size(), &packed_vertex_normals[0]);
    }

    if (!packed_vertex_colors.empty())
    {
      vertex_colors_var = var_area->createArray(vertex_color_bytes);
      if (!vertex_colors_var) throw Error(getNameStr() + ": Couldn't create colors VAR");
      vertex_colors_var->updateColors(0, (long)packed_vertex_colors.size(), &packed_vertex_colors[0]);
    }

    if (!packed_vertex_texcoords.empty())
    {
      vertex_texcoords_var = var_area->createArray(vertex_texcoord_bytes);
      if (!vertex_texcoords_var) throw Error(getNameStr() + ": Couldn't create texcoords VAR");
      vertex_texcoords_var->updateVectors(0, (long)packed_vertex_texcoords.size(), &packed_vertex_texcoords[0]);
    }

#ifndef THEA_GENERAL_MESH_NO_INDEX_ARRAY
    if (!packed_tris.empty())
    {
      tris_var = var_area->createArray(tri_bytes);
      if (!tris_var) throw Error(getNameStr() + ": Couldn't create triangle indices VAR");
      tris_var->updateIndices(0, (long)packed_tris.size(), &packed_tris[0]);
    }

    if (!packed_quads.empty())
    {
      quads_var = var_area->createArray(quad_bytes);
      if (!quads_var) throw Error(getNameStr() + ": Couldn't create quad indices VAR");
      quads_var->updateIndices(0, (long)packed_quads.size(), &packed_quads[0]);
    }

    if (!packed_edges.empty())
    {
      edges_var = var_area->createArray(edge_bytes);
      if (!edges_var) throw Error(getNameStr() + ": Couldn't create edge indices VAR");
      edges_var->updateIndices(0, (long)packed_edges.size(), &packed_edges[0]);
    }
#endif
  }
  else
  {
    if (!gpuBufferIsValid(BufferID::VERTEX_POSITION) && !vertices.empty())
    {
      packVertexPositions();
      vertex_positions_var->updateVectors(0, (long)packed_vertex_positions.size(), &packed_vertex_positions[0]);
    }

    if (!gpuBufferIsValid(BufferID::VERTEX_NORMAL) && !vertices.empty())
    {
      packVertexNormals();
      vertex_normals_var->updateVectors (0, (long)packed_vertex_normals.size(), &packed_vertex_normals[0]);
    }

    if (!gpuBufferIsValid(BufferID::VERTEX_COLOR) && hasVertexColors())
    {
      packVertexColors<Vertex>();
      vertex_colors_var->updateColors(0, (long)packed_vertex_colors.size(), &packed_vertex_colors[0]);
    }

    if (!gpuBufferIsValid(BufferID::VERTEX_TEXCOORD) && hasVertexTexCoords())
    {
      packVertexTexCoords<Vertex>();
      vertex_texcoords_var->updateVectors(0, (long)packed_vertex_texcoords.size(), &packed_vertex_texcoords[0]);
    }
  }

  allGPUBuffersAreValid();
}

template <typename V, typename E, typename F, template <typename T> class A>
inline void
GeneralMesh<V, E, F, A>::drawBuffered(RenderSystem & render_system, RenderOptions const & options) const
{
  if (options.drawEdges() && !buffered_wireframe)
    throw Error(getNameStr() + ": Can't draw mesh edges with GPU-buffered wireframe disabled");

  const_cast<GeneralMesh *>(this)->uploadToGraphicsSystem(render_system);

  if (!vertex_positions_var) return;
  if (!options.drawFaces() && !options.drawEdges()) return;
  if (!options.drawFaces() && edges.size() <= 0) return;
  if (!options.drawEdges() && faces.size() <= 0) return;

  render_system.beginIndexedPrimitives();

    render_system.setVertexArray(vertex_positions_var);
    if (options.sendNormals() && vertex_normals_var)
      render_system.setNormalArray(vertex_normals_var);
    else
      render_system.setNormalArray(NULL);

    if (options.sendColors() && vertex_colors_var)
      render_system.setColorArray(vertex_colors_var);
    else
      render_system.setColorArray(NULL);

    if (options.sendTexCoords() && vertex_texcoords_var)
      render_system.setTexCoordArray(0, vertex_texcoords_var);
    else
      render_system.setTexCoordArray(0, NULL);

    if (options.drawFaces())
    {
      if (options.drawEdges())
      {
        render_system.pushShapeFlags();
        render_system.setPolygonOffset(true, 2);
      }

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
        if (num_tri_indices > 0)
          render_system.sendIndices(RenderSystem::Primitive::TRIANGLES, num_tri_indices, &packed_tris[0]);

        if (num_quad_indices > 0)
          render_system.sendIndices(RenderSystem::Primitive::QUADS, num_quad_indices, &packed_quads[0]);
#else
        if (num_tri_indices > 0)
        {
          render_system.setIndexArray(tris_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::TRIANGLES, 0, num_tri_indices);
        }

        if (num_quad_indices > 0)
        {
          render_system.setIndexArray(quads_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::QUADS, 0, num_quad_indices);
        }
#endif

        // Finish off with all larger polygons
        if (has_large_polys)
        {
          for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
            if (fi->numEdges() > 4)
            {
              render_system.beginPrimitive(RenderSystem::Primitive::POLYGON);
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

        render_system.setShader(NULL);
        render_system.setColorArray(NULL);
        render_system.setTexCoordArray(0, NULL);
        render_system.setNormalArray(NULL);
        render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)

#ifdef THEA_GENERAL_MESH_NO_INDEX_ARRAY
        if (!edges.empty())
          render_system.sendIndices(RenderSystem::Primitive::LINES, 2 * (long)edges.size(), &packed_edges[0]);
#else
        if (!edges.empty())
        {
          render_system.setIndexArray(edges_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::LINES, 0, 2 * (long)edges.size());
        }
#endif

      render_system.popColorFlags();
      render_system.popShader();
    }

  render_system.endIndexedPrimitives();
}

template <typename V, typename E, typename F, template <typename T> class A>
inline void
GeneralMesh<V, E, F, A>::drawImmediate(RenderSystem & render_system, RenderOptions const & options) const
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
    render_system.beginPrimitive(RenderSystem::Primitive::TRIANGLES);
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isTriangle()) drawFace(*fi, render_system, options);
    render_system.endPrimitive();

    // Now render all quads
    render_system.beginPrimitive(RenderSystem::Primitive::QUADS);
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if (fi->isQuad()) drawFace(*fi, render_system, options);
    render_system.endPrimitive();

    // Finish off with all larger polygons
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->numEdges() > 4)
      {
        render_system.beginPrimitive(RenderSystem::Primitive::POLYGON);
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

      render_system.setShader(NULL);
      render_system.setColor(options.edgeColor());  // set default edge color

      render_system.beginPrimitive(RenderSystem::Primitive::LINES);
        for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        {
          ei->attr().draw(render_system, options);
          ei->getEndpoint(0)->drawPosition(render_system, options);
          ei->getEndpoint(1)->drawPosition(render_system, options);
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }
}

} // namespace Graphics
} // namespace Thea

#endif
