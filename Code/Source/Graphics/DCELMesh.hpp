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

#ifndef __Thea_Graphics_DCELMesh_hpp__
#define __Thea_Graphics_DCELMesh_hpp__

// #define THEA_DCELMESH_VERBOSE

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Colors.hpp"
#include "../NamedObject.hpp"
#include "../Set.hpp"
#include "../UnorderedSet.hpp"
#include "DCELFace.hpp"
#include "DCELVertex.hpp"
#include "DCELHalfedge.hpp"
#include "GraphicsAttributes.hpp"
#include "IncrementalDCELMeshBuilder.hpp"
#include "DefaultMeshCodecs.hpp"
#include "DrawableObject.hpp"
#include <cmath>
#include <limits>

#ifdef THEA_DCELMESH_VERBOSE
#  include "../UnorderedMap.hpp"
#endif

namespace Thea {
namespace Graphics {

/**
 * Mesh based on a doubly-connected edge list (or halfedge data structure).
 *
 * Adapted from: DCELMesh class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 *
 * @todo Automatically invalidate appropriate GPU buffers on all modifications.
 */
template < typename VertexAttribute    =  Graphics::NullAttribute,
           typename HalfedgeAttribute  =  Graphics::NullAttribute,
           typename FaceAttribute      =  Graphics::NullAttribute >
class /* THEA_API */ DCELMesh : public virtual NamedObject, public DrawableObject
{
  public:
    THEA_DEF_POINTER_TYPES(DCELMesh, shared_ptr, weak_ptr)

    /** Mesh type tag. */
    struct DCEL_MESH_TAG {};

    typedef DCELVertex  <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Vertex;    ///< Vertex of the mesh (3-vector).
    typedef DCELHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute>  Halfedge;  ///< Halfedge of the mesh.
    typedef DCELFace    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Face;      ///< Face of the mesh.

  private:
    struct HalfedgeComparator
    {
      bool operator()(Halfedge const * e0, Halfedge const * e1) const
      {
        return e0 ? (e1 ? e0->index < e1->index : false) : (bool)e1;
      }
    };

    typedef TheaUnorderedSet<Vertex *>               VertexSet;
    typedef TheaSet<Halfedge *, HalfedgeComparator>  HalfedgeSet;  // store halfedges in sequence so twins are consecutive
    typedef TheaUnorderedSet<Face *>                 FaceSet;

    /** Iterate over edges (every other halfedge). */
    template <typename BaseIterT>
    struct EdgeIterTmpl : public BaseIterT
    {
      /** Default constructor. */
      EdgeIterTmpl() {}

      /** General copy constructor. */
      template <typename T> EdgeIterTmpl(T const & src) : BaseIterT(src) {}

      /** Assignment. Use with caution: \a src must have even index. */
      EdgeIterTmpl & operator=(BaseIterT const & src) { BaseIterT::operator=(src); }

      /** Pre-increment. */
      EdgeIterTmpl & operator++()
      {
        BaseIterT::operator++();
        BaseIterT::operator++();
        return *this;
      }

      /** Pre-decrement. */
      EdgeIterTmpl & operator--()
      {
        BaseIterT::operator--();
        BaseIterT::operator--();
        return *this;
      }

      /** Post-increment. */
      EdgeIterTmpl operator++(int)
      {
        EdgeIterTmpl ret = *this;
        BaseIterT::operator++();
        BaseIterT::operator++();
        return ret;
      }

      /** Post-decrement. */
      EdgeIterTmpl operator--(int)
      {
        EdgeIterTmpl ret = *this;
        BaseIterT::operator--();
        BaseIterT::operator--();
        return ret;
      }
    };

  public:
    typedef typename VertexSet::iterator          VertexIterator;         ///< Iterator over vertices.
    typedef typename VertexSet::const_iterator    VertexConstIterator;    ///< Const iterator over vertices.
    typedef typename HalfedgeSet::iterator        HalfedgeIterator;       ///< Iterator over halfedges.
    typedef typename HalfedgeSet::const_iterator  HalfedgeConstIterator;  ///< Const iterator over halfedges.
    typedef typename FaceSet::iterator            FaceIterator;           ///< Iterator over faces.
    typedef typename FaceSet::const_iterator      FaceConstIterator;      ///< Const iterator over faces.

    /** Iterator over edges (alternate halfedges starting from the first). */
    typedef EdgeIterTmpl<HalfedgeIterator> EdgeIterator;

    /** Const iterator over edges (alternate halfedges starting from the first). */
    typedef EdgeIterTmpl<HalfedgeConstIterator> EdgeConstIterator;

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
    DCELMesh(std::string const & name = "AnonymousMesh")
    : NamedObject(name),
      next_halfedge_index(0),
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
     * Copy constructor. Creates a deep copy of the mesh (including copies of the attributes). <b>Currently not
     * implemented.</b>
     */
    DCELMesh(DCELMesh const & src)
    : NamedObject(src), bounds(src.bounds)
    {
      throw Error("DCELMesh: Copy constructor not currently implemented");
    }

    ~DCELMesh() { clear(); }

    /** Get an iterator pointing to the first vertex. */
    VertexConstIterator verticesBegin() const { return vertices.begin(); }

    /** Get an iterator pointing to the first vertex. */
    VertexIterator verticesBegin() { return vertices.begin(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexConstIterator verticesEnd() const { return vertices.end(); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexIterator verticesEnd() { return vertices.end(); }

    /** Get an iterator pointing to the first halfedge. */
    HalfedgeConstIterator halfedgesBegin() const { return halfedges.begin(); }

    /** Get an iterator pointing to the first halfedge. */
    HalfedgeIterator halfedgesBegin() { return halfedges.begin(); }

    /** Get an iterator pointing to the position beyond the last halfedge. */
    HalfedgeConstIterator halfedgesEnd() const { return halfedges.end(); }

    /** Get an iterator pointing to the position beyond the last halfedge. */
    HalfedgeIterator halfedgesEnd() { return halfedges.end(); }

    /** Get an iterator pointing to the first edge. */
    EdgeConstIterator edgesBegin() const { return halfedgesBegin(); }

    /** Get an iterator pointing to the first edge. */
    EdgeIterator edgesBegin() { return halfedgesBegin(); }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeConstIterator edgesEnd() const
    {
      debugAssertM(numHalfedges() % 2 == 0, getNameStr() + ": Number of halfedges is odd");
      return halfedgesEnd();
    }

    /** Get an iterator pointing to the position beyond the last edge. */
    EdgeIterator edgesEnd()
    {
      debugAssertM(numHalfedges() % 2 == 0, getNameStr() + ": Number of halfedges is odd");
      return halfedgesEnd();
    }

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
      for (VertexIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
        delete *vi;

      for (HalfedgeIterator ei = halfedges.begin(); ei != halfedges.end(); ++ei)
        delete *ei;

      for (FaceIterator fi = faces.begin(); fi != faces.end(); ++fi)
        delete *fi;

      vertices.clear();
      halfedges.clear();
      faces.clear();

      next_halfedge_index = 0;
      bounds = AxisAlignedBox3();

      invalidateGPUBuffers();
    }

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const { return vertices.empty() && faces.empty() && halfedges.empty(); }

    /** Get the number of vertices. */
    long numVertices() const { return (long)vertices.size(); }

    /** Get the number of halfedges. */
    long numHalfedges() const { return (long)halfedges.size(); }

    /** Get the number of edges. */
    long numEdges() const
    {
      debugAssertM(halfedges.size() % 2 == 0, getNameStr() + ": Number of halfedges is odd");
      return (long)(halfedges.size() / 2);
    }

    /** Get the number of faces. */
    long numFaces() const { return (long)faces.size(); }

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
     * Add a vertex to the mesh, with an optional precomputed normal. The index of the vertex will be the value of
     * getNumVertices() before the addition. Automatically updates the bounding box of the mesh, so you don't need to call
     * updateBounds() after calling this function.
     *
     * @return A pointer to the new vertex on success, null on failure.
     */
    Vertex * addVertex(Vector3 const & point, Vector3 const * normal = NULL)
    {
      if (vertices.empty()) bounds.set(point, point);
      else                  bounds.merge(point);

      Vertex * vertex = normal ? new Vertex(point, *normal) : new Vertex(point);
      vertices.insert(vertex);

#ifdef THEA_DCELMESH_VERBOSE
      vertex_indices[vertex] = vertices.size() - 1;
      THEA_CONSOLE << "Added vertex " << vertices.size() - 1 << " at " << vertex->getPosition();
#endif

      invalidateGPUBuffers();

      return vertex;
    }

    /**
     * Add a face to the mesh, specified by the sequence of vertices obtained by dereferencing [vbegin, vend).
     * VertexInputIterator must dereference to a pointer to a Vertex. Unless the mesh is already in an inconsistent state,
     * failure to add the face will not affect the mesh.
     *
     * @return A pointer to the new created face, or null on error.
     */
    template <typename VertexInputIterator>
    Face * addFace(VertexInputIterator vbegin, VertexInputIterator vend)
    {
      if (vbegin == vend)
      {
        THEA_WARNING << getName() << ": Skipping face -- it has no vertices";
        return NULL;
      }

#ifdef THEA_DCELMESH_VERBOSE
      std::cout << "Adding face:";
#endif

      // Read the vertex pointers into an internal array
      if (face_vertices.size() != 256) face_vertices.resize(256);  // default size, should ensure not too many resizes
      array_size_t num_verts = 0;
      for (VertexInputIterator vi = vbegin; vi != vend; ++vi, ++num_verts)
      {
        debugAssertM(*vi, getNameStr() + ": Null vertex pointer specified for new face");

        if (num_verts >= face_vertices.size())
          face_vertices.resize(2 * face_vertices.size() + 1);

        face_vertices[num_verts] = *vi;

#ifdef THEA_DCELMESH_VERBOSE
        std::cout << ' ' << vertex_indices[*vi];
#endif
      }

#ifdef THEA_DCELMESH_VERBOSE
      std::cout << std::endl;
#endif

      if (num_verts < 3)
      {
        THEA_WARNING << getName() << ": Skipping face -- too few vertices (" << num_verts << ')';
        return NULL;
      }

      // Compute the face normal, assume it is consistent across the face
      Vector3 e1 = face_vertices[0]->getPosition() - face_vertices[1]->getPosition();
      Vector3 e2 = face_vertices[2]->getPosition() - face_vertices[1]->getPosition();
      Vector3 normal = e2.cross(e1).unit();  // counter-clockwise

      Face * face = addFace(num_verts, &face_vertices[0], normal);

      invalidateGPUBuffers();

      return face;
    }

    /**
     * Split an edge along its length, in the ratio given by \a frac.
     *
     * @return The new vertex created by the split operation, or null on error.
     */
    Vertex * splitEdge(Halfedge * edge, Real frac)
    {
      alwaysAssertM(frac >= 0 && frac <= 1, getNameStr() + ": Edge split fraction should be between 0 and 1")
      if (edge)
      {
        THEA_ERROR << getName() << "Can't split null edge";
        return NULL;
      }

      Vector3 p = (1 - frac) * edge->getOrigin()->getPosition() + frac * edge->getEnd()->getPosition();
      Vector3 n = ((1 - frac) * edge->getOrigin()->getNormal() + frac * edge->getEnd()->getNormal()).unit();
      Vertex * new_vx = addVertex(p, &n);
      if (!new_vx)
        return NULL;

      if (!splitEdge(edge, new_vx))  // should generally never happen
      {
        removeVertex(new_vx);
        return NULL;
      }

      return new_vx;
    }

    /**
     * Split an edge along its length at a given position \a pos.
     *
     * @return The new vertex created by the split operation, or null on error.
     */
    Vertex * splitEdge(Halfedge * edge, Vector3 const & pos)
    {
      Real s = (pos - edge->getOrigin()->getPosition()).length();
      Real t = (pos - edge->getEnd()->getPosition()).length();
      Vector3 n = (t * edge->getOrigin()->getNormal() + s * edge->getEnd()->getNormal()).unit();
      Vertex * new_vx = addVertex(pos, &n);
      if (!new_vx)
        return NULL;

      if (!splitEdge(edge, new_vx))  // should generally never happen
      {
        removeVertex(new_vx);
        return NULL;
      }

      return new_vx;
    }

    /**
     * Stitch two boundary loops together, closing the gap between them. Each loop is specified by a boundary edge on it.
     * Starting from these edges, pairs of edges are successively glued together. The loops must have exactly the same number of
     * edges. This function steps <em>forward</em> along the first loop (in the direction of the next pointers) and
     * <em>backward</em> in the second loop (against the direction of the next pointers).
     *
     * @param loop0_start Starting edge of the first loop. Must be a boundary edge (no linked face). This is glued to
     *   \a loop1_start.
     * @param loop1_start Starting edge of the second loop. Must be a boundary edge (no linked face). This is glued to
     *   \a loop0_start.
     *
     * @return True if the loops were successfully glued, else false.
     */
    bool stitchBoundaryLoops(Halfedge * loop0_start, Halfedge * loop1_start)
    {
      // In the first pass, find the length of the loops. This makes the stopping condition easier to evaluate since we don't
      // need to ensure that we're not examining a deleted halfedge, and also allows us to check for some errors before we make
      // any changes.
      Halfedge * e0 = loop0_start, * e1 = loop1_start;
      long n = 0;
      do
      {
        if (!e0->isBoundary() || !e1->isBoundary())
        {
          THEA_ERROR << getName() << ": Can't stitch loops with non-boundary edges";
          return false;
        }

        if (++n > numEdges())
        {
          THEA_ERROR << getName() << ": Loop is not closed";
          return false;
        }

        e0 = e0->next();
        e1 = e1->next();

      } while (e0 != loop0_start);

      if (e1 != loop1_start)
      {
        THEA_ERROR << getName() << ": Loops are not of same length";
        return false;
      }

      // Now do the actual stitching
      e0 = loop0_start;  // we should be here already, but let's play safe
      e1 = loop1_start;
      Halfedge * next0 = NULL, * next1 = NULL;
      for (long i = 0; i < n; ++i, e0 = next0, e1 = next1)
      {
        if (i < n - 1)
        {
          // Store in advance where we're going to go next, since e0 and e1 are going to be deleted in this iteration
          next0 = e0->next();
          next1 = e1->origin->leaving;  // we're stepping backwards around the second loop
          do
          {
            if (next1->twin_he->next() == e1)
              break;

            next1 = next1->nextAroundOrigin();

          } while (next1 != e1->origin->leaving);

          next1 = next1->twin_he;
          if (next1->next() != e1)
          {
            THEA_ERROR << getName() << ": Predecessor of halfedge in second loop not found";
            return false;
          }
        }

        // Shift all halfedges originating at e1's origin to the new origin at the tip of e0
        Vertex * old_origin = e1->origin, * new_origin = e0->twin_he->origin;
        Halfedge * he = old_origin->leaving;  // this is guaranteed to not be a deleted loop edge, since we're stepping
        do                                    // backwards in the second loop
        {
          he->origin = new_origin;
          he = he->nextAroundOrigin();  // this is also guaranteed to be a valid operation, the relevant next and twin pointers
                                        // are ok
        } while (he != old_origin->leaving);

        // Delete the old origin
        removeVertex(old_origin);

        // Seal the edge by deleting the loop-side halfedges, and making the face-side halfedges twins
        Halfedge * t0 = e0->twin_he;
        Halfedge * t1 = e1->twin_he;
        long index = e0->index;

        removeHalfedge(e0);
        removeHalfedge(e1);

        // Update the index of the face-side halfedge in the second loop, so it is stored right after its new twin
        halfedges.erase(t1);
        t1->index = index;
        halfedges.insert(t1);

        // Update the twin pointers to complete the pairing
        t0->twin_he = t1;
        t1->twin_he = t0;
      }

      return true;
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
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
        bounds.merge((*vi)->getPosition());
    }

    AxisAlignedBox3 const & getBounds() const { return bounds; }

  private:
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
     * Add a face to the mesh, specified by a sequence of boundary vertices.
     *
     * @return A pointer to the newly created face.
     */
    Face * addFace(array_size_t num_verts, Vertex ** verts, Vector3 const & normal)
    {
      // Try to locate a boundary edge that is already in the mesh
      array_size_t origin = 0;
      Halfedge * e = findExistingEdge(num_verts, verts, origin);

      Face * face = NULL;
      face = addFace(num_verts, verts, e, origin, normal);
      if (!face)
      {
        THEA_WARNING << "Adding isolated face";
        TheaArray<Vertex *> iso_verts;
        array_size_t existing0 = 0, existing1 = 0;
        if (e && e->isBoundaryEdge())
        {
          existing0 = origin;
          existing1 = (origin + 1) % num_verts;
        }

        if (!addIsolatedVertices(num_verts, verts, existing0, existing1, iso_verts))
          return NULL;

        if (e && e->isBoundaryEdge())
          face = addFace(num_verts, const_cast<Vertex **>(&iso_verts[0]), e, origin, normal);
        else
          face = addFace(num_verts, const_cast<Vertex **>(&iso_verts[0]), NULL, 0, false, normal);
      }

      return face;
    }

    /** Utility function for addFace(array_size_t, Vertex **, Vector3 const &). */
    Face * addFace(array_size_t num_verts, Vertex ** verts, Halfedge * first, array_size_t origin, Vector3 const & normal)
    {
      if (first)
      {
        if (first->isBoundary())
        {
          // All good... keep adding starting here
          return addFace(num_verts, verts, first, origin, false, normal);
        }
        else
        {
          if (first->twin()->isBoundary())
          {
            // Try adding in reverse
            THEA_DEBUG << getName() << ": Trying to add face by reversing the order of vertices";
            return addFace(num_verts, verts, first->twin(), (origin + 1) % num_verts, true, -normal);  // normal got flipped
          }
          else
          {
            THEA_WARNING << getName() << ": Face has edge already adjoining two faces";
            return NULL;
          }
        }
      }
      else
        return addFace(num_verts, verts, NULL, 0, false, normal);
    }

    /**
     * Add a sequence of isolated vertices, optionally including one existing edge (if existing0 != existing1), and return
     * pointers to the new vertices.
     */
    bool addIsolatedVertices(array_size_t num_verts, Vertex ** verts, array_size_t existing0, array_size_t existing1,
                             TheaArray<Vertex *> & iso_verts)
    {
      bool has_existing_edge = (existing0 != existing1);
      for (array_size_t i = 0; i < num_verts; ++i)
      {
        if (has_existing_edge && (i == existing0 || i == existing1))
          iso_verts.push_back(verts[i]);
        else
        {
          Vertex * v = addVertex(verts[i]->getPosition());
          if (!v)
          {
            THEA_WARNING << "Could not add isolated vertex";
            return false;
          }

          v->setAttr(verts[i]->attr());
          iso_verts.push_back(v);
        }
      }

      return true;
    }

    /**
     * Add a face to the mesh, specified by a sequence of boundary vertices. The boundary edges are added starting with a
     * specified first one. Optionally, the boundary is traversed in reverse order.
     *
     * @param num_verts Number of boundary vertices.
     * @param verts Boundary vertices, sequentially.
     * @param first The first existing halfedge of the face. If non-null, it must be on the mesh boundary.
     * @param origin The index (in indices) of the starting vertex of the first halfedge.
     * @param reverse Insert edges in reverse order?
     * @param normal The precomputed unit face normal.
     *
     * @return A pointer to the newly created face.
     */
    Face * addFace(array_size_t num_verts, Vertex ** verts, Halfedge * first, array_size_t origin, bool reverse,
                   Vector3 const & normal)
    {
      // For each vertex, check that the halfedge emanating from it exists as a boundary halfedge, or can be successfully added.
      // Store either the existing edge to the next vertex, or the preceding boundary halfedge _into_ the vertex.
      TheaArray<Halfedge *> edges(num_verts);
      Halfedge * prev = first, * e;
      array_size_t start_index = reverse ? (origin > 0 ? origin - 1 : num_verts - 1)
                                         : (origin < num_verts - 1 ? origin + 1 : 0);
      array_size_t i = start_index, next;
      do
      {
        // See if there is an existing edge
        next = reverse ? (i > 0 ? i - 1 : num_verts - 1) : (i < num_verts - 1 ? i + 1 : 0);
        e = verts[i]->getEdgeTo(verts[next]);

        if (e)  // edge exists, needs to be boundary edge
        {
          if (!e->isBoundary())
          {
            THEA_WARNING << getName() << ": Can't stitch a face to a halfedge that already adjoins a face";
            return NULL;
          }

          if (prev && prev->next() != e)
          {
            THEA_WARNING << getName() << ": Face breaks halfedge order at vertex";
            return NULL;
          }

          edges[i] = e;
          prev = e;
        }
        else  // edge does not exist
        {
          if (prev)
          {
            // Just to be sure, check that prev->next() is a boundary edge
            debugAssertM(prev->isBoundary(),
                         getNameStr() + ": Previous halfedge of new face is not currently on the mesh boundary");
            debugAssertM(prev->next() && prev->next()->isBoundary(),
                         getNameStr()
                       + ": Mesh has internal inconsistency (boundary halfedge not followed by boundary halfedge)");

            // Twin of prev is predecessor of (i, next) around vertex i.
            edges[i] = prev;
          }
          else
          {
            // Check that the edge can be added, i.e. the vertex is on the boundary, and find its predecessor
            Halfedge * e = findPrevAroundVertex(verts[i], verts[next], normal);
            if (e)
            {
              if (e->twin()->isBoundary())
                edges[i] = e->twin();
              else
              {
                THEA_WARNING << getName() << ": Can't stitch a face to a halfedge that already adjoins a face";
                return NULL;
              }
            }
            else
              edges[i] = NULL;
          }

          prev = NULL;  // no current edge between i and next
        }

        i = next;

      } while (i != start_index);

      // Now create the new face
      Face * face = new Face;
      face->num_edges = (int)num_verts;
      face->setNormal(normal);

      // Now actually create and add the edges. In this loop we directly access private members (ha ha) of the edges to bypass
      // any consistency checks until we've finished adding the face.
      Halfedge * new_e, * new_twin, * next_e, * last = NULL;
      Vertex * vi, * vnext;
      long index0, index1;
      bool added_canonical_edge_to_face = false;

#ifdef THEA_DCELMESH_VERBOSE
      bool last_existed = false, this_exists = false;
#endif

      i = start_index;
      vi = verts[i];

      do
      {
        next = reverse ? (i > 0 ? i - 1 : num_verts - 1) : (i < num_verts - 1 ? i + 1 : 0);
        vnext = verts[next];

        e = edges[i];
        if (e && e->origin == vi)  // edge exists, no problem
        {
#ifdef THEA_DCELMESH_VERBOSE
          this_exists = true;
#endif

          debugAssertM(e->getOrigin() == vi && e->getEnd() == vnext, getNameStr() + ": Edge has wrong endpoints");
          debugAssertM(e->isBoundary(),
                       getNameStr() + ": Can't stitch a face to a halfedge that already adjoins a face");
          new_e = e;

          // This edge is necessarily the predecessor of the edge we will add at the next vertex
          next_e = edges[next];

          debugAssertM(next_e, getNameStr()
                             + ": Vertex has an existing edge, yet it was not found in search for predecessor");

          if (next_e->origin != vnext)  // no existing next edge
            edges[next] = e;
        }
        else  // edge doesn't exist, predecessor is stored
        {
#ifdef THEA_DCELMESH_VERBOSE
          this_exists = false;
#endif

          // Create a new pair of halfedges
          nextHalfedgeIndices(index0, index1);

          new_e = new Halfedge(index0);
          halfedges.insert(new_e);

          new_twin = new Halfedge(index1);
          halfedges.insert(new_twin);

          // Each is the twin of the other
          new_e->twin_he    = new_twin;
          new_twin->twin_he = new_e;

          // Set the halfedge origins
          new_e->origin = vi;
          new_twin->origin = vnext;

          // The new edge is (duh) the new edge connecting this vertex to the next
          edges[i] = new_e;

          // Insert the edge at this vertex
          if (e)
          {
            new_twin->next_he = e->next_he;
            e->next_he = new_e;
          }
          else
          {
            new_twin->next_he = new_e;  // vertex has degree 1 for now (will be fixed when the boundary loops round)
            vi->leaving = new_e;
          }

          // Insert the edge at the next vertex (the new edge is the new predecessor at the next vertex)
          next_e = edges[next];
          if (next_e)
          {
            if (next_e->origin == vnext)  // next edge exists
            {
              new_e->next_he = next_e;
              prev = next_e->prevAroundOrigin();  // this search can be avoided if the first edge was newly created, and we
                                                  // have looped back to it (since we can reuse the predecessor that was
                                                  // initially stored); however, in the interest of clarity we'll omit this
                                                  // optimization
              prev->twin_he->next_he = new_twin;
            }
            else  // predecessor is directly stored (note that it points _into_ the next vertex, so it corresponds to prev->twin
                  // of the if block)
            {
              new_e->next_he  = next_e->next_he;
              next_e->next_he = new_twin;

              edges[next] = new_e;
            }
          }
          else  // isolated next vertex
          {
            new_e->next_he = new_twin;  // next vertex has degree 1 for now (will be fixed by next edge)
            vnext->leaving = new_twin;

            edges[next] = new_e;
          }
        }

        // Update the face pointer of the halfedge (the twin, if newly created, becomes a boundary halfedge)
        new_e->face = face;
        if (!added_canonical_edge_to_face)
        {
          face->halfedge = new_e;
          added_canonical_edge_to_face = true;
        }

        debugAssertM(!last || last->next_he == new_e,
                     getNameStr() + ": Next pointers on face boundary not consistent");
        last = new_e;

#ifdef THEA_DCELMESH_VERBOSE
        THEA_CONSOLE << "Last existed = " << last_existed << ", this exists = " << this_exists;
        last_existed = this_exists;
#endif

        // Update the normal at the vertex
        vi->addFaceNormal(normal);  // weight by face area?

        i = next;
        vi = vnext;

      } while (i != start_index);

      faces.insert(face);
      return face;
    }

    /**
     * Get indices for a pair of new twin halfedges. If the indices are getting too large, this function reindexes the set of
     * halfedges.
     */
    void nextHalfedgeIndices(long & index0, long & index1)
    {
      static long const THRESHOLD = std::numeric_limits<long>::max() - 4;
      if (next_halfedge_index > THRESHOLD)
      {
        if ((long)halfedges.size() > THRESHOLD)  // too many halfedges, we can't do anything!
          throw Error(getNameStr() + ": Too many halfedges, can't assign indices!");

        next_halfedge_index = 0;
        for (HalfedgeIterator ei = halfedges.begin(); ei != halfedges.end(); ++ei)
          (*ei)->index = next_halfedge_index++;
      }

      index0 = next_halfedge_index++;
      index1 = next_halfedge_index++;
    }

    /**
     * Find an edge that is already in the mesh, from a loop of edges specified by the sequence of their vertices. This function
     * will preferentially return a boundary edge if it finds one.
     *
     * @return A pointer to an existing half-edge, if found, else null. The index of the originating vertex of this half-edge
     *   (respective to the input array of vertices) is stored in origin_index.
     */
    Halfedge * findExistingEdge(array_size_t num_verts, Vertex ** verts, array_size_t & origin_index)
    {
      Halfedge * e  = NULL;
      array_size_t last = num_verts - 1;
      for (array_size_t i = 0; i < num_verts; ++i)
      {
        e = verts[last]->getEdgeTo(verts[i]);
        if (e)
        {
          origin_index = last;

          if (e->isBoundary() || e->twin()->isBoundary())
            break;
        }

        last = i;
      }

      return e;
    }

    /**
     * Get the angle going counter-clockwise from \a u to \a v w.r.t. a given up direction, in the range [0, 2 pi]. Assumes \a u
     * is orthogonal to \a unit_up.
     */
    static float ccwAngle(Vector3 const & u, Vector3 const & v, Vector3 const & unit_up)
    {
      Vector3 v_proj = v - (v.dot(unit_up) * unit_up);
      float s = u.cross(v_proj).fastLength();
      float c = u.dot(v_proj);
      float ang = std::atan2(s, c);

      return ang < 0 ? 2 * M_PI - ang : ang;
    }

    /**
     * Find the predecessor of a yet-to-be-added halfedge around a vertex. The predecessor, if found, will be the twin of a
     * boundary edge.
     *
     * @param u The starting vertex of the halfedge.
     * @param v The end vertex of the halfedge.
     * @param normal The precomputed unit normal of the face adjoining the halfedge from u to v.
     *
     * @return The predecessor, if found, else null.
     */
    Halfedge * findPrevAroundVertex(Vertex * u, Vertex * v, Vector3 const & normal)
    {
      Halfedge * first = u->getHalfedge();
      if (!first)
        return NULL;

      // Compute the new normal, which defines the plane of projection
      Vector3 new_normal = u->estimateUpdatedNormal(normal);
      Halfedge * e = first, * best = first;
      float ang, best_ang = 100;  // anything > 2 pi
      Vector3 edge = v->getPosition() - u->getPosition();
      edge = (edge - (edge.dot(new_normal) * new_normal));
      do
      {
        if (e->twin()->isBoundary())
        {
          Vector3 test_edge = (e->getEnd()->getPosition() - u->getPosition());
          ang = ccwAngle(edge, test_edge, new_normal);
          if (ang < best_ang)
          {
            best_ang = ang;
            best = e;
          }
        }

        e = e->nextAroundOrigin();  // the next halfedge around u, counter-clockwise when looking from the "outside"

      } while (e != first);

      return best;
    }

    /** Delete a vertex from the mesh. No pointers are updated, the vertex is just deleted from storage. */
    void removeVertex(Vertex * vertex)
    {
      vertices.erase(vertex);
      delete vertex;
    }

    /** Delete a halfedge from the mesh. No pointers are updated, the halfedge is just deleted from storage. */
    void removeHalfedge(Halfedge * halfedge)
    {
      halfedges.erase(halfedge);
      delete halfedge;
    }

    /** Delete a face from the mesh. No pointers are updated, the face is just deleted from storage. */
    void removeFace(Face * face)
    {
      faces.erase(face);
      delete face;
    }

    /** Split an edge along its length at a given vertex location. The vertex is assumed to have been newly added. */
    bool splitEdge(Halfedge * edge, Vertex * vertex)
    {
      Halfedge * twin = edge->twin();

      long twin_index = twin->index;
      long index0 = -1, index1 = -1;
      nextHalfedgeIndices(index0, index1);

      // Update the twin's index
      halfedges.erase(twin);
      twin->index = index0;
      halfedges.insert(twin);

      // Create a new pair of halfedges
      Halfedge * e0 = new Halfedge(twin_index);
      Halfedge * e1 = new Halfedge(index1);
      e0->origin = e1->origin = vertex;
      e0->face = twin->face;
      e1->face = edge->face;

      halfedges.insert(e0);
      halfedges.insert(e1);

      // Update twin pairings
      edge->twin_he = e0;
      e0->twin_he = edge;

      twin->twin_he = e1;
      e1->twin_he = twin;

      // Update next pointers
      e1->next_he = edge->next_he;
      edge->next_he = e1;

      e0->next_he = twin->next_he;
      twin->next_he = e0;

      // Mark one of the new halfedges (doesn't matter which) as leaving the new vertex
      vertex->leaving = e0;

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

      Halfedge const * first = face.getHalfedge();
      Halfedge const * e = first;
      do
      {
        Vertex const * v = e->getOrigin();

        if (options.useVertexNormals() && options.sendNormals())
          v->drawNormal(render_system, options);

        if (options.useVertexData())  // vertex attributes (per-vertex color, texcoord etc)
          v->attr().draw(render_system, options);

        v->drawPosition(render_system, options);  // finally send the position of the vertex

        e = e->next();

      } while (e != first);
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
        packed_vertex_positions[i] = (*vi)->getPosition();
    }

    /** Pack vertex positions densely in an array. */
    void packVertexNormals()
    {
      packed_vertex_normals.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_normals[i] = (*vi)->getNormal();
    }

    /** Pack vertex colors densely in an array. */
    template <typename VertexT>
    void packVertexColors(typename boost::enable_if< HasColor<VertexT> >::type * dummy = NULL)
    {
      packed_vertex_colors.resize(vertices.size());
      array_size_t i = 0;
      for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_colors[i] = ColorRGBA((*vi)->attr().getColor());
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
        packed_vertex_texcoords[i] = (*vi)->attr().getTexCoord();
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
        (*vi)->setPackingIndex(index++);

      has_large_polys = false;
      for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
      {
        Face * face = *fi;
        if (face->isTriangle())
        {
          Halfedge const * he = face->getHalfedge();
          packed_tris.push_back(he->getOrigin()->getPackingIndex()); he = he->next();
          packed_tris.push_back(he->getOrigin()->getPackingIndex()); he = he->next();
          packed_tris.push_back(he->getOrigin()->getPackingIndex());
        }
        else if (face->isQuad())
        {
          Halfedge const * he = face->getHalfedge();
          packed_tris.push_back(he->getOrigin()->getPackingIndex()); he = he->next();
          packed_tris.push_back(he->getOrigin()->getPackingIndex()); he = he->next();
          packed_tris.push_back(he->getOrigin()->getPackingIndex()); he = he->next();
          packed_tris.push_back(he->getOrigin()->getPackingIndex());
        }
        else
        {
          has_large_polys = true;
          THEA_DEBUG << getName() << ": Mesh has polygons with 5 or more vertices. These will not be drawn with GPU buffers.";
        }
      }

      if (buffered_wireframe)
      {
        packed_edges.resize(halfedges.size());  // i.e. 2 * numEdges()
        array_size_t i = 0;
        for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei, i += 2)
        {
          packed_edges[i    ] = (*ei)->getOrigin()->getPackingIndex();
          packed_edges[i + 1] = (*ei)->getEnd()->getPackingIndex();
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

    FaceSet          faces;
    VertexSet        vertices;
    HalfedgeSet      halfedges;
    long             next_halfedge_index;
    AxisAlignedBox3  bounds;

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

#ifdef THEA_DCELMESH_VERBOSE
    typedef TheaUnorderedMap<Vertex *, array_size_t> VertexIndexMap;
    VertexIndexMap vertex_indices;
#endif

    mutable TheaArray<Vertex *> face_vertices;  // internal cache for vertex pointers for a face
};

template <typename V, typename E, typename F>
inline void
DCELMesh<V, E, F>::uploadToGraphicsSystem(RenderSystem & render_system)
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

    if (vertices.empty() || (faces.empty() && halfedges.empty()))
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

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
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

#ifndef THEA_DCEL_MESH_NO_INDEX_ARRAY
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

template <typename V, typename E, typename F>
inline void
DCELMesh<V, E, F>::drawBuffered(RenderSystem & render_system, RenderOptions const & options) const
{
  if (options.drawEdges() && !buffered_wireframe)
    throw Error(getNameStr() + ": Can't draw mesh edges with GPU-buffered wireframe disabled");

  const_cast<DCELMesh *>(this)->uploadToGraphicsSystem(render_system);

  if (!vertex_positions_var) return;
  if (!options.drawFaces() && !options.drawEdges()) return;
  if (!options.drawFaces() && halfedges.size() <= 0) return;
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

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
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
            if ((*fi)->numEdges() > 4)
            {
              render_system.beginPrimitive(RenderSystem::Primitive::POLYGON);
                drawFace(**fi, render_system, options);
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

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
        if (!edges.empty())
          render_system.sendIndices(RenderSystem::Primitive::LINES, (long)halfedges.size(), &packed_edges[0]);
#else
        if (!halfedges.empty())
        {
          render_system.setIndexArray(edges_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::LINES, 0, (long)halfedges.size());
        }
#endif

      render_system.popColorFlags();
      render_system.popShader();
    }

  render_system.endIndexedPrimitives();
}

template <typename V, typename E, typename F>
inline void
DCELMesh<V, E, F>::drawImmediate(RenderSystem & render_system, RenderOptions const & options) const
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
        if ((*fi)->isTriangle()) drawFace(**fi, render_system, options);
    render_system.endPrimitive();

    // Now render all quads
    render_system.beginPrimitive(RenderSystem::Primitive::QUADS);
      for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
        if ((*fi)->isQuad()) drawFace(**fi, render_system, options);
    render_system.endPrimitive();

    // Finish off with all larger polygons
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if ((*fi)->numEdges() > 4)
      {
        render_system.beginPrimitive(RenderSystem::Primitive::POLYGON);
          drawFace(**fi, render_system, options);
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
          Halfedge * edge = *ei;
          edge->attr().draw(render_system, options);
          render_system.sendVertex(edge->getOrigin()->getPosition());
          render_system.sendVertex(edge->getEnd()->getPosition());
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }
}

} // namespace Graphics
} // namespace Thea

#endif
