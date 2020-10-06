//============================================================================ //
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

#ifndef __Thea_Graphics_DcelMesh_hpp__
#define __Thea_Graphics_DcelMesh_hpp__

// #define THEA_DCELMESH_VERBOSE

#include "../Common.hpp"
#include "../Algorithms/IteratorModifiers.hpp"
#include "../Array.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Colors.hpp"
#include "../Math.hpp"
#include "../MatrixWrapper.hpp"
#include "../NamedObject.hpp"
#include "../Set.hpp"
#include "../UnorderedSet.hpp"
#include "IMesh.hpp"
#include "DcelFace.hpp"
#include "DcelVertex.hpp"
#include "DcelHalfedge.hpp"
#include "DefaultMeshCodecs.hpp"
#include "GraphicsAttributes.hpp"
#include "IncrementalDcelMeshBuilder.hpp"
#include <limits>

#ifdef THEA_DCELMESH_VERBOSE
#  include "../UnorderedMap.hpp"
#endif

namespace Thea {
namespace Graphics {

/**
 * Mesh based on a doubly-connected edge list (or halfedge data structure). Optionally allows GPU-buffered rendering (see
 * important note below).
 *
 * Adapted from: DcelMesh class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 *
 * @note When using GPU-buffered rendering, any methods of this class which change the mesh will automatically re-initialize the
 * buffers. <b>However</b>, if <i>external</i> methods change the mesh, such as methods of the DcelVertex, DcelEdge and DcelFace
 * classes, then the user <i>must</i> manually indicate that the mesh needs to be resynchronized with the GPU. The
 * invalidateGpuBuffers() function should be used for this.
 */
template < typename VertexAttribute    =  Graphics::NullAttribute,
           typename HalfedgeAttribute  =  Graphics::NullAttribute,
           typename FaceAttribute      =  Graphics::NullAttribute >
class /* THEA_API */ DcelMesh : public virtual NamedObject, public IMesh
{
  public:
    THEA_DECL_SMART_POINTERS(DcelMesh)

    /** Mesh type tag. */
    struct DCEL_MESH_TAG {};

    typedef DcelVertex  <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Vertex;    ///< Vertex of the mesh (3-vector).
    typedef DcelHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute>  Halfedge;  ///< Halfedge of the mesh.
    typedef Halfedge                                                         Edge;      /**< Typedef for interoperability with
                                                                                             GeneralMesh. */
    typedef DcelFace    <VertexAttribute, HalfedgeAttribute, FaceAttribute>  Face;      ///< Face of the mesh.

  private:
    struct HalfedgeComparator
    {
      bool operator()(Halfedge const * e0, Halfedge const * e1) const
      {
        return e0 ? (e1 ? e0->index < e1->index : false) : (bool)e1;
      }
    };

    typedef UnorderedSet<Vertex *>               VertexSet;
    typedef Set<Halfedge *, HalfedgeComparator>  HalfedgeSet;  // store halfedges in sequence so twins are consecutive
    typedef UnorderedSet<Face *>                 FaceSet;

  public:
    typedef Algorithms::RefIterator<typename VertexSet::iterator>          VertexIterator;         ///< Vertex iterator.
    typedef Algorithms::RefIterator<typename VertexSet::const_iterator>    VertexConstIterator;    ///< Const vertex iterator.
    typedef Algorithms::RefIterator<typename HalfedgeSet::iterator>        HalfedgeIterator;       ///< Halfedge iterator.
    typedef Algorithms::RefIterator<typename HalfedgeSet::const_iterator>  HalfedgeConstIterator;  ///< Const halfedge iterator.
    typedef Algorithms::RefIterator<typename FaceSet::iterator>            FaceIterator;           ///< Face iterator.
    typedef Algorithms::RefIterator<typename FaceSet::const_iterator>      FaceConstIterator;      ///< Const face iterator.

    // Generic typedefs, each mesh class must define these for builder and codec compatibility
    typedef Vertex        *  VertexHandle;       ///< Handle to a mesh vertex.
    typedef Vertex const  *  VertexConstHandle;  ///< Handle to an immutable mesh vertex.
    typedef Face          *  FaceHandle;         ///< Handle to a mesh face.
    typedef Face   const  *  FaceConstHandle;    ///< Handle to an immutable mesh face.

    /** Iterator over edges (alternate halfedges starting from the first). */
    typedef DCELInternal::BidirEdgeIterator<HalfedgeIterator> EdgeIterator;

    /** Const iterator over edges (alternate halfedges starting from the first). */
    typedef DCELInternal::BidirEdgeIterator<HalfedgeConstIterator> EdgeConstIterator;

    /** Identifiers for the various buffers (enum class). */
    struct BufferId
    {
      /** Supported values. */
      enum Value
      {
        ALL              =  0xFFFF,  ///< The set of all GPU buffers.
        VERTEX_POSITION  =  0x0001,  ///< IBuffer containing vertex positions.
        VERTEX_NORMAL    =  0x0002,  ///< IBuffer containing vertex normals.
        VERTEX_COLOR     =  0x0004,  ///< IBuffer containing vertex colors.
        VERTEX_TEXCOORD  =  0x0008,  ///< IBuffer containing vertex texture coordinates.
        FACE_NORMAL      =  0x0010,  ///< IBuffer containing face normals (currently not used).
        TOPOLOGY         =  0x0020   ///< IBuffer(s) containing face indices.
      };

      THEA_ENUM_CLASS_BODY(BufferId)

    }; // struct BufferId

    /** Constructor. */
    DcelMesh(std::string const & name = "AnonymousMesh")
    : NamedObject(name),
      next_halfedge_index(0),
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
     * Copy constructor. Creates a deep copy of the mesh (including copies of the attributes). <b>Currently not
     * implemented.</b>
     */
    DcelMesh(DcelMesh const & src) : NamedObject(src)
    {
      throw Error("DcelMesh: Copy constructor not currently implemented");
    }

    ~DcelMesh() { clear(); }

    /**
     * Make an exact copy of the mesh, optionally returning mapping from source to destination vertices/edges/faces. Previous
     * data in the maps is <b>not</b> cleared. <b>Currently not implemented.</b>
     */
    void copyTo(DcelMesh & dst,
                UnorderedMap<Vertex const *, Vertex *> * vertex_map = nullptr,
                UnorderedMap<Edge const *, Edge *> * edge_map = nullptr,
                UnorderedMap<Face const *, Face *> * face_map = nullptr) const
    {
      throw Error("DcelMesh: Copy function not currently implemented");
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
    VertexConstIterator verticesBegin() const { return Algorithms::makeRefIterator(vertices.begin()); }

    /** Get an iterator pointing to the first vertex. */
    VertexIterator verticesBegin() { return Algorithms::makeRefIterator(vertices.begin()); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexConstIterator verticesEnd() const { return Algorithms::makeRefIterator(vertices.end()); }

    /** Get an iterator pointing to the position beyond the last vertex. */
    VertexIterator verticesEnd() { return Algorithms::makeRefIterator(vertices.end()); }

    /** Get an iterator pointing to the first halfedge. */
    HalfedgeConstIterator halfedgesBegin() const { return Algorithms::makeRefIterator(halfedges.begin()); }

    /** Get an iterator pointing to the first halfedge. */
    HalfedgeIterator halfedgesBegin() { return Algorithms::makeRefIterator(halfedges.begin()); }

    /** Get an iterator pointing to the position beyond the last halfedge. */
    HalfedgeConstIterator halfedgesEnd() const { return Algorithms::makeRefIterator(halfedges.end()); }

    /** Get an iterator pointing to the position beyond the last halfedge. */
    HalfedgeIterator halfedgesEnd() { return Algorithms::makeRefIterator(halfedges.end()); }

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
    FaceConstIterator facesBegin() const { return Algorithms::makeRefIterator(faces.begin()); }

    /** Get an iterator pointing to the first face. */
    FaceIterator facesBegin() { return Algorithms::makeRefIterator(faces.begin()); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceConstIterator facesEnd() const { return Algorithms::makeRefIterator(faces.end()); }

    /** Get an iterator pointing to the position beyond the last face. */
    FaceIterator facesEnd() { return Algorithms::makeRefIterator(faces.end()); }

    /** Deletes all data in the mesh and resets automatic element indexing. */
    void clear()
    {
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
        delete *vi;

      for (auto ei = halfedges.begin(); ei != halfedges.end(); ++ei)
        delete *ei;

      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        delete *fi;

      vertices.clear();
      halfedges.clear();
      faces.clear();

      next_halfedge_index = 0;
      max_vertex_index = -1;
      max_face_index = -1;
      bounds = AxisAlignedBox3();

      invalidateGpuBuffers();
    }

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const { return vertices.empty() && faces.empty() && halfedges.empty(); }

    /** Get the number of vertices. */
    intx numVertices() const { return (intx)vertices.size(); }

    /** Get the number of halfedges. */
    intx numHalfedges() const { return (intx)halfedges.size(); }

    /** Get the number of edges. */
    intx numEdges() const
    {
      debugAssertM(halfedges.size() % 2 == 0, getNameStr() + ": Number of halfedges is odd");
      return (intx)(halfedges.size() / 2);
    }

    /** Get the number of faces. */
    intx numFaces() const { return (intx)faces.size(); }

    /** Compute the number of triangles in the mesh. */
    intx numTriangles() const
    {
      intx rval = 0;
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if ((*fi)->isTriangle())
          rval++;

      return rval;
    }

    /** Compute the number of quads in the mesh. */
    intx numQuads() const
    {
      intx rval = 0;
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if ((*fi)->isQuad())
          rval++;

      return rval;
    }

    /** Recompute and cache the bounding box for the mesh. Make sure this has been called before calling getBounds(). */
    void updateBounds()
    {
      bounds = AxisAlignedBox3();
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
        bounds.merge((*vi)->getPosition());
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
     * Add a vertex to the mesh, with optional precomputed normal and index. If the index is negative, a new, unique index is
     * generated for the vertex. Automatically updates the bounding box of the mesh, so you don't need to call updateBounds()
     * after calling this function.
     *
     * @return A pointer to the new vertex on success, null on failure.
     */
    Vertex * addVertex(Vector3 const & point, intx index = -1, Vector3 const * normal = nullptr)
    {
      if (vertices.empty()) bounds.set(point, point);
      else                  bounds.merge(point);

      Vertex * vertex = normal ? new Vertex(point, *normal) : new Vertex(point);
      vertices.insert(vertex);

      if (index < 0)
        index = (++max_vertex_index);
      else if (index > max_vertex_index)
        max_vertex_index = index;

      vertex->setIndex(index);

#ifdef THEA_DCELMESH_VERBOSE
      THEA_CONSOLE << "Added vertex " << index << " at " << vertex->getPosition();
#endif

      invalidateGpuBuffers();
      return vertex;
    }

    /**
     * Add a face to the mesh, specified by the sequence of vertices obtained by dereferencing [vbegin, vend) and an optional
     * index. VertexInputIterator must dereference to a pointer to a Vertex. If the index is negative, a new, unique index is
     * generated for the face. Unless the mesh is already in an inconsistent state, failure to add the face will not affect the
     * mesh.
     *
     * @return A pointer to the new created face, or null on error.
     */
    template <typename VertexInputIterator>
    Face * addFace(VertexInputIterator vbegin, VertexInputIterator vend, intx index = -1)
    {
      if (vbegin == vend)
      {
        THEA_WARNING << getName() << ": Skipping face -- it has no vertices";
        return nullptr;
      }

#ifdef THEA_DCELMESH_VERBOSE
      std::cout << "Adding face:";
#endif

      // Read the vertex pointers into an internal array
      if (face_vertices.size() != 256) face_vertices.resize(256);  // default size, should ensure not too many resizes
      size_t num_verts = 0;
      for (auto vi = vbegin; vi != vend; ++vi, ++num_verts)
      {
        debugAssertM(*vi, getNameStr() + ": Null vertex pointer specified for new face");

        if (num_verts >= face_vertices.size())
          face_vertices.resize(2 * face_vertices.size() + 1);

        face_vertices[num_verts] = *vi;

#ifdef THEA_DCELMESH_VERBOSE
        std::cout << ' ' << (*vi)->getIndex();
#endif
      }

#ifdef THEA_DCELMESH_VERBOSE
      std::cout << std::endl;
#endif

      if (num_verts < 3)
      {
        THEA_WARNING << getName() << ": Skipping face -- too few vertices (" << num_verts << ')';
        return nullptr;
      }

      // Compute the face normal, assume it is consistent across the face
      Vector3 e1 = face_vertices[0]->getPosition() - face_vertices[1]->getPosition();
      Vector3 e2 = face_vertices[2]->getPosition() - face_vertices[1]->getPosition();
      Vector3 normal = e2.cross(e1).normalized();  // counter-clockwise

      Face * face = addFace(num_verts, &face_vertices[0], normal);  // invalidates GPU buffers
      if (face)
      {
        if (index < 0)
          index = (++max_face_index);
        else if (index > max_face_index)
          max_face_index = index;

        face->setIndex(index);
      }

      return face;
    }

    /**
     * Split an edge along its length, in the ratio given by \a frac.
     *
     * @return The new vertex created by the split operation, or null on error.
     */
    Vertex * splitEdge(Halfedge * edge, Real frac)
    {
      alwaysAssertM(frac >= 0 && frac <= 1, getNameStr() + ": Edge split fraction should be between 0 and 1");
      if (edge)
      {
        THEA_ERROR << getName() << "Can't split null edge";
        return nullptr;
      }

      Vector3 p = (1 - frac) * edge->getOrigin()->getPosition() + frac * edge->getEnd()->getPosition();
      Vector3 n = ((1 - frac) * edge->getOrigin()->getNormal() + frac * edge->getEnd()->getNormal()).normalized();
      Vertex * new_vx = addVertex(p, &n);
      if (!new_vx)
        return nullptr;

      if (!splitEdge(edge, new_vx))  // invalidates GPU buffers
      {
        // We should generally never get here
        removeVertex(new_vx);
        return nullptr;
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
      Real s = (pos - edge->getOrigin()->getPosition()).norm();
      Real t = (pos - edge->getEnd()->getPosition()).norm();
      Vector3 n = (t * edge->getOrigin()->getNormal() + s * edge->getEnd()->getNormal()).normalized();
      Vertex * new_vx = addVertex(pos, &n);
      if (!new_vx)
        return nullptr;

      if (!splitEdge(edge, new_vx))  // invalidates GPU buffers
      {
        // We should generally never get here
        removeVertex(new_vx);
        return nullptr;
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
      intx n = 0;
      do
      {
        if (!e0->isBoundaryHalfedge() || !e1->isBoundaryHalfedge())
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
      Halfedge * next0 = nullptr, * next1 = nullptr;
      for (intx i = 0; i < n; ++i, e0 = next0, e1 = next1)
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
        intx index = e0->index;

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

      invalidateGpuBuffers();
      return true;
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
    /** Set vertex color. */
    template < typename VertexT, typename std::enable_if< HasColor<VertexT>::value, int >::type = 0 >
    void setVertexColor(VertexT * vertex, ColorRgba const & color)
    {
      vertex->attr().setColor(color);
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
    }

    /** Set vertex texture coordinates (no-op, called if vertex does not have texture coordinate attribute). */
    template < typename VertexT, typename std::enable_if< !HasTexCoord<VertexT>::value, int >::type = 0 >
    void setVertexTexCoord(VertexT * vertex, Vector2 const & texcoord)
    {}

    /**
     * Add a face to the mesh, specified by a sequence of boundary vertices.
     *
     * @return A pointer to the newly created face.
     */
    Face * addFace(size_t num_verts, Vertex ** verts, Vector3 const & normal)
    {
      // Try to locate a boundary edge that is already in the mesh
      size_t origin = 0;
      Halfedge * e = findExistingEdge(num_verts, verts, origin);

      Face * face = nullptr;
      face = addFace(num_verts, verts, e, origin, normal);
      if (!face)
      {
        THEA_WARNING << "Adding isolated face";
        Array<Vertex *> iso_verts;
        size_t existing0 = 0, existing1 = 0;
        if (e && e->isBoundaryEdge())
        {
          existing0 = origin;
          existing1 = (origin + 1) % num_verts;
        }

        if (!addIsolatedVertices(num_verts, verts, existing0, existing1, iso_verts))
          return nullptr;

        if (e && e->isBoundaryEdge())
          face = addFace(num_verts, const_cast<Vertex **>(&iso_verts[0]), e, origin, normal);
        else
          face = addFace(num_verts, const_cast<Vertex **>(&iso_verts[0]), nullptr, 0, false, normal);
      }

      return face;
    }

    /** Utility function for addFace(size_t, Vertex **, Vector3 const &). */
    Face * addFace(size_t num_verts, Vertex ** verts, Halfedge * first, size_t origin, Vector3 const & normal)
    {
      if (first)
      {
        if (first->isBoundaryHalfedge())
        {
          // All good... keep adding starting here
          return addFace(num_verts, verts, first, origin, false, normal);
        }
        else
        {
          if (first->twin()->isBoundaryHalfedge())
          {
            // Try adding in reverse
            THEA_DEBUG << getName() << ": Trying to add face by reversing the order of vertices";
            return addFace(num_verts, verts, first->twin(), (origin + 1) % num_verts, true, -normal);  // normal got flipped
          }
          else
          {
            THEA_WARNING << getName() << ": Face has edge already adjoining two faces";
            return nullptr;
          }
        }
      }
      else
        return addFace(num_verts, verts, nullptr, 0, false, normal);
    }

    /**
     * Add a sequence of isolated vertices, optionally including one existing edge (if existing0 != existing1), and return
     * pointers to the new vertices.
     */
    bool addIsolatedVertices(size_t num_verts, Vertex ** verts, size_t existing0, size_t existing1,
                             Array<Vertex *> & iso_verts)
    {
      bool has_existing_edge = (existing0 != existing1);
      for (size_t i = 0; i < num_verts; ++i)
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
    Face * addFace(size_t num_verts, Vertex ** verts, Halfedge * first, size_t origin, bool reverse,
                   Vector3 const & normal)
    {
      // For each vertex, check that the halfedge emanating from it exists as a boundary halfedge, or can be successfully added.
      // Store either the existing edge to the next vertex, or the preceding boundary halfedge _into_ the vertex.
      Array<Halfedge *> edges(num_verts);
      Halfedge * prev = first, * e;
      size_t start_index = reverse ? (origin > 0 ? origin - 1 : num_verts - 1)
                                         : (origin < num_verts - 1 ? origin + 1 : 0);
      size_t i = start_index, next;
      do
      {
        // See if there is an existing edge
        next = reverse ? (i > 0 ? i - 1 : num_verts - 1) : (i < num_verts - 1 ? i + 1 : 0);
        e = verts[i]->getEdgeTo(verts[next]);

        if (e)  // edge exists, needs to be boundary edge
        {
          if (!e->isBoundaryHalfedge())
          {
            THEA_WARNING << getName() << ": Can't stitch a face to a halfedge that already adjoins a face";
            return nullptr;
          }

          if (prev && prev->next() != e)
          {
            THEA_WARNING << getName() << ": Face breaks halfedge order at vertex";
            return nullptr;
          }

          edges[i] = e;
          prev = e;
        }
        else  // edge does not exist
        {
          if (prev)
          {
            // Just to be sure, check that prev->next() is a boundary edge
            debugAssertM(prev->isBoundaryHalfedge(),
                         getNameStr() + ": Previous halfedge of new face is not currently on the mesh boundary");
            debugAssertM(prev->next() && prev->next()->isBoundaryHalfedge(),
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
              if (e->twin()->isBoundaryHalfedge())
                edges[i] = e->twin();
              else
              {
                THEA_WARNING << getName() << ": Can't stitch a face to a halfedge that already adjoins a face";
                return nullptr;
              }
            }
            else
              edges[i] = nullptr;
          }

          prev = nullptr;  // no current edge between i and next
        }

        i = next;

      } while (i != start_index);

      // Now create the new face
      Face * face = new Face;
      face->num_edges = (int)num_verts;
      face->setNormal(normal);

      // Now actually create and add the edges. In this loop we directly access private members (ha ha) of the edges to bypass
      // any consistency checks until we've finished adding the face.
      Halfedge * new_e, * new_twin, * next_e, * last = nullptr; (void)last;  // Squash -Wunused-but-set-variable
      Vertex * vi, * vnext;
      intx index0, index1;
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
          debugAssertM(e->isBoundaryHalfedge(),
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

        debugAssertM(!last || last->next_he == new_e, getNameStr() + ": Next pointers on face boundary not consistent");
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

      invalidateGpuBuffers();
      return face;
    }

    /**
     * Get indices for a pair of new twin halfedges. If the indices are getting too large, this function reindexes the set of
     * halfedges.
     */
    void nextHalfedgeIndices(intx & index0, intx & index1)
    {
      static intx const THRESHOLD = std::numeric_limits<intx>::max() - 4;
      if (next_halfedge_index > THRESHOLD)
      {
        if ((intx)halfedges.size() > THRESHOLD)  // too many halfedges, we can't do anything!
          throw Error(getNameStr() + ": Too many halfedges, can't assign indices!");

        next_halfedge_index = 0;
        for (auto ei = halfedges.begin(); ei != halfedges.end(); ++ei)
          (*ei)->index = next_halfedge_index++;
      }

      index0 = next_halfedge_index++;
      index1 = next_halfedge_index++;

      invalidateGpuBuffers(BufferId::TOPOLOGY);
    }

    /**
     * Find an edge that is already in the mesh, from a loop of edges specified by the sequence of their vertices. This function
     * will preferentially return a boundary edge if it finds one.
     *
     * @return A pointer to an existing half-edge, if found, else null. The index of the originating vertex of this half-edge
     *   (respective to the input array of vertices) is stored in origin_index.
     */
    Halfedge * findExistingEdge(size_t num_verts, Vertex ** verts, size_t & origin_index) const
    {
      Halfedge * e  = nullptr;
      size_t last = num_verts - 1;
      for (size_t i = 0; i < num_verts; ++i)
      {
        e = verts[last]->getEdgeTo(verts[i]);
        if (e)
        {
          origin_index = last;

          if (e->isBoundaryHalfedge() || e->twin()->isBoundaryHalfedge())
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
    static Real ccwAngle(Vector3 const & u, Vector3 const & v, Vector3 const & unit_up)
    {
      Vector3 v_proj = v - (v.dot(unit_up) * unit_up);
      Real s = u.cross(v_proj).norm();
      Real c = u.dot(v_proj);
      Real ang = std::atan2(s, c);

      return ang < 0 ? 2 * Math::pi() - ang : ang;
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
    Halfedge * findPrevAroundVertex(Vertex * u, Vertex * v, Vector3 const & normal) const
    {
      Halfedge * first = u->getHalfedge();
      if (!first)
        return nullptr;

      // Compute the new normal, which defines the plane of projection
      Vector3 new_normal = u->estimateUpdatedNormal(normal);
      Halfedge * e = first, * best = first;
      Real ang, best_ang = 100;  // anything > 2 pi
      Vector3 edge = v->getPosition() - u->getPosition();
      edge = (edge - (edge.dot(new_normal) * new_normal));
      do
      {
        if (e->twin()->isBoundaryHalfedge())
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
      invalidateGpuBuffers();
    }

    /** Delete a halfedge from the mesh. No pointers are updated, the halfedge is just deleted from storage. */
    void removeHalfedge(Halfedge * halfedge)
    {
      halfedges.erase(halfedge);
      delete halfedge;
      invalidateGpuBuffers();
    }

    /** Delete a face from the mesh. No pointers are updated, the face is just deleted from storage. */
    void removeFace(Face * face)
    {
      faces.erase(face);
      delete face;
      invalidateGpuBuffers();
    }

    /** Split an edge along its length at a given vertex location. The vertex is assumed to have been newly added. */
    bool splitEdge(Halfedge * edge, Vertex * vertex)
    {
      Halfedge * twin = edge->twin();

      intx twin_index = twin->index;
      intx index0 = -1, index1 = -1;
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
        packed_vertex_positions[i] = (*vi)->getPosition();

      setPackedArrayIsValid(BufferId::VERTEX_POSITION);
    }

    /** Pack vertex positions densely in an array. */
    void packVertexNormals() const
    {
      if (packedArrayIsValid(BufferId::VERTEX_NORMAL)) return;

      packed_vertex_normals.resize(vertices.size());
      size_t i = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi, ++i)
        packed_vertex_normals[i] = (*vi)->getNormal();

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
        packed_vertex_colors[i] = ColorRgba((*vi)->attr().getColor());

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
        packed_vertex_texcoords[i] = (*vi)->attr().getTexCoord();

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

      uint32 index = 0;
      for (auto vi = vertices.begin(); vi != vertices.end(); ++vi)
        (*vi)->setPackingIndex(index++);

      has_large_polys = false;
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
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
        size_t i = 0;
        for (auto ei = edgesBegin(); ei != edgesEnd(); ++ei, i += 2)
        {
          packed_edges[i    ] = ei->getOrigin()->getPackingIndex();
          packed_edges[i + 1] = ei->getEnd()->getPackingIndex();
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

    FaceSet      faces;                ///< Set of mesh faces.
    VertexSet    vertices;             ///< Set of mesh vertices.
    HalfedgeSet  halfedges;            ///< Set of mesh halfedges.
    intx         next_halfedge_index;  ///< Index to be assigned to the next added halfedge.

    intx max_vertex_index;  ///< The largest index of a vertex in the mesh.
    intx max_face_index;    ///< The largest index of a face in the mesh.

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

    mutable Array<Vertex *> face_vertices;  // internal cache for vertex pointers for a face

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
};

template <typename V, typename E, typename F>
bool
DcelMesh<V, E, F>::uploadToGraphicsSystem(IRenderSystem & render_system)
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

    if (vertices.empty() || (faces.empty() && halfedges.empty()))
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

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
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

#ifndef THEA_DCEL_MESH_NO_INDEX_ARRAY
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
     && !vertex_texcoords_buf->updateAttributes(0, (int64)packed_vertex_texcoords.size(), 2, NumericType::UINT32,
                                                &packed_vertex_texcoords[0])) return false;
  }

  setAllGpuBuffersAreValid();

  return true;
}

template <typename V, typename E, typename F>
int8
DcelMesh<V, E, F>::drawBuffered(IRenderSystem & render_system, IRenderOptions const & options) const
{
  if (options.drawEdges() && !buffered_wireframe)
  { THEA_ERROR << getName() << ": Can't draw mesh edges with GPU-buffered wireframe disabled"; return false; }

  if (!const_cast<DcelMesh *>(this)->uploadToGraphicsSystem(render_system)) return false;

  if (!vertex_positions_buf) return true;
  if (!options.drawFaces() && !options.drawEdges()) return true;
  if (!options.drawFaces() && halfedges.size() <= 0) return true;
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

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
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
          for (auto fi = faces.begin(); fi != faces.end(); ++fi)
            if ((*fi)->numEdges() > 4)
            {
              render_system.beginPrimitive(IRenderSystem::Primitive::POLYGON);
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

        render_system.setShader(nullptr);
        render_system.setColorBuffer(nullptr);
        render_system.setTexCoordBuffer(0, nullptr);
        render_system.setNormalBuffer(nullptr);
        render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)

#ifdef THEA_DCEL_MESH_NO_INDEX_ARRAY
        if (!edges.empty())
          render_system.sendIndices(IRenderSystem::Primitive::LINES, (int64)halfedges.size(), &packed_edges[0]);
#else
        if (!halfedges.empty())
        {
          render_system.setIndexBuffer(edges_buf);
          render_system.sendIndicesFromBuffer(IRenderSystem::Primitive::LINES, 0, (int64)halfedges.size());
        }
#endif

      render_system.popColorFlags();
      render_system.popShader();
    }

  render_system.endIndexedPrimitives();

  char const * err = nullptr;
  if ((err = render_system.getAndClearError()))
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

template <typename V, typename E, typename F>
int8
DcelMesh<V, E, F>::drawImmediate(IRenderSystem & render_system, IRenderOptions const & options) const
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
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if ((*fi)->isTriangle()) drawFace(**fi, render_system, options);
    render_system.endPrimitive();

    // Now render all quads
    render_system.beginPrimitive(IRenderSystem::Primitive::QUADS);
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
        if ((*fi)->isQuad()) drawFace(**fi, render_system, options);
    render_system.endPrimitive();

    // Finish off with all larger polygons
    for (auto fi = faces.begin(); fi != faces.end(); ++fi)
      if ((*fi)->numEdges() > 4)
      {
        render_system.beginPrimitive(IRenderSystem::Primitive::POLYGON);
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

      render_system.setShader(nullptr);
      render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)

      render_system.beginPrimitive(IRenderSystem::Primitive::LINES);
        for (auto ei = edgesBegin(); ei != edgesEnd(); ++ei)
        {
          ei->attr().draw(render_system, options);
          render_system.sendVertex(3, ei->getOrigin()->getPosition().data());
          render_system.sendVertex(3, ei->getEnd()->getPosition().data());
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }

  char const * err = nullptr;
  if ((err = render_system.getAndClearError()))
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

} // namespace Graphics
} // namespace Thea

#endif
