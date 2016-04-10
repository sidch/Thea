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

#ifndef __Thea_Graphics_CGALMesh_hpp__
#define __Thea_Graphics_CGALMesh_hpp__

#ifdef THEA_ENABLE_CGAL

#include "../Common.hpp"
#include "../Array.hpp"
#include "../AttributedObject.hpp"
#include "../NamedObject.hpp"
#include "../UnorderedMap.hpp"
#include "DrawableObject.hpp"
#include "GraphicsAttributes.hpp"
#include "IncrementalCGALMeshBuilder.hpp"
#include "DefaultMeshCodecs.hpp"
#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace Thea {
namespace Graphics {

// [Internal] A custom traits class for an attributed type.
template <typename PointT>
struct /* THEA_DLL_LOCAL */ _CGALMeshPolyhedronTraits_3
{
  typedef PointT Point_3;

  struct Plane_3 {};  // to avoid syncing problems, force recomputation of plane equations for faces as needed (maybe add a
                      // user-controlled caching mechanism later?)

  struct Construct_opposite_plane_3
  {
    Plane_3 operator()(Plane_3 plane) { return plane; }
  };

  Construct_opposite_plane_3 construct_opposite_plane_3_object() { return Construct_opposite_plane_3(); }

}; // class _CGALMeshPolyhedronTraits_3

// [Internal] Specification of vertex, halfedge and face types for a CGALMesh.
template <typename VertexAttributeT, typename HalfedgeAttributeT, typename FaceAttributeT>
struct /* THEA_DLL_LOCAL */ _CGALMeshPolyhedronItems_3 : public CGAL::Polyhedron_items_3
{
  // Wrapper for a vertex with position and attribute.
  template <typename Refs, typename Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    struct Vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>,
                    public AttributedObject<VertexAttributeT>
    {
      typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> VertexBaseT;

      Vertex() {};
      Vertex(Point const & p) : VertexBaseT(p) {}
    };
  };

  // Wrapper for an edge with an attribute.
  template <typename Refs, typename Traits>
  struct Halfedge_wrapper
  {
    struct Halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>, public AttributedObject<HalfedgeAttributeT> {};
  };

  // Wrapper for a face with an attribute.
  template <typename Refs, typename Traits>
  struct Face_wrapper
  {
    typedef typename Traits::Plane_3 Plane3;
    struct Face : public CGAL::HalfedgeDS_face_base<Refs, CGAL::Tag_true, Plane3>, public AttributedObject<FaceAttributeT>
    {};
  };

}; // class _CGALMeshPolyhedronItems_3

/**
 * A general CGAL mesh object, inherits CGAL::Polyhedron_3. The base interface is documented in the CGAL manual at
 * http://www.cgal.org/Manual . This wrapper is intended to make customizing the mesh with user-specified attributes a little
 * easier.
 *
 * Attributes (in addition to vertex locations, specified as points) may be attached to vertices, half-edges and faces, all of
 * which derive from AttributedObject. An attribute class must have the member function <code>%draw(RenderSystem &,
 * RenderOptions const &) const</code> defined (preferably, use a class derived from GraphicsAttribute). This function will be
 * called <em>within</em> a RenderSystem::beginPrimitive() / RenderSystem::endPrimitive() block, so it should be valid in this
 * context. Typically, attributes with non-trivial drawing functions include normals (call RenderSystem::sendNormal()) and
 * texture coordinates (call RenderSystem::sendTexCoord()). Some standard graphics attributes are defined in
 * GraphicsAttributes.hpp .
 *
 * This template is typically instantiated as
 * <code>CGALMesh< Vector3, NullAttribute, NormalTexCoordAttribute < Vector3, Vector3 > ></code>. This defines a 3D mesh with a
 * normal and a texture coordinate for each vertex-face incidence (so that discontinuities in smoothness and color can be
 * represented).
 */
template < typename PointT              =  Vector3,
           typename VertexAttributeT    =  Graphics::NullAttribute,
           typename HalfedgeAttributeT  =  Graphics::NullAttribute,
           typename FaceAttributeT      =  Graphics::NullAttribute >
class CGALMesh : public CGAL::Polyhedron_3< _CGALMeshPolyhedronTraits_3<PointT>,
                                            _CGALMeshPolyhedronItems_3<VertexAttributeT, HalfedgeAttributeT, FaceAttributeT> >,
                 public virtual NamedObject,
                 public DrawableObject
{
  public:
    THEA_DEF_POINTER_TYPES(CGALMesh, shared_ptr, weak_ptr)

    /** Mesh type tag. */
    typedef int CGAL_MESH_TAG;

    // The base mesh type.
    typedef CGAL::Polyhedron_3< _CGALMeshPolyhedronTraits_3<PointT>,
                                _CGALMeshPolyhedronItems_3<VertexAttributeT, HalfedgeAttributeT, FaceAttributeT> >
            BaseT;

    /** Incrementally builds the mesh. */
    typedef CGAL::Polyhedron_incremental_builder_3<typename BaseT::HalfedgeDS> Builder;
    typedef shared_ptr<Builder> BuilderPtr;

    typedef PointT              Point;              ///< Point type.
    typedef VertexAttributeT    VertexAttribute;    ///< Vertex attribute type.
    typedef HalfedgeAttributeT  HalfedgeAttribute;  ///< Half-edge attribute type.
    typedef FaceAttributeT      FaceAttribute;      ///< Face attribute type.
    typedef FaceAttributeT      FacetAttribute;     ///< Face attribute type (synonym for FaceAttribute).

    /** Create an empty mesh. */
    CGALMesh(std::string const & name = "AnonymousMesh") : NamedObject(name) {}

    /** Copy constructor. */
    CGALMesh(CGALMesh const & src) : NamedObject(src), BaseT(src), bounds(src.bounds) {}

    /** Create a mesh from a set of faces. FacetInputIterator must double-dereference to Facet. */
    template <typename FacetInputIterator>
    CGALMesh(FacetInputIterator src_facets_begin, FacetInputIterator src_facets_end, std::string const & name = "AnonymousMesh")
    : NamedObject(name)
    {
      // Create an array to hold the set of vertex indices for a face
      TheaArray<size_t> indices(256);  // default size, should ensure not too many resizes

      // Create a mesh builder
      BuilderPtr bp = createBuilder();
      Builder & b = *bp;  // avoid repeated dereferencing of shared pointer

      b.begin_surface(0, 0, 0, Builder::ABSOLUTE_INDEXING);
      if (b.error()) throw Error(getNameStr() + ": Builder error at start of construction from facets");

        // Track the indices of copied vertices
        typedef TheaUnorderedMap<typename BaseT::Vertex const *, size_t> VertexIndexMap;
        VertexIndexMap vertex_indices;

        // Track the images of copied edges
        typedef TheaUnorderedMap<typename BaseT::Halfedge const *, typename BaseT::Halfedge *> HalfedgeMap;
        HalfedgeMap halfedge_map;

        for (FacetInputIterator fi = src_facets_begin; fi != src_facets_end; ++fi)
        {
          if ((*fi)->facet_degree() > (size_t)indices.size())
            indices.resize((array_size_t)(*fi)->facet_degree());

          array_size_t num_vertices = 0;
          typename BaseT::Halfedge_around_facet_const_circulator hc = (*fi)->facet_begin();
          do
          {
            typename BaseT::Vertex const * vp = &(*hc->vertex());
            typename VertexIndexMap::const_iterator existing_index = vertex_indices.find(vp);
            if (existing_index == vertex_indices.end())
            {
              typename BaseT::Vertex_handle new_vertex = b.add_vertex(vp->point());
              if (b.error()) throw Error(getNameStr() + ": Builder error adding vertex");

              // Copy vertex attribute
              new_vertex->setAttr(vp->attr());

              size_t index = vertex_indices.size();
              indices[num_vertices++] =  index;
              vertex_indices[vp]      =  index;
            }
            else
              indices[num_vertices++] = existing_index->second;

          } while (++hc != (*fi)->facet_begin());

          typename BaseT::Halfedge_handle new_facet = b.add_facet(indices.begin(), indices.begin() + num_vertices);
          if (b.error()) throw Error(getNameStr() + ": Builder error adding facet");

          // Copy facet attribute
          new_facet->facet()->setAttr((*fi)->attr());

          // Map source halfedges to their images, for batched attribute copy at the end
          debugAssertM(hc == (*fi)->facet_begin(), getNameStr() + ": Circulator has not looped back to beginning of facet");

          typename BaseT::Halfedge_around_facet_circulator new_hc = new_facet->facet_begin();
          do
          {
            halfedge_map[&(*hc)] = &(*new_hc);

            if (new_hc->opposite()->is_border())  // reduces some overwriting
              halfedge_map[&(*hc->opposite())] = &(*new_hc->opposite());

          } while (++new_hc != new_facet->facet_begin());
        }

        // Copy halfedge attributes
        for (typename HalfedgeMap::iterator hi = halfedge_map.begin(); hi != halfedge_map.end(); ++hi)
          hi->second->setAttr(hi->first->attr());

      b.end_surface();
      if (b.error()) throw Error(getNameStr() + ": Builder error constructing component");

      updateBounds();
    }

    /** Assignment operator. */
    CGALMesh & operator=(CGALMesh const & src)
    {
      DrawableObject::operator=(src);
      BaseT::operator=(src);
      bounds = src.bounds;

      return *this;
    }

    /** Check if the mesh is empty (contains no vertices). */
    bool isEmpty() const { return BaseT::size_of_vertices() <= 0; }

    /** Create an incremental mesh builder for this mesh. */
    BuilderPtr createBuilder() { return BuilderPtr(new Builder(BaseT::hds, true)); }

    void uploadToGraphicsSystem(RenderSystem & render_system) {}

    void draw(RenderSystem & render_system, RenderOptions const & options = RenderOptions::defaults()) const
    {
      if (options.drawFaces())
      {
        if (options.drawEdges())
        {
          render_system.pushShapeFlags();
          render_system.setPolygonOffset(true, 1);
        }

        // For now, ignore all facets with fewer than 3 edges (i.e. no Primitive::LINES). Also note that is_pure_triangle(),
        // is_pure_quad() etc are (currently: 10/2009) implemented by looping over all facets to verify the predicate, so they
        // can't be used as fast tests to eliminate some of the loops below when the polyhedron consists only of triangles or
        // only of quads.

        // Three separate passes over the facet list is probably faster (TODO: profile) than using Primitive::POLYGON for each
        // facet.

        // First try to render as much stuff using triangles as possible
        render_system.beginPrimitive(RenderSystem::Primitive::TRIANGLES);
          for (typename BaseT::Facet_const_iterator fi = BaseT::facets_begin(); fi != BaseT::facets_end(); ++fi)
            if (fi->is_triangle()) drawFace(*fi, render_system, options);
        render_system.endPrimitive();

        // Now render all quads
        render_system.beginPrimitive(RenderSystem::Primitive::QUADS);
          for (typename BaseT::Facet_const_iterator fi = BaseT::facets_begin(); fi != BaseT::facets_end(); ++fi)
            if (fi->is_quad()) drawFace(*fi, render_system, options);
        render_system.endPrimitive();

        // Finish off with all larger polygons
        for (typename BaseT::Facet_const_iterator fi = BaseT::facets_begin(); fi != BaseT::facets_end(); ++fi)
          if (fi->facet_degree() > 4)
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
            for (typename BaseT::Edge_const_iterator ei = BaseT::edges_begin(); ei != BaseT::edges_end(); ++ei)
            {
              ei->attr().draw(render_system, options);  // edge attributes (e.g. color)
              render_system.sendVertex(ei->vertex()->point());
              render_system.sendVertex(ei->opposite()->vertex()->point());
            }
          render_system.endPrimitive();

        render_system.popColorFlags();
        render_system.popShader();
      }
    }

    /**
     * {@inheritDoc}
     *
     * The default definition works for all types that can be used to construct and expand an AxisAlignedBox3 (typically Real,
     * Vector2, Vector3). You should add a specialization for any other Point type that you use.
     */
    void updateBounds()
    {
      bounds = AxisAlignedBox3();
      for (typename BaseT::Point_const_iterator pi = BaseT::points_begin(); pi != BaseT::points_end(); ++pi)
        bounds.merge(*pi);
    }

    AxisAlignedBox3 const & getBounds() const { return bounds; }

  private:
    /**
     * Utility function to draw a facet. Must be enclosed in the appropriate
     * RenderSystem::beginPrimitive()/RenderSystem::endPrimitive() block.
     */
    void drawFace(typename BaseT::Facet const & face, RenderSystem & render_system, RenderOptions const & options) const
    {
      face.attr().draw(render_system, options);  // send any attributes for the face first (e.g. a face normal)
      typename BaseT::Halfedge_around_facet_const_circulator hc = face.facet_begin();
      do
      {
        hc->vertex()->attr().draw(render_system, options);  // vertex attributes (e.g. per-vertex normal or texcoord)
        render_system.sendVertex(hc->vertex()->point());  // finally send the position of the vertex

      } while (++hc != face.facet_begin());
    }

    AxisAlignedBox3 bounds;

}; // class CGALMesh

} // namespace Graphics
} // namespace Thea

#endif

#endif
