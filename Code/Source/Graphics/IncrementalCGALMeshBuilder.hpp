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

#ifndef __Thea_Graphics_IncrementalCGALMeshBuilder_hpp__
#define __Thea_Graphics_IncrementalCGALMeshBuilder_hpp__

#ifdef THEA_ENABLE_CGAL

#include "IncrementalMeshBuilder.hpp"
#include "MeshType.hpp"
#include "../UnorderedMap.hpp"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Graphics {

/**
 * Incrementally constructs a CGAL mesh (CGALMesh, CGAL::Polyhedron_3, or any suitable derived class) from vertex and face data.
 */
template <typename MeshT>
class IncrementalMeshBuilder<MeshT, typename boost::enable_if< IsCGALMesh<MeshT> >::type>
{
  public:
    typedef MeshT Mesh;  ///< Type of mesh being built.
    typedef typename Mesh::Vertex_handle VertexHandle;  ///< Handle to a mesh vertex.
    typedef typename Mesh::Face_handle FaceHandle;  ///< Handle to a mesh face.

    THEA_DEF_POINTER_TYPES(IncrementalMeshBuilder, shared_ptr, weak_ptr)

  private:
    typedef TheaUnorderedMap<typename Mesh::Vertex *, size_t> VertexIndexMap;

    /** Delegate class which and holds a modifiable reference to the mesh and wraps a Polyhedron_3 builder. */
    template <typename HDS>
    class CGALBuilder : public CGAL::Modifier_base<HDS>
    {
      private:
        typedef CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
        Builder * builder;
        VertexIndexMap * vertex_indices;
        size_t next_index;

      public:
        CGALBuilder(VertexIndexMap * vertex_indices_) : builder(NULL), vertex_indices(vertex_indices_) {}
        ~CGALBuilder() { delete builder; }

        void operator()(HDS & hds)
        {
          next_index = hds.size_of_vertices();
          builder = new Builder(hds, true);
        }

        void begin() { builder->begin_surface(0, 0, 0, Builder::ABSOLUTE_INDEXING); }

        VertexHandle addVertex(Vector3 const & pos)
        {
          VertexHandle ref = builder->add_vertex(typename Mesh::Point_3(pos.x(), pos.y(), pos.z()));
          if (builder->error()) throw Error("IncrementalCGALMeshBuilder: Error adding vertex");

          vertex_indices->insert(typename VertexIndexMap::value_type(&(*ref), next_index));
          next_index++;
          return ref;
        }

        template <typename VertexInputIterator> FaceHandle addFace(VertexInputIterator begin, VertexInputIterator end)
        {
          FaceHandle ref = builder->begin_facet();

          for (VertexInputIterator vi = begin; vi != end; ++vi)
          {
            typename VertexIndexMap::const_iterator existing = vertex_indices->find(&(**vi));
            if (existing != vertex_indices->end())
              builder->add_vertex_to_facet(existing->second);
            else
              throw Error("IncrementalCGALMeshBuilder: Vertex reference cannot be mapped to a known index");
          }

          builder->end_facet();
          if (builder->error()) throw Error("IncrementalCGALMeshBuilder: Error adding vertex");

          return ref;
        }

        void end() { builder->end_surface(); }

    }; // class CGALBuilder

  public:
    /** Construct from a raw mesh pointer. Ensure the mesh exists till you've finished using this builder object. */
    IncrementalMeshBuilder(Mesh * mesh_) : mesh(mesh_), builder(NULL), num_vertices(0), num_faces(0), building(false)
    {
      alwaysAssertM(mesh, "IncrementalMeshBuilder: Mesh pointer cannot be null");
      builder = new CGALBuilder<typename Mesh::HalfedgeDS>(&vertex_indices);
      mesh->delegate(*builder);
    }

    /**
     * Construct from a shared mesh pointer. The builder takes shared ownership of the mesh so it survives till the builder is
     * finished.
     */
    IncrementalMeshBuilder(shared_ptr<Mesh> mesh_)
    : mp(mesh_), mesh(mesh_.get()), builder(NULL), num_vertices(0), num_faces(0), building(false)
    {
      alwaysAssertM(mesh, "IncrementalMeshBuilder: Mesh pointer cannot be null");
      builder = new CGALBuilder<typename Mesh::HalfedgeDS>(&vertex_indices);
      mesh->delegate(*builder);
    }

    /**
     * Start building the mesh. A mesh may be incrementally built in piecewise fashion through multiple begin() / end()
     * blocks. For every begin there must be a corresponding end(). Blocks may not be nested.
     *
     * @see end()
     */
    void begin()
    {
      alwaysAssertM(!building, "IncrementalMeshBuilder: begin/end not matched");
      builder->begin();
      building = true;
    }

    /** Add a vertex to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    VertexHandle addVertex(Vector3 const & pos, Vector3 const * normal = NULL, ColorRGBA const * color = NULL,
                           Vector2 const * texcoord = NULL)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A vertex cannot be added outside a begin/end block");

      // CGALMesh doesn't accept normals/colors/texcoords by default
      VertexHandle ref = builder->addVertex(pos);
      num_vertices++;
      return ref;
    }

    /** Add a face to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    template <typename VertexInputIterator> FaceHandle addFace(VertexInputIterator begin, VertexInputIterator end)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A face cannot be added outside a begin/end block");

      FaceHandle ref = builder->addFace(begin, end);
      num_faces++;
      return ref;
    }

    /** Get the number of vertices added so far. */
    long numVertices() const { return num_vertices; }

    /** Get the number of faces added so far. */
    long numFaces() const { return num_faces; }

    /**
     * Complete the current build process. Must be matched to begin().
     *
     * @see begin()
     */
    void end()
    {
      alwaysAssertM(building, "IncrementalMeshBuilder: begin/end not matched");
      builder->end();
      building = false;
    }

  private:
    shared_ptr<Mesh> mp;  // used for shared ownership
    Mesh * mesh;
    VertexIndexMap vertex_indices;
    CGALBuilder<typename Mesh::HalfedgeDS> * builder;
    long num_vertices, num_faces;
    bool building;

}; // class IncrementalMeshBuilder<CGALMesh>

} // namespace Graphics
} // namespace Thea

#endif

#endif
