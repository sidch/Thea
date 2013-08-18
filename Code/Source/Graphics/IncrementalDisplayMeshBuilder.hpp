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

#ifndef __Thea_Graphics_IncrementalDisplayMeshBuilder_hpp__
#define __Thea_Graphics_IncrementalDisplayMeshBuilder_hpp__

#include "IncrementalMeshBuilder.hpp"
#include "MeshType.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Graphics {

/**
 * Incrementally constructs a display mesh from vertex and face data.
 *
 * @see DisplayMesh
 */
template <typename MeshT>
class IncrementalMeshBuilder<MeshT, typename boost::enable_if< IsDisplayMesh<MeshT> >::type>
{
  public:
    THEA_DEF_POINTER_TYPES(IncrementalMeshBuilder, shared_ptr, weak_ptr)

    typedef MeshT Mesh;                      ///< Type of mesh being built.
    typedef long VertexHandle;               ///< Handle to a mesh vertex.
    typedef typename Mesh::Face FaceHandle;  ///< Handle to a mesh face.

    /** Construct from a raw mesh pointer. Ensure the mesh exists till you've finished using this builder object. */
    IncrementalMeshBuilder(Mesh * mesh_) : mesh(mesh_), num_vertices(0), num_faces(0), building(false)
    {
      alwaysAssertM(mesh, "IncrementalMeshBuilder: Mesh pointer cannot be null");
    }

    /**
     * Construct from a shared mesh pointer. The builder takes shared ownership of the mesh so it survives till the builder is
     * finished.
     */
    IncrementalMeshBuilder(shared_ptr<Mesh> mesh_)
    : mp(mesh_), mesh(mesh_.get()), num_vertices(0), num_faces(0), building(false)
    {
      alwaysAssertM(mesh, "IncrementalMeshBuilder: Mesh pointer cannot be null");
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
      building = true;
    }

    /** Add a vertex to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    VertexHandle addVertex(Vector3 const & pos, Vector3 const * normal = NULL, ColorRGBA const * color = NULL,
                           Vector2 const * texcoord = NULL)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A vertex cannot be added outside a begin/end block");

      VertexHandle ref = mesh->addVertex(pos, normal, color, texcoord);
      num_vertices++;
      return ref;
    }

    /** Add a face to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    template <typename VertexInputIterator> FaceHandle addFace(VertexInputIterator begin, VertexInputIterator end)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A face cannot be added outside a begin/end block");

      FaceHandle ref = mesh->addFace(begin, end);
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
      building = false;
    }

  private:
    shared_ptr<Mesh> mp;  // used for shared ownership
    Mesh * mesh;
    long num_vertices, num_faces;
    bool building;

}; // class IncrementalMeshBuilder<DisplayMesh>

} // namespace Graphics
} // namespace Thea

#endif
