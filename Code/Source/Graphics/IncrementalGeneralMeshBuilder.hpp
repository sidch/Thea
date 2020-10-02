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

#ifndef __Thea_Graphics_IncrementalGeneralMeshBuilder_hpp__
#define __Thea_Graphics_IncrementalGeneralMeshBuilder_hpp__

#include "IncrementalMeshBuilder.hpp"
#include "MeshType.hpp"
#include <type_traits>

namespace Thea {
namespace Graphics {

/**
 * Incrementally constructs a general mesh from vertex and face data.
 *
 * @see GeneralMesh
 */
template <typename MeshT>
class IncrementalMeshBuilder<MeshT, typename std::enable_if< IsGeneralMesh<MeshT>::value >::type>
{
  public:
    THEA_DECL_SMART_POINTERS(IncrementalMeshBuilder)

    typedef MeshT                         Mesh;          ///< Type of mesh being built.
    typedef typename MeshT::VertexHandle  VertexHandle;  ///< Handle to a mesh vertex.
    typedef typename MeshT::FaceHandle    FaceHandle;    ///< Handle to a mesh face.

    /** Construct from a raw mesh pointer. Ensure the mesh exists till you've finished using this builder object. */
    IncrementalMeshBuilder(Mesh * mesh_) : mesh(mesh_), num_vertices(0), num_faces(0), building(false)
    {
      alwaysAssertM(mesh, "IncrementalMeshBuilder: Mesh pointer cannot be null");
    }

    /**
     * Construct from a shared mesh pointer. The builder takes shared ownership of the mesh so it survives till the builder is
     * finished.
     */
    IncrementalMeshBuilder(std::shared_ptr<Mesh> mesh_)
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
    VertexHandle addVertex(Vector3 const & pos, intx index = -1, Vector3 const * normal = nullptr,
                           ColorRgba const * color = nullptr, Vector2 const * texcoord = nullptr)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A vertex cannot be added outside a begin/end block");

      VertexHandle ref = mesh->addVertex(pos, index, normal, color, texcoord);
      num_vertices++;
      return ref;
    }

    /** Add a face to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    template <typename VertexInputIterator>
    FaceHandle addFace(VertexInputIterator begin, VertexInputIterator end, intx index = -1)
    {
      debugAssertM(building, "IncrementalMeshBuilder: A face cannot be added outside a begin/end block");

      FaceHandle ref = mesh->addFace(begin, end, index);
      num_faces++;
      return ref;
    }

    /** Get the number of vertices added so far. */
    intx numVertices() const { return num_vertices; }

    /** Get the number of faces added so far. */
    intx numFaces() const { return num_faces; }

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
    std::shared_ptr<Mesh> mp;  // used for shared ownership
    Mesh * mesh;
    intx num_vertices, num_faces;
    bool building;

}; // class IncrementalMeshBuilder<GeneralMesh>

} // namespace Graphics
} // namespace Thea

#endif
