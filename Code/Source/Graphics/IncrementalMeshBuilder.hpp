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

#ifndef __Thea_Graphics_IncrementalMeshBuilder_hpp__
#define __Thea_Graphics_IncrementalMeshBuilder_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Graphics {

/** Incrementally constructs a mesh from vertex and face data. */
template <typename MeshT, typename Enable = void>
class IncrementalMeshBuilder
{
  private:
    typedef int VertexHandle;  // dummy
    typedef int FaceHandle;    // dummy

  public:
    THEA_DECL_SMART_POINTERS(IncrementalMeshBuilder)

    typedef MeshT Mesh;  ///< Type of mesh being built.

    /** Get the mesh being built. */
    Mesh const * getMesh() const;

    /** Get the mesh being built. */
    Mesh * getMesh();

    /**
     * Start building the mesh. For every begin there must be a corresponding end(). Depending on the type of mesh, you should
     * not try to modify the mesh through other means until the corresponding call to end(). A mesh may be incrementally built
     * in piecewise fashion through multiple begin() / end() blocks. Blocks may not be nested.
     *
     * @see end()
     */
    void begin();

    /**
     * Add a vertex to the mesh and return a handle to it. Must be called within a begin() / end() block. The return value can
     * be checked for errors with isValidVertexHandle().
     */
    VertexHandle addVertex(Vector3 const & pos, intx index = -1, Vector3 const * normal = nullptr,
                           ColorRgba const * color = nullptr, Vector2 const * texcoord = nullptr);

    /**
     * Add a face to the mesh and return a handle to it. Must be called within a begin() / end() block. The return value can be
     * checked for errors with isValidFaceHandle().
     */
    template <typename IndexIterator> FaceHandle addFace(IndexIterator begin, IndexIterator end, intx index = -1);

    /**
     * Check if a vertex handle has a valid value and is not a null pointer (for pointer-based mesh classes) or a negative
     * integer (for index-based mesh classes).
     */
    static bool isValidVertexHandle(VertexHandle handle);

    /**
     * Check if a face handle has a valid value and is not a null pointer (for pointer-based mesh classes) or a negative integer
     * (for index-based mesh classes).
     */
    static bool isValidFaceHandle(FaceHandle handle);

    /** Get the number of vertices added so far. */
    intx numVertices() const;

    /** Get the number of faces added so far. */
    intx numFaces() const;

    /**
     * Complete the current build process. Must be matched to begin().
     *
     * @see begin()
     */
    void end();

}; // class IncrementalMeshBuilder

} // namespace Graphics
} // namespace Thea

#endif
