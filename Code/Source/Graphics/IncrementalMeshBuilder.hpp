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

#ifndef __Thea_Graphics_IncrementalMeshBuilder_hpp__
#define __Thea_Graphics_IncrementalMeshBuilder_hpp__

#include "../Common.hpp"
#include "../Vector3.hpp"

namespace Thea {
namespace Graphics {

/** Incrementally constructs a mesh from vertex and face data. */
template <typename MeshT, typename Enable = void>
class IncrementalMeshBuilder
{
  private:
    typedef int VertexHandle;  // dummy
    typedef int FaceHandle;  // dummy

  public:
    THEA_DEF_POINTER_TYPES(IncrementalMeshBuilder, shared_ptr, weak_ptr)

    typedef MeshT Mesh;  ///< Type of mesh being built.

    /**
     * Start building the mesh. For every begin there must be a corresponding end(). Depending on the type of mesh, you should
     * not try to modify the mesh through other means until the corresponding call to end(). A mesh may be incrementally built
     * in piecewise fashion through multiple begin() / end() blocks. Blocks may not be nested.
     *
     * @see end()
     */
    void begin();

    /** Add a vertex to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    VertexHandle addVertex(Vector3 const & pos, Vector3 const * normal = NULL, Color4 const * color = NULL,
                           Vector2 const * texcoord = NULL);

    /** Add a face to the mesh and return a handle to it. Must be called within a begin() / end() block. */
    template <typename IndexIterator> FaceHandle addFace(IndexIterator begin, IndexIterator end);

    /** Get the number of vertices added so far. */
    long numVertices() const;

    /** Get the number of faces added so far. */
    long numFaces() const;

    /**
     * Complete the current build process. Must be matched to begin().
     *
     * @see begin()
     */
    void end() {}

}; // class IncrementalMeshBuilder

} // namespace Graphics
} // namespace Thea

#endif
