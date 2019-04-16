//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_MeshKDTree_hpp__
#define __Thea_Algorithms_MeshKDTree_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "KDTreeN.hpp"
#include "MeshTriangles.hpp"

namespace Thea {
namespace Algorithms {

/**
 * A kd-tree on mesh triangles. Implemented for general, DCEL and display meshes.
 *
 * @see GeneralMesh, DCELMesh, DisplayMesh
 */
template <typename MeshT, typename NodeAttributeT = NullAttribute>
class MeshKDTree : public Algorithms::KDTreeN< Triangle3< MeshVertexTriple<MeshT> >, 3, Real, NodeAttributeT >
{
  private:
    typedef KDTreeN< Triangle3< MeshVertexTriple<MeshT> >, 3, Real, NodeAttributeT > BaseT;
    typedef MeshTriangles<MeshT> Triangles;

  public:
    THEA_DEF_POINTER_TYPES(MeshKDTree, std::shared_ptr, std::weak_ptr)

    typedef MeshT Mesh;                                       ///< The mesh type.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;              ///< A group of meshes.
    typedef typename Triangles::VertexTriple VertexTriple;    ///< A triple of mesh vertices.
    typedef typename Triangles::Triangle Triangle;            ///< The triangle defined by a triple of mesh vertices.
    typedef typename Triangles::TriangleArray TriangleArray;  ///< An array of mesh triangles.

    /**
     * Add a mesh to the kd-tree. The mesh is converted to triangles which are cached internally. The tree is <b>not</b>
     * actually constructed until you call init().
     */
    void add(Mesh & mesh)
    {
      tris.add(mesh);
    }

    /**
     * Add a group of meshes to the kd-tree. The meshes are converted to triangles which are cached internally. The tree is
     * <b>not</b> actually constructed until you call init().
     */
    void add(MeshGroup & mg)
    {
      tris.add(mg);
    }

    /** Add a mesh face to the kd-tree. The tree is <b>not</b> updated until you call init(). */
    void addFace(Mesh & mesh, typename Mesh::Face & face)
    {
      tris.addFace(mesh, face);
    }

    /** Add a single triangle to the kd-tree. The tree is <b>not</b> updated until you call init(). */
    void addTriangle(Triangle const & tri)
    {
      tris.addTriangle(tri);
    }

    /**
     * Add a set of triangles to the kd-tree. TriangleIterator must dereference to Triangle. The tree is <b>not</b> updated
     * until you call init().
     */
    template <typename TriangleIterator> void addTriangles(TriangleIterator tris_begin, TriangleIterator tris_end)
    {
      tris.addTriangles(tris_begin, tris_end);
    }

    /**
     * Get direct access to the triangles cached by the tree, which will be used to build the tree on the next call to init().
     */
    TriangleArray const & getTriangles() const { return tris.getTriangles(); }

    /**
     * Get direct access to the triangles cached by the tree, which will be used to build the tree on the next call to init().
     */
    TriangleArray & getTriangles() { return tris.getTriangles(); }

    /**
     * Compute the kd-tree from the added meshes. You <b>must</b> call this function to construct (or recompute) the tree after
     * any addMesh() or addMeshGroup() calls. Clears the triangle cache.
     */
    void init(int max_depth = -1, int max_elems_in_leaf = -1, bool save_memory = false, bool deallocate_previous_memory = true)
    {
      TriangleArray const & tri_array = tris.getTriangles();
      BaseT::init(tri_array.begin(), tri_array.end(), max_depth, max_elems_in_leaf, save_memory, deallocate_previous_memory);
      tris.clear();
    }

    /**
     * Clear the tree. If \a deallocate_all_memory is false, memory allocated in pools is held to be reused if possible by the
     * next init() operation.
     */
    void clear(bool deallocate_all_memory = true)
    {
      BaseT::clear(deallocate_all_memory);
      tris.clear();
    }

  private:
    Triangles tris;  ///< Internal cache of triangles used to initialize the tree.

}; // class MeshKDTree

} // namespace Algorithms
} // namespace Thea

#endif
