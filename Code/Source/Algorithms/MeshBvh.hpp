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
// First version: 2011
//
//============================================================================

#ifndef __Thea_Algorithms_MeshBvh_hpp__
#define __Thea_Algorithms_MeshBvh_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "BvhN.hpp"
#include "MeshTriangles.hpp"

namespace Thea {
namespace Algorithms {

/**
 * A bounding volume hierarchy on mesh triangles. Implemented for general, DCEL and display meshes.
 *
 * @see GeneralMesh, DcelMesh, DisplayMesh
 */
template < typename MeshT,
           typename NodeAttributeT = NullAttribute,
           int MaxDegree = 2,
           typename BoundingVolumeT = AxisAlignedBox3 >
class MeshBvh
: public Algorithms::BvhN< TriangleN< 3, MeshVertexTriple<MeshT>, Real >, 3, Real, NodeAttributeT, MaxDegree, AxisAlignedBox3 >
{
  private:
    typedef BvhN< TriangleN< 3, MeshVertexTriple<MeshT>, Real >, 3, Real, NodeAttributeT, MaxDegree, AxisAlignedBox3 > BaseT;
    typedef MeshTriangles<MeshT> Triangles;

  public:
    THEA_DECL_SMART_POINTERS(MeshBvh)

    typedef MeshT Mesh;                                       ///< The mesh type.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;              ///< A group of meshes.
    typedef typename Triangles::VertexTriple VertexTriple;    ///< A triple of mesh vertices.
    typedef typename Triangles::Triangle Triangle;            ///< The triangle defined by a triple of mesh vertices.
    typedef typename Triangles::TriangleArray TriangleArray;  ///< An array of mesh triangles.

    /**
     * Add a mesh to the BVH. The mesh is converted to triangles which are cached internally. The tree is <b>not</b>
     * actually constructed until you call init().
     */
    void add(Mesh & mesh)
    {
      tris.add(mesh);
    }

    /**
     * Add a group of meshes to the BVH. The meshes are converted to triangles which are cached internally. The tree is
     * <b>not</b> actually constructed until you call init().
     */
    void add(MeshGroup & mg)
    {
      tris.add(mg);
    }

    /** Add a mesh face to the BVH. The tree is <b>not</b> updated until you call init(). */
    void addFace(Mesh & mesh, typename Mesh::Face & face)
    {
      tris.addFace(mesh, face);
    }

    /** Add a single triangle to the BVH. The tree is <b>not</b> updated until you call init(). */
    void addTriangle(Triangle const & tri)
    {
      tris.addTriangle(tri);
    }

    /**
     * Add a set of triangles to the BVH. TriangleIterator must dereference to Triangle. The tree is <b>not</b> updated
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
     * Compute the BVH from the added meshes. You <b>must</b> call this function to construct (or recompute) the tree after
     * any addMesh() or addMeshGroup() calls. Clears the triangle cache.
     */
    void init(typename BaseT::Method method = BaseT::Method::AUTO, int max_depth = -1, int max_elems_in_leaf = -1,
              bool save_memory = false, bool deallocate_previous_memory = true)
    {
      TriangleArray const & tri_array = tris.getTriangles();
      BaseT::init(tri_array.begin(), tri_array.end(), method,
                  max_depth, max_elems_in_leaf, save_memory, deallocate_previous_memory);
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

    void updateElements(typename BaseT::Node const * start = nullptr)
    {
      // Recompute triangle cached properties in case vertex positions have changed
      if (!start || start == this->getRoot())  // non-recursive loop avoids double-processing elements spanning multiple leaves
      {
        auto const * all_elems = this->getElements();
        for (auto ei = all_elems, all_elems_end = all_elems + BaseT::numElements(); ei != all_elems_end; ++ei)
          ei->update();
      }
      else
        updateElementsRecursive(start);
    }

  protected:
    void samplePointsFromElements(intx num_samples, typename BaseT::ElementSample * samples) const
    {
      for (intx i = 0; i < num_samples; ++i)
        samples[i].position = samples[i].element->randomPoint();
    }

  private:
    /** Recursive helper function for updateElements(). */
    void updateElementsRecursive(typename BaseT::Node const * start)
    {
      if (start->isLeaf())
      {
        auto const * all_elems = BaseT::getElements();
        for (auto ei = start->elementIndicesBegin(), ei_end = start->elementIndicesEnd(); ei != ei_end; ++ei)
          all_elems[*ei].update();
      }
      else
      {
        for (intx i = 0; i < start->maxDegree(); ++i)
        {
          auto c = start->getChild(i);
          if (c) { updateElementsRecursive(c); }
        }
      }
    }

    Triangles tris;  ///< Internal cache of triangles used to initialize the tree.

}; // class MeshBvh

} // namespace Algorithms
} // namespace Thea

#endif
