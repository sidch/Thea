//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_MeshSampler_hpp__
#define __Thea_Algorithms_MeshSampler_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "MeshKDTree.hpp"
#include "../Random.hpp"
#include "../Vector3.hpp"

namespace Thea {
namespace Algorithms {

/** Sample points from a mesh. */
template <typename MeshT>
class MeshSampler
{
  public:
    typedef MeshT                      Mesh;      ///< The mesh class.
    typedef MeshKDTree<Mesh>           KDTree;    ///< A kd-tree on the mesh.
    typedef typename KDTree::Triangle  Triangle;  ///< A kd-tree triangle belonging to a face of the mesh.

    /**
     * Constructs the object to compute sample points on a given mesh. The mesh must persist as long as this object does.
     * Initializes internal data structures that do not need to be recomputed for successive calls to functions that generate
     * samples on this mesh.
     */
    MeshSampler(Mesh const & mesh) : kdtree(new KDTree), owns_kdtree(true)
    {
      kdtree->add(const_cast<Mesh &>(mesh));  // safe -- the kd-tree won't be used to modify the mesh
      // don't init() the tree, we only need access to the triangle array;
    }

    /**
     * Constructs the object to compute sample points on a given mesh group. The mesh group must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to functions that
     * generate samples on this mesh.
     */
    MeshSampler(Graphics::MeshGroup<Mesh> const & mesh_group) : kdtree(new KDTree), owns_kdtree(true)
    {
      kdtree->add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));  // safe -- the kd-tree won't be used to modify the mesh
      // don't init() the tree, we only need access to the triangle array;
    }

    /**
     * Constructs the object to compute sample points on a mesh with a precomputed kd-tree. The kd-tree must persist as long as
     * this object does.
     */
    MeshSampler(KDTree const * kdtree_) : kdtree(const_cast<KDTree *>(kdtree_)), owns_kdtree(false)
    {
      alwaysAssertM(kdtree, "Curvature: Precomputed KD-tree cannot be null");
    }

    /** Destructor. */
    ~MeshSampler()
    {
      if (owns_kdtree)
        delete kdtree;
    }

    /**
     * Samples points from the mesh evenly by area.
     *
     * @param approx_num_samples The approximate number of samples to compute. This need <em>not</em> be exactly the number
     *   actually computed (and returned by the function).
     * @param positions Used to return the sample positions.
     * @param face_normals If not null, used to return the normals of the parent faces of the samples.
     * @param triangles If not null, used to return the kd-tree triangles from which the samples were selected.
     *
     * @return The number of samples computed.
     */
    long sampleEvenlyByArea(long approx_num_samples,
                            TheaArray<Vector3> & positions,
                            TheaArray<Vector3> * face_normals = NULL,
                            TheaArray<Triangle const *> * triangles = NULL) const
    {
      typedef typename KDTree::Triangle Triangle;
      typedef typename KDTree::TriangleArray TriangleArray;

      positions.clear();
      if (face_normals) face_normals->clear();
      if (triangles) triangles->clear();

      TriangleArray const & all_tris = kdtree->getTriangles();  // computed even without call to kdtree->init()
      if (all_tris.empty())
        return 0;

      double total_area = 0;
      for (array_size_t i = 0; i < all_tris.size(); ++i)
        total_area += all_tris[i].getArea();

      if (total_area <= 0)
        return 0;

      for (array_size_t i = 0; i < all_tris.size(); ++i)
      {
        double frac_samples = (approx_num_samples * all_tris[i].getArea()) / total_area;
        while (frac_samples > 0)
        {
          // Last (possible) sample?
          if (frac_samples < 1 && Random::common().uniform01() > frac_samples)
            break;

          positions.push_back(all_tris[i].randomPoint());
          if (face_normals) face_normals->push_back(all_tris[i].getNormal());
          if (triangles) triangles->push_back(&all_tris[i]);

          frac_samples -= 1.0;
        }
      }

      return (long)positions.size();
    }

  private:
    KDTree * kdtree;  ///< KD-tree on the mesh.
    bool owns_kdtree;  ///< Did this object create the kd-tree, or does it simply wrap a precomputed one?

}; // class MeshSampler

} // namespace Algorithms
} // namespace Thea

#endif
