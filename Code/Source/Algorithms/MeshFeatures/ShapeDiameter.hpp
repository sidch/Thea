//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_MeshFeatures_ShapeDiameter_hpp__
#define __Thea_Algorithms_MeshFeatures_ShapeDiameter_hpp__

#include "../../Common.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../MeshKDTree.hpp"
#include "../MetricL2.hpp"
#include "../RayIntersectionTester.hpp"
#include "../../Math.hpp"
#include "../../Matrix3.hpp"
#include "../../Vector3.hpp"
#include <algorithm>

namespace Thea {
namespace Algorithms {

/** Namespace for classes that compute features on a mesh. */
namespace MeshFeatures {

/**
 * Computes the shape diameter function (SDF) at a given set of positions on a mesh.
 *
 * Gal, Shamir and Cohen-Or, "Pose-Oblivious Shape Signature", IEEE TVCG 2007.
 */
template < typename MeshT, typename ExternalKDTreeT = MeshKDTree<MeshT> >
class ShapeDiameter
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.
    typedef ExternalKDTreeT ExternalKDTree;  ///< A precomputed kd-tree on the mesh.

  private:
    typedef MeshKDTree<Mesh> KDTree;  ///< A kd-tree on the mesh.

  public:
    /**
     * Constructs the object to compute shape diameter at sample points on a given mesh. The mesh must persist as long as this
     * object does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     */
    ShapeDiameter(Mesh const & mesh) : kdtree(new KDTree), precomp_kdtree(NULL), scale(0)
    {
      kdtree->add(const_cast<Mesh &>(mesh));  // safe -- the kd-tree won't be used to modify the mesh
      kdtree->init();
      scale = kdtree->getBounds().getExtent().length();
    }

    /**
     * Constructs the object to compute the shape diameter at sample points on a given mesh group. The mesh group must persist
     * as long as this object does. Initializes internal data structures that do not need to be recomputed for successive calls
     * to compute().
     */
    ShapeDiameter(Graphics::MeshGroup<Mesh> const & mesh_group) : kdtree(new KDTree), precomp_kdtree(NULL), scale(0)
    {
      kdtree->add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));  // safe -- the kd-tree won't be used to modify the mesh
      kdtree->init();
      scale = kdtree->getBounds().getExtent().length();
    }

    /**
     * Constructs the object to compute the shape diameter at sample points of a shape with a precomputed kd-tree. The kd-tree
     * must persist as long as this object does.
     */
    ShapeDiameter(ExternalKDTree const * kdtree_) : kdtree(NULL), precomp_kdtree(kdtree_), scale(0)
    {
      alwaysAssertM(precomp_kdtree, "ShapeDiameter: Precomputed KD-tree cannot be null");
      scale = precomp_kdtree->getBounds().getExtent().length();
    }

    /** Destructor. */
    ~ShapeDiameter()
    {
      delete kdtree;
    }

    /**
     * Get the normalization scale. This is the length by which ray intersection distances will be divided to get the normalized
     * shape diameter values.
     */
    Real getNormalizationScale() const { return scale; }

    /**
     * Compute the shape diameter function at a given set of sample points on the mesh. This explicitly computes the normal at
     * each point -- the other version of the function should be used if the normals are known in advance.
     */
    void compute(TheaArray<Vector3> const & positions, TheaArray<Real> & sdf_values) const
    {
      TheaArray<Vector3> const & normals(positions.size());
      for (array_size_t i = 0; i < positions.size(); ++i)
      {
        long nn_index = precomp_kdtree ? precomp_kdtree->template closestElement<MetricL2>(positions[i])
                                       : kdtree->template closestElement<MetricL2>(positions[i]);
        if (nn_index < 0)
        {
          THEA_WARNING << "ShapeDiameter: Query point cannot be mapped to mesh, all SDF values set to zero";
          sdf_values.resize(positions.size());
          std::fill(sdf_values.begin(), sdf_values.end(), 0);
          return;
        }

        normals[i] = precomp_kdtree ? precomp_kdtree->getElements()[(array_size_t)nn_index].getNormal()
                                    : kdtree->getElements()[(array_size_t)nn_index].getNormal();
      }

      compute(positions, normals, sdf_values);
    }

    /** Compute the shape diameter function at a given set of sample points with known normals on the mesh. */
    void compute(TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals, TheaArray<Real> & sdf_values) const
    {
      alwaysAssertM(positions.size() == normals.size(), "ShapeDiameter: Number of sample positions and normals do not match");
      alwaysAssertM(kdtree || precomp_kdtree, "ShapeDiameter: KD-tree not initialized for mesh");

      sdf_values.resize(positions.size());

      for (array_size_t i = 0; i < positions.size(); ++i)
        sdf_values[i] = computeSDF(positions[i], normals[i]);
    }

  private:
    /** Compute the SDF at a point. */
    Real computeSDF(Vector3 const & p, Vector3 const & n) const
    {
      Vector3 in = -n;
      Vector3 u, v;
      getTwoMutuallyPerpendicularAxes(in, u, v);
      Matrix3 rot(u[0], v[0], in[0],
                  u[1], v[1], in[1],
                  u[2], v[2], in[2]);

      Vector3 offset = 0.001f * scale * in;

      static int const NUM_RAYS = 30;
      static Vector3 const CONE_DIRS[NUM_RAYS] = {
        Vector3( 0.407497f,  0.134086f,  0.903309f),
        Vector3( 0.128722f, -0.822715f,  0.553688f),
        Vector3( 0.359298f,  0.686129f,  0.632560f),
        Vector3( 0.813318f, -0.153956f,  0.561081f),
        Vector3(-0.707169f,  0.296987f,  0.641647f),
        Vector3(-0.478783f, -0.543441f,  0.689520f),
        Vector3( 0.561083f,  0.523634f,  0.641087f),
        Vector3(-0.380898f, -0.655789f,  0.651811f),
        Vector3( 0.722588f,  0.472155f,  0.504912f),
        Vector3(-0.535205f,  0.421235f,  0.732200f),
        Vector3( 0.577593f, -0.477896f,  0.661817f),
        Vector3( 0.393652f, -0.144484f,  0.907834f),
        Vector3( 0.141508f, -0.660777f,  0.737122f),
        Vector3(-0.274967f, -0.141703f,  0.950954f),
        Vector3(-0.304935f,  0.474085f,  0.825989f),
        Vector3( 0.711150f,  0.301383f,  0.635164f),
        Vector3(-0.696119f, -0.112341f,  0.709082f),
        Vector3(-0.248327f,  0.381035f,  0.890588f),
        Vector3(-0.0772394f, 0.0462501f, 0.995939f),
        Vector3( 0.438024f,  0.656463f,  0.614159f),
        Vector3(-0.483384f, -0.542477f,  0.687065f),
        Vector3( 0.498569f, -0.697466f,  0.514752f),
        Vector3(-0.128902f, -0.0173429f, 0.991506f),
        Vector3( 0.338759f, -0.318772f,  0.885227f),
        Vector3(-0.574364f, -0.463887f,  0.674474f),
        Vector3(-0.459207f,  0.260125f,  0.849390f),
        Vector3(-0.628057f, -0.528928f,  0.570771f),
        Vector3( 0.232916f, -0.328273f,  0.915416f),
        Vector3( 0.231667f, -0.796882f,  0.557951f),
        Vector3(-0.270612f, -0.809654f,  0.520797f),
      };

      Real values[NUM_RAYS];
      Real weights[NUM_RAYS];
      int num_values = 0;
      for (int i = 0; i < NUM_RAYS; ++i)
      {
        Vector3 dir = rot * CONE_DIRS[i];
        Ray3 ray(p + offset, dir);
        RayStructureIntersection3 isec = precomp_kdtree
                                       ? precomp_kdtree->template rayStructureIntersection<RayIntersectionTester>(ray)
                                       : kdtree->template rayStructureIntersection<RayIntersectionTester>(ray);

        if (isec.isValid() && isec.getNormal().dot(dir) >= 0)
        {
          values[num_values] = isec.getTime();
          weights[num_values] = CONE_DIRS[i][2];  // cos(angle) is just the z-component
          num_values++;
        }
      }

      if (num_values <= 0)
        return 0;  // assumes flat plate

      // Outlier rejection: reject all values more than one standard deviation from the median
      int mid = num_values / 2;  // integer division takes floor
      std::nth_element(values, values + mid, values + num_values);

      // Compute variance
      Real sum_values = 0, sum_squares = 0;
      for (int i = 0; i < num_values; ++i)
      {
        sum_values += values[i];
        sum_squares += (values[i] * values[i]);
      }

      Real avg = sum_values / num_values;
      Real var = sum_squares / num_values - avg * avg;

      sum_values = 0;
      Real sum_weights = 0;
      for (int i = 0; i < num_values; ++i)
        if (Math::square(values[i] - values[mid]) <= var)
        {
          sum_values += (weights[i] * values[i]);
          sum_weights += weights[i];
        }

      if (sum_weights > 0)
        return Math::clamp(sum_values / (sum_weights * scale), (Real)0, (Real)1);
      else  // should never happen
        return Math::clamp(values[mid] / scale, (Real)0, (Real)1);
    }

    /** Get two unit vectors perpendicular to a given vector and to each other. */
    static void
    getTwoMutuallyPerpendicularAxes(Vector3 const & dir, Vector3 & u, Vector3 & v)
    {
      u = (v.maxAbsAxis() == 0) ? Vector3(v.y(), -v.x(), 0) : Vector3(0, v.z(), -v.y());
      u.unitize();
      v = dir.cross(u).unit();
    }

    KDTree * kdtree;  ///< Self-owned KD-tree on the mesh for computing ray intersections.
    ExternalKDTree const * precomp_kdtree;  ///< Precomputed KD-tree on the mesh for computing ray intersections.
    Real scale;  ///< The normalization length.

}; // class ShapeDiameter

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
