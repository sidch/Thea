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
template <typename MeshT>
class ShapeDiameter
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.

  private:
    typedef MeshKDTree<Mesh> KDTree;  ///< Ray-intersection search structure.

  public:
    /**
     * Constructs the object to compute shape diameter at sample points on a given mesh. The mesh must persist as long as this
     * object does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     */
    ShapeDiameter(Mesh const & mesh) : scale(0)
    {
      kdtree.add(mesh);
      kdtree.init();
      scale = kdtree.getBounds().getExtent().length();
    }

    /**
     * Constructs the object to compute shape diameter at sample points on a given mesh group. The mesh group must persist as
     * long as this object does. Initializes internal data structures that do not need to be recomputed for successive calls to
     * compute().
     */
    ShapeDiameter(Graphics::MeshGroup<Mesh> const & mesh_group) : scale(0)
    {
      kdtree.add(mesh_group);
      kdtree.init();
      scale = kdtree.getBounds().getExtent().length();
    }

    /** Destructor. */
    ~ShapeDiameter() { delete kdtree; }

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
        long nn_index = kdtree.closestElement<MetricL2>(positions[i]);
        if (nn_index < 0)
        {
          THEA_WARNING << "ShapeDiameter: Query point cannot be mapped to mesh, all SDF values set to zero";
          sdf_values.resize(positions.size());
          std::fill(sdf_values.begin(), sdf_values.end(), 0);
          return;
        }

        normals[i] = kdtree.getElements()[(array_size_t)nn_index].getNormal();
      }

      compute(positions, normals, sdf_values);
    }

    /** Compute the shape diameter function at a given set of sample points with known normals on the mesh. */
    void compute(TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals, TheaArray<Real> & sdf_values) const
    {
      alwaysAssertM(positions.size() == normals.size(), "ShapeDiameter: Number of sample positions and normals do not match");
      alwaysAssertM(kdtree, "ShapeDiameter: KD-tree not initialized for mesh");

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
        Vector3( 0.407497,  0.134086,  0.903309),
        Vector3( 0.128722, -0.822715,  0.553688),
        Vector3( 0.359298,  0.686129,  0.632560),
        Vector3( 0.813318, -0.153956,  0.561081),
        Vector3(-0.707169,  0.296987,  0.641647),
        Vector3(-0.478783, -0.543441,  0.689520),
        Vector3( 0.561083,  0.523634,  0.641087),
        Vector3(-0.380898, -0.655789,  0.651811),
        Vector3( 0.722588,  0.472155,  0.504912),
        Vector3(-0.535205,  0.421235,  0.732200),
        Vector3( 0.577593, -0.477896,  0.661817),
        Vector3( 0.393652, -0.144484,  0.907834),
        Vector3( 0.141508, -0.660777,  0.737122),
        Vector3(-0.274967, -0.141703,  0.950954),
        Vector3(-0.304935,  0.474085,  0.825989),
        Vector3( 0.711150,  0.301383,  0.635164),
        Vector3(-0.696119, -0.112341,  0.709082),
        Vector3(-0.248327,  0.381035,  0.890588),
        Vector3(-0.0772394, 0.0462501, 0.995939),
        Vector3( 0.438024,  0.656463,  0.614159),
        Vector3(-0.483384, -0.542477,  0.687065),
        Vector3( 0.498569, -0.697466,  0.514752),
        Vector3(-0.128902, -0.0173429, 0.991506),
        Vector3( 0.338759, -0.318772,  0.885227),
        Vector3(-0.574364, -0.463887,  0.674474),
        Vector3(-0.459207,  0.260125,  0.849390),
        Vector3(-0.628057, -0.528928,  0.570771),
        Vector3( 0.232916, -0.328273,  0.915416),
        Vector3( 0.231667, -0.796882,  0.557951),
        Vector3(-0.270612, -0.809654,  0.520797),
      };

      Real values[NUM_RAYS];
      Real weights[NUM_RAYS];
      int num_values = 0;
      for (int i = 0; i < NUM_RAYS; ++i)
      {
        Vector3 dir = rot * CONE_DIRS[i];
        Ray3 ray(p + offset, dir);
        RayStructureIntersection3 isec = kdtree->rayStructureIntersection<RayIntersectionTester>(ray);

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
    void
    getTwoMutuallyPerpendicularAxes(Vector3 const & dir, Vector3 & u, Vector3 & v)
    {
      u = (v.maxAbsAxis() == 0) ? Vector3(v.y(), -v.x(), 0) : Vector3(0, v.z(), -v.y());
      u.unitize();
      v = dir.cross(u).unit();
    }

    KDTree kdtree;  ///< KD-tree on the mesh for computing ray intersections.
    Real scale;  ///< The normalization length.

}; // class ShapeDiameter

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
