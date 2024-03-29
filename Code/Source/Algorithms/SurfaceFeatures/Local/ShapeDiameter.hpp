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
// First version: 2012
//
//============================================================================

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_ShapeDiameter_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_ShapeDiameter_hpp__

#include "../../../Common.hpp"
#include "../../../Graphics/MeshGroup.hpp"
#include "../../BestFitSphere3.hpp"
#include "../../MeshBvh.hpp"
#include "../../MetricL2.hpp"
#include "../../PointCollectorN.hpp"
#include "../../RayIntersectionTester.hpp"
#include "../../../Math.hpp"
#include "../../../MatVec.hpp"
#include <algorithm>

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {

/** Namespace for classes that compute local features on a surface. */
namespace Local {

/**
 * Compute the shape diameter function (SDF) at a point on a surface.
 *
 * Gal, Shamir and Cohen-Or, "Pose-Oblivious Shape Signature", IEEE TVCG 2007.
 */
template < typename MeshT, typename ExternalBvhT = MeshBvh<MeshT> >
class ShapeDiameter
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.
    typedef ExternalBvhT ExternalBvh;  ///< A precomputed BVH on the mesh.

  private:
    typedef MeshBvh<Mesh> Bvh;  ///< A BVH on the mesh.

  public:
    /**
     * Constructs the object to compute shape diameter at sample points on a given mesh. The mesh must persist as long as this
     * object does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     * @param normalization_scale The scale of the shape, used to normalize shape diameters to [0, 1]. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    ShapeDiameter(Mesh const & mesh, Real normalization_scale = -1)
    : bvh(new Bvh), precomp_bvh(nullptr), scale(normalization_scale)
    {
      bvh->add(const_cast<Mesh &>(mesh));  // safe -- the BVH won't be used to modify the mesh
      bvh->init();

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Constructs the object to compute the shape diameter at sample points on a given mesh group. The mesh group must persist
     * as long as this object does. Initializes internal data structures that do not need to be recomputed for successive calls
     * to compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param normalization_scale The scale of the shape, used to normalize shape diameters to [0, 1]. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    ShapeDiameter(Graphics::MeshGroup<Mesh> const & mesh_group, Real normalization_scale = -1)
    : bvh(new Bvh), precomp_bvh(nullptr), scale(normalization_scale)
    {
      bvh->add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));  // safe -- the BVH won't be used to modify the mesh
      bvh->init();

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh_group);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Constructs the object to compute the shape diameter at sample points of a shape with a precomputed BVH. The BVH
     * must persist as long as this object does.
     *
     * @param bvh_ The precomputed BVH representing the shape.
     * @param normalization_scale The scale of the shape, used to normalize shape diameters to [0, 1]. If <= 0, the bounding
     *   box diagonal will be used.
     *
     * @warning This function uses the <b>bounding box diagonal</b> as the default normalization scale, instead of the bounding
     *   sphere diameter as in the other constructors. This is because the latter cannot currently be computed from only a
     *   BVH. If you want to use the bounding sphere diameter (or other value) as the normalization scale, you must compute
     *   it separately and pass it as a parameter to this function.
     */
    ShapeDiameter(ExternalBvh const * bvh_, Real normalization_scale = -1)
    : bvh(nullptr), precomp_bvh(bvh_), scale(normalization_scale)
    {
      alwaysAssertM(precomp_bvh, "ShapeDiameter: Precomputed BVH cannot be null");

      if (scale <= 0)
        scale = precomp_bvh->getBounds().getExtent().norm();
    }

    /** Destructor. */
    ~ShapeDiameter()
    {
      delete bvh;
    }

    /**
     * Get the normalization scale. This is the length by which ray intersection distances will be divided to get the normalized
     * shape diameter values.
     */
    Real getScale() const { return scale; }

    /**
     * Compute the shape diameter function at a query point on the mesh. This explicitly computes the normal at the sample point
     * point -- the other version of the function should be used if the normal is known in advance. The shape diameter will be
     * normalized to [0, 1] by dividing by the mesh scale, as returned by getScale(). If absolutely no query ray intersects the
     * object, a negative value is returned.
     *
     * @param position The position of the query point.
     * @param only_hit_interior_surfaces Only consider ray intersections with surfaces whose normals are in the same direction
     *   as the ray.
     */
    double compute(Vector3 const & position, bool only_hit_interior_surfaces = true) const
    {
      intx nn_index = precomp_bvh ? precomp_bvh->template closestElement<MetricL2>(position)
                                  : bvh->template closestElement<MetricL2>(position);
      if (nn_index < 0)
      {
        THEA_WARNING << "ShapeDiameter: Query point cannot be mapped to mesh, returning negative SDF value";
        return -1.0;
      }

      // Use the face normal and not the smooth normal, to handle sharp edged slabs etc
      Vector3 normal = precomp_bvh ? precomp_bvh->getElements()[(size_t)nn_index].getNormal()
                                   : bvh->getElements()[(size_t)nn_index].getNormal();

      return compute(position, normal, only_hit_interior_surfaces);
    }

    /**
     * Compute the shape diameter function at a query point with a known (outwards-pointing) normal on the mesh. The shape
     * diameter will be normalized to [0, 1] by dividing by the mesh scale, as returned by getScale(). If absolutely no query
     * ray intersects the object, a negative value is returned.
     *
     * @param position The position of the query point.
     * @param normal The normal of the query point.
     * @param only_hit_interior_surfaces Only consider ray intersections with surfaces whose normals are in the same direction
     *   as the ray.
     */
    double compute(Vector3 const & position, Vector3 const & normal, bool only_hit_interior_surfaces = true) const
    {
      Vector3 in = -normal.normalized();
      Matrix3 rot = Math::orthonormalBasis(in);
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

      double values[NUM_RAYS];
      double weights[NUM_RAYS];
      int num_values = 0;
      for (int i = 0; i < NUM_RAYS; ++i)
      {
        Vector3 dir = rot * CONE_DIRS[i];
        Ray3 ray(position + offset, dir);
        RayStructureIntersection3 isec = precomp_bvh
                                       ? precomp_bvh->template rayStructureIntersection<RayIntersectionTester>(ray)
                                       : bvh->template rayStructureIntersection<RayIntersectionTester>(ray);

        if (isec.isValid() && (!only_hit_interior_surfaces || isec.getNormal().dot(dir) >= 0))
        {
          values[num_values] = isec.getTime();
          weights[num_values] = CONE_DIRS[i][2];  // cos(angle) is just the z-component
          num_values++;
        }
      }

      if (num_values <= 0)
        return -1.0;  // either the normal is in the wrong direction or this is a 2D surface

      // Outlier rejection: reject all values more than one standard deviation from the median
      int mid = num_values / 2;  // integer division takes floor
      std::nth_element(values, values + mid, values + num_values);

      // Compute variance
      double sum_values = 0, sum_squares = 0;
      for (int i = 0; i < num_values; ++i)
      {
        sum_values += values[i];
        sum_squares += (values[i] * values[i]);
      }

      double avg = sum_values / num_values;
      double var = sum_squares / num_values - avg * avg;

      sum_values = 0;
      double sum_weights = 0;
      for (int i = 0; i < num_values; ++i)
        if (Math::square(values[i] - values[mid]) <= var)
        {
          sum_values += (weights[i] * values[i]);
          sum_weights += weights[i];
        }

      if (sum_weights > 0)
        return Math::clamp((double)(sum_values / (sum_weights * scale)), 0.0, 1.0);
      else  // should never happen
        return Math::clamp((double)(values[mid] / scale), 0.0, 1.0);
    }

  private:
    Bvh * bvh;  ///< Self-owned BVH on the mesh for computing ray intersections.
    ExternalBvh const * precomp_bvh;  ///< Precomputed BVH on the mesh for computing ray intersections.
    Real scale;  ///< The normalization length.

}; // class ShapeDiameter

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
