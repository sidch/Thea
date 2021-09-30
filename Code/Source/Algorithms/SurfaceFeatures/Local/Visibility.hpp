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
// First version: 2019
//
//============================================================================

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_Visibility_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_Visibility_hpp__

#include "../../../Common.hpp"
#include "../../../Graphics/MeshGroup.hpp"
#include "../../MeshBvh.hpp"
#include "../../RayIntersectionTester.hpp"
#include "../../../Math.hpp"
#include "../../../MatVec.hpp"
#include "../../../Random.hpp"
#include "../../../Ray3.hpp"

namespace Thea {
namespace Algorithms {

/** Namespace for classes that compute mesh features. */
namespace SurfaceFeatures {

/** Namespace for classes that compute local features on a mesh. */
namespace Local {

/** Compute the external visibility of a point on a mesh, by shooting rays outwards from it and checking how many escape. */
template < typename MeshT, typename ExternalBvhT = MeshBvh<MeshT> >
class Visibility
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.
    typedef ExternalBvhT ExternalBvh;  ///< A precomputed BVH on the mesh.

  private:
    typedef MeshBvh<Mesh> Bvh;  ///< A BVH on the mesh.

  public:
    /**
     * Constructs the object to evaluate visibility of points on a given mesh. The mesh must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     */
    Visibility(Mesh const & mesh)
    : bvh(new Bvh), precomp_bvh(nullptr)
    {
      bvh->add(const_cast<Mesh &>(mesh));  // safe -- the BVH won't be used to modify the mesh
      bvh->init();
      scale = bvh->getBounds().getExtent().norm();
    }

    /**
     * Constructs the object to evaluate visibility of points on a given mesh. The mesh must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh_group The mesh group representing the shape.
     */
    Visibility(Graphics::MeshGroup<Mesh> const & mesh_group)
    : bvh(new Bvh), precomp_bvh(nullptr)
    {
      bvh->add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));  // safe -- the BVH won't be used to modify the mesh
      bvh->init();
      scale = bvh->getBounds().getExtent().norm();
    }

    /**
     * Constructs the object to evaluate visibility of points on a shape with a precomputed BVH. The BVH must persist as
     * long as this object does.
     *
     * @param bvh_ The precomputed BVH representing the shape.
     */
    Visibility(ExternalBvh const * bvh_)
    : bvh(nullptr), precomp_bvh(bvh_)
    {
      alwaysAssertM(precomp_bvh, "Visibility: Precomputed BVH cannot be null");

      scale = precomp_bvh->getBounds().getExtent().norm();
    }

    /** Destructor. */
    ~Visibility()
    {
      delete bvh;
    }

    /**
     * Compute the external visibility of a query point on the shape.
     *
     * @param position The position of the query point.
     * @param num_rays The number of rays to sample randomly from a sphere (negative for default).
     *
     * @return The fraction of rays from the query point that do not intersect the shape again.
     */
    double compute(Vector3 const & position, intx num_rays = -1) const
    {
      static intx const DEFAULT_NUM_RAYS = 100;

      if (num_rays <= 0)  // 0 is probably user error, snap it to default as well though this should not be relied upon
        num_rays = DEFAULT_NUM_RAYS;

      Vector3 u;
      intx num_escaped = 0;
      for (intx i = 0; i < num_rays; ++i)
      {
        Random::common().sphere(u[0], u[1], u[2]);
        Vector3 offset = 0.001f * scale * u;
        Ray3 ray = Ray3(position + offset, u);

        bool hit = precomp_bvh ? precomp_bvh->template rayIntersects<RayIntersectionTester>(ray)
                               : bvh->template rayIntersects<RayIntersectionTester>(ray);
        if (!hit)
          num_escaped++;
      }

      return (double)num_escaped / num_rays;
    }

  private:
    Bvh * bvh;  ///< Self-owned BVH on the mesh for computing ray intersections.
    ExternalBvh const * precomp_bvh;  ///< Precomputed BVH on the mesh for computing ray intersections.
    Real scale;  ///< The normalization scale for offsetting ray origins.

}; // class Visibility

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
