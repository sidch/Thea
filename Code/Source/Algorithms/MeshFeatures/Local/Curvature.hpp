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
// First version: 2013
//
//============================================================================

#ifndef __Thea_Algorithms_MeshFeatures_Local_Curvature_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_Curvature_hpp__

#include "../../../Common.hpp"
#include "../SampledSurface.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Noncopyable.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/** Compute the curvature at a point on a shape. */
template < typename ExternalSampleKdTreeT = KdTreeN<MeshFeatures::SurfaceSample, 3> >
class Curvature : public SampledSurface<ExternalSampleKdTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKdTreeT> BaseT;  ///< Base class.
    static intx const DEFAULT_NUM_SAMPLES = 100000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute curvature at sample points on a given mesh. Initializes internal data structures that do
     * not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    Curvature(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute curvature at sample points on a given mesh group. Initializes internal data structures
     * that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    Curvature(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute curvature of a shape with a precomputed kd-tree on points densely sampled from the
     * shape. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    Curvature(ExternalSampleKdTreeT const * sample_kdtree_, Real normalization_scale = -1)
    : BaseT(sample_kdtree_, normalization_scale)
    {}

    /**
     * Compute the <em>projected</em> curvature at a query point on the mesh. The projected curvature is an approximation to the
     * actual curvature, obtained by projecting sample points in the neighborhood of the query point onto the normal.
     *
     * This version of the function explicitly computes the normal at the query point -- the other version of the function
     * should be used if the normal is known in advance. The normal is computed from samples, so <b>may be quite inaccurate</b>
     * especially in thin areas.
     *
     * @param position The position at which to compute curvature.
     * @param nbd_radius The size of the neighborhood over which to compute curvature, specified as a multiple of the shape
     *   scale. A negative argument selects a default size.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    double computeProjectedCurvature(Vector3 const & position, Real nbd_radius = -1) const
    {
      intx nn_index = this->hasExternalKdTree() ? this->getExternalKdTree()->template closestElement<MetricL2>(position)
                                                : this->getInternalKdTree()->template closestElement<MetricL2>(position);
      if (nn_index < 0)
      {
        THEA_WARNING << "Curvature: Query point cannot be mapped to mesh, curvature value set to zero";
        return 0.0;
      }

      return computeProjectedCurvature(position, this->getSampleNormal(nn_index), nbd_radius);
    }

    /**
     * Compute the <em>projected</em> curvature at a query point with a known normal on the mesh. The projected curvature is an
     * approximation to the actual curvature, obtained by projecting sample points in the neighborhood of the query point onto
     * the normal.
     *
     * @param position The position at which to compute curvature.
     * @param normal The normal at this position.
     * @param nbd_radius The size of the neighborhood over which to compute curvature, specified as a multiple of the shape
     *   scale. A negative argument selects a default size.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    double computeProjectedCurvature(Vector3 const & position, Vector3 const & normal, Real nbd_radius = -1) const
    {
      if (nbd_radius <= 0)
        nbd_radius = 0.1f;

      nbd_radius *= this->getNormalizationScale();

      ProjectedCurvatureFunctor func(position, normal);
      Ball3 range(position, nbd_radius);

      if (this->hasExternalKdTree())
        this->getMutableExternalKdTree()->template processRangeUntil<IntersectionTester>(range, std::ref(func));
      else
        this->getMutableInternalKdTree()->template processRangeUntil<IntersectionTester>(range, std::ref(func));

      return func.getCurvature();
    }

  private:
    /** Called for each point in the neighborhood. */
    struct ProjectedCurvatureFunctor : public Noncopyable
    {
      ProjectedCurvatureFunctor(Vector3 const & p, Vector3 const & n)
      : position(p), normal(n), num_offsets(0), sum_offsets(Vector3::Zero())
      {}

      template <typename SampleT> bool operator()(intx index, SampleT & t)
      {
        if (NormalTraits<SampleT>::getNormal(t).dot(normal) > -1.0e-05f)  // ignore points on hidden side
        {
          Vector3 offset = (PointTraitsN<SampleT, 3>::getPosition(t) - position).normalized();
          sum_offsets += offset;
          num_offsets++;
        }

        return false;
      }

      double getCurvature() const
      {
        return num_offsets > 0 ? -sum_offsets.dot(normal) / num_offsets : 0;
      }

      Vector3 position, normal;
      intx num_offsets;
      Vector3 sum_offsets;

    }; // struct ProjectedCurvatureFunctor

}; // class Curvature

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
