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
// First version: 2015
//
//============================================================================

#ifndef __Thea_Algorithms_MeshFeatures_Local_SpinImage_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_SpinImage_hpp__

#include "../../../Common.hpp"
#include "../SampledSurface.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Math.hpp"
#include "../../../MatVec.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/**
 * Compute the spin image at a point on a shape. The spin image at a point is a histogram of the rest of the surface, where each
 * bin corresponds to an annular ring with the normal at the central point as its axis. The bins are parametrized by radial
 * distance from this axis, and height above/below the tangent plane.
 *
 * A. Johnson and M. Hebert, "Using Spin Images for Efficient Object Recognition in Cluttered 3D Scenes",
 * IEEE Trans. PAMI 21(5), 433-449, 1999.
 */
template < typename ExternalSampleKdTreeT = KdTreeN<MeshFeatures::SurfaceSample, 3> >
class SpinImage : public SampledSurface<ExternalSampleKdTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKdTreeT> BaseT;  ///< Base class.
    static intx const DEFAULT_NUM_SAMPLES = 50000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute spin images of a shape with a precomputed set of surface samples.
     *
     * @param num_samples The number of precomputed samples.
     * @param positions The positions of the samples.
     * @param normals The normals of the samples.
     * @param normalization_scale The scale of the shape, used to define the extents of the spin image. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    SpinImage(intx num_samples, Vector3 const * positions, Vector3 const & normals, Real normalization_scale = -1)
    : BaseT(num_samples, positions, normals, normalization_scale)
    {}

    /**
     * Constructs the object to compute spin images at sample points on a given mesh. Initializes internal data structures that
     * do not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the extents of the spin image. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    template <typename MeshT>
    SpinImage(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute spin images at sample points on a given mesh group. Initializes internal data structures
     * that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the extents of the spin image. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    template <typename MeshT>
    SpinImage(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute spin images of a shape with a precomputed kd-tree on points densely sampled from the
     * shape. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define the extents of the spin image. If <= 0, the bounding
     *   sphere diameter will be used.
     */
    SpinImage(ExternalSampleKdTreeT const * sample_kdtree_, Real normalization_scale = -1)
    : BaseT(sample_kdtree_, normalization_scale)
    {}

    /**
     * Compute the spin image at a query point on the mesh.
     *
     * This version of the function explicitly computes the normal at the query point -- the other version of the function
     * should be used if the normal is known in advance. The normal is computed from samples, so <b>may be quite inaccurate</b>
     * especially in thin areas.
     *
     * @param position The point at which the spin image is to be calculated.
     * @param num_radial_bins The number of bins in the radial direction (distance from axis).
     * @param num_height_bins The number of bins in the height direction (distance from tangent plane).
     * @param spin_image Used to return the computed spin image, with rows corresponding to radial divisions and columns to
     *   height divisions.
     */
    void compute(Vector3 const & position, int num_radial_bins, int num_height_bins, MatrixX<double> & spin_image) const
    {
      intx nn_index = this->hasExternalKdTree() ? this->getExternalKdTree()->template closestElement<MetricL2>(position)
                                                : this->getInternalKdTree()->template closestElement<MetricL2>(position);
      if (nn_index < 0)
      {
        THEA_WARNING << "SpinImage: Query point cannot be mapped to mesh, spin image set to zero";
        spin_image.resize(num_radial_bins, num_height_bins);
        spin_image.setZero();
        return;
      }

      compute(position, this->getSampleNormal(nn_index), num_radial_bins, num_height_bins, spin_image);
    }

    /**
     * Compute the spin image at a query point on the mesh with a known normal.
     *
     * @param position The point at which the spin image is to be calculated.
     * @param normal The normal at this point.
     * @param num_radial_bins The number of bins in the radial direction (distance from axis).
     * @param num_height_bins The number of bins in the height direction (distance from tangent plane).
     * @param spin_image Used to return the computed spin image, with rows corresponding to radial divisions and columns to
     *   height divisions.
     */
    void compute(Vector3 const & position, Vector3 const & normal, int num_radial_bins, int num_height_bins,
                 MatrixX<double> & spin_image) const
    {
      alwaysAssertM(num_radial_bins > 0, "SpinImage: Number of radial bins must be positive");
      alwaysAssertM(num_height_bins > 0, "SpinImage: Number of height bins must be positive");

      // Guess suitable limiting extents of the spin image, to have enough non-zero bins while not having too much of the
      // surface outside these extents
      Real max_radius = 0.75f * this->getNormalizationScale();  // radial limits: [0, max_radius]
      Real max_height = 0.75f * this->getNormalizationScale();  // height limits: [-max_height, max_height]

      spin_image.resize(num_radial_bins, num_height_bins);
      spin_image.setZero();

      Vector3 axis = normal.normalized();
      intx n = this->numSamples();
      for (intx i = 0; i < n; ++i)
      {
        Vector3 offset = this->getSamplePosition(i) - position;
        Real height = offset.dot(axis);
        Real radius = (offset - height * axis).norm();

        int radial_bin = (int)Math::clamp((int)std::floor(num_radial_bins * (radius / max_radius)), 0, num_radial_bins - 1);
        int height_bin = (int)Math::clamp((int)std::floor(num_height_bins * (0.5f * (height / max_height + 1))),
                                                          0, num_height_bins - 1);

        spin_image(radial_bin, height_bin) += 1;
      }
    }

}; // class SpinImage

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
