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

#include "SpinImage.hpp"
#include "../../../Math.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

SpinImage::SpinImage(PointSet3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "SpinImage: Cannot construct with a null surface");
}

void
SpinImage::compute(Vector3 const & position, int num_radial_bins, int num_height_bins, MatrixX<double> & spin_image) const
{
  alwaysAssertM(surf->hasNormals(), "SpinImage: Cannot infer query normal from surface lacking normals");

  intx nn_index = surf->getBvh().closestElement<MetricL2>(position);
  if (nn_index < 0)
  {
    THEA_WARNING << "SpinImage: Query point cannot be mapped to mesh, spin image set to zero";
    spin_image.resize(num_radial_bins, num_height_bins);
    spin_image.setZero();
    return;
  }

  compute(position, surf->getSample(nn_index).getNormal(), num_radial_bins, num_height_bins, spin_image);
}

void
SpinImage::compute(Vector3 const & position, Vector3 const & normal, int num_radial_bins, int num_height_bins,
                   MatrixX<double> & spin_image) const
{
  alwaysAssertM(num_radial_bins > 0, "SpinImage: Number of radial bins must be positive");
  alwaysAssertM(num_height_bins > 0, "SpinImage: Number of height bins must be positive");

  // Guess suitable limiting extents of the spin image, to have enough non-zero bins while not having too much of the
  // surface outside these extents
  Real max_radius = 0.75f * surf->getScale();  // radial limits: [0, max_radius]
  Real max_height = 0.75f * surf->getScale();  // height limits: [-max_height, max_height]

  spin_image.resize(num_radial_bins, num_height_bins);
  spin_image.setZero();

  Vector3 axis = normal.stableNormalized();
  intx n = surf->numSamples();
  for (intx i = 0; i < n; ++i)
  {
    Vector3 offset = surf->getSample(i).getPosition() - position;
    Real height = offset.dot(axis);
    Real radius = (offset - height * axis).norm();

    int radial_bin = (int)Math::clamp((int)std::floor(num_radial_bins * (radius / max_radius)), 0, num_radial_bins - 1);
    int height_bin = (int)Math::clamp((int)std::floor(num_height_bins * (0.5f * (height / max_height + 1))),
                                                      0, num_height_bins - 1);

    spin_image(radial_bin, height_bin) += 1;
  }
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
