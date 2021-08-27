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

#include "Curvature.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../NormalTraitsN.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Noncopyable.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

Curvature::Curvature(PointSet3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "Curvature: Cannot construct with a null surface");
  alwaysAssertM(surf_->hasNormals(), "Curvature: Cannot construct with a surface whose samples lack normals");
}

namespace CurvatureInternal {

// Called for each point in the neighborhood.
struct ProjectedCurvatureFunctor : private Noncopyable
{
  ProjectedCurvatureFunctor(Vector3 const & p, Vector3 const & n)
  : position(p), normal(n), num_offsets(0), sum_offsets(Vector3::Zero())
  {}

  template <typename SampleT> bool operator()(intx index, SampleT const & t)
  {
    if (NormalTraitsN<SampleT, 3>::getNormal(t).dot(normal) > -1.0e-05f)  // ignore points on hidden side
    {
      Vector3 offset = (PointTraitsN<SampleT, 3>::getPosition(t) - position).stableNormalized();
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

} // namespace CurvatureInternal

double
Curvature::computeProjectedCurvature(Vector3 const & position, Real nbd_radius) const
{
  intx nn_index = surf->getKdTree().closestElement<MetricL2>(position);
  if (nn_index < 0)
  {
    THEA_WARNING << "Curvature: Query point cannot be mapped to mesh, curvature value set to zero";
    return 0.0;
  }

  return computeProjectedCurvature(position, surf->getSample(nn_index).getNormal(), nbd_radius);
}

double
Curvature::computeProjectedCurvature(Vector3 const & position, Vector3 const & normal, Real nbd_radius) const
{
  if (nbd_radius < 0)
    nbd_radius = 0.1f;

  nbd_radius *= surf->getScale();

  CurvatureInternal::ProjectedCurvatureFunctor func(position, normal);
  Ball3 range(position, nbd_radius);
  const_cast<PointSet3::SampleKdTree &>(surf->getKdTree()).processRangeUntil<IntersectionTester>(range, std::ref(func));

  return func.getCurvature();
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
