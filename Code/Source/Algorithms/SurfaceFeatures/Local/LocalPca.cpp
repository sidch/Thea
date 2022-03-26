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
// First version: 2014
//
//============================================================================

#include "LocalPca.hpp"
#include "../../IntersectionTester.hpp"
#include "../../PcaN.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

LocalPca::LocalPca(PointSet3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "LocalPca: Cannot construct with a null surface");
}

Vector3
LocalPca::LocalPcaFunctor::getPcaFeatures(Vector3 * eigenvectors) const
{
  Real eval[3];
  Vector3 evec[3];
  PcaN<Vector3, 3>::compute(nbd_pts.begin(), nbd_pts.end(), eval, evec);  // returns ordered by decreasing eigenvalue

  if (eigenvectors)
  {
    eigenvectors[0] = evec[0];
    eigenvectors[1] = evec[1];
    eigenvectors[2] = evec[2];
  }

  return Vector3(eval[0], eval[1], eval[2]);
}

Vector3
LocalPca::compute(Vector3 const & position, Vector3 * eigenvectors, Real nbd_radius) const
{
  if (nbd_radius < 0)
    nbd_radius = 0.1f;

  nbd_radius *= surf->getScale();

  Ball3 range(position, nbd_radius);
  func.reset();
  surf->getBvh().processRange<IntersectionTester>(range, std::ref(func));

  return func.getPcaFeatures(eigenvectors);
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
