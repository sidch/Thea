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
// First version: 2009
//
//============================================================================

#ifndef __Thea_RayIntersectable3_hpp__
#define __Thea_RayIntersectable3_hpp__

#include "Common.hpp"
#include "RayIntersectableN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API RayIntersectionN<3, Real>;
  template class THEA_API RayIntersectableN<3, Real>;
#endif

// The default ray intersection class in real 3-space.
typedef RayIntersectionN<3, Real> RayIntersection3;

// The default ray intersectable interface in real 3-space.
typedef RayIntersectableN<3, Real> RayIntersectable3;

} // namespace Thea

#endif
