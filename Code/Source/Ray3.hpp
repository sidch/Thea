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

#ifndef __Thea_Ray3_hpp__
#define __Thea_Ray3_hpp__

#include "Common.hpp"
#include "RayN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API RayN<3, Real>;
#endif

// The default ray class in 3-dimensional real space.
typedef RayN<3, Real> Ray3;

} // namespace Thea

#endif
