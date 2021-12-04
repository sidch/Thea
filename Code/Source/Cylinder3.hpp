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
// First version: 2021
//
//============================================================================

#ifndef __Thea_Cylinder3_hpp__
#define __Thea_Cylinder3_hpp__

#include "Common.hpp"
#include "CylinderN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API CylinderN<3, Real>;
#endif

/** The default cylinder class in real 3-space. */
typedef CylinderN<3, Real> Cylinder3;

} // namespace Thea

#endif
