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

#ifndef __Thea_Cone3_hpp__
#define __Thea_Cone3_hpp__

#include "Common.hpp"
#include "ConeN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API ConeN<3, Real>;
#endif

/** The default cone class in real 3-space. */
typedef ConeN<3, Real> Cone3;

} // namespace Thea

#endif
