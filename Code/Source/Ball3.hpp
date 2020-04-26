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

#ifndef __Thea_Ball3_hpp__
#define __Thea_Ball3_hpp__

#include "Common.hpp"
#include "BallN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API BallN<3, Real>;
#endif

/** The default ball class in real 3-space. */
typedef BallN<3, Real> Ball3;

} // namespace Thea

#endif
