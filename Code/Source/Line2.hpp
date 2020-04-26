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
// First version: 2011
//
//============================================================================

#ifndef __Thea_Line2_hpp__
#define __Thea_Line2_hpp__

#include "Common.hpp"
#include "LineN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API LineN<2, Real>;
#endif

/** The default straight line class in 2-dimensional real space. */
typedef LineN<2, Real> Line2;

} // namespace Thea

#endif
