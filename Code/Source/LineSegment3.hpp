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

#ifndef __Thea_LineSegment3_hpp__
#define __Thea_LineSegment3_hpp__

#include "Common.hpp"
#include "LineSegmentN.hpp"

namespace Thea {

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API LineSegmentN<3, Real>;
#endif

/** The default line segment class in 3-dimensional real space. */
typedef LineSegmentN<3, Real> LineSegment3;

} // namespace Thea

#endif
