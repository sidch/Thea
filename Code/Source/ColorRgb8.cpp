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

/*
 ORIGINAL HEADER

 @file Color3uint8.cpp

 @author Morgan McGuire, http://graphics.cs.williams.edu

 @created 2003-04-07
 @edited  2006-01-07
*/

#include "ColorRgb8.hpp"
#include "ColorRgb.hpp"
#include "ColorRgba.hpp"
#include "Math.hpp"

namespace Thea {

ColorRgb8::ColorRgb8(ColorRgb const & src)
{
  c[0] = (uint8)Math::clamp((Real)Math::round(src.r() * 255), (Real)0, (Real)255);
  c[1] = (uint8)Math::clamp((Real)Math::round(src.g() * 255), (Real)0, (Real)255);
  c[2] = (uint8)Math::clamp((Real)Math::round(src.b() * 255), (Real)0, (Real)255);
}

ColorRgb8::ColorRgb8(ColorRgba const & src)
{
  *this = ColorRgb8(src.rgb());
}

ColorRgb8
ColorRgb8::fromARGB(uint32 argb)
{
  return ColorRgb8((uint8)((argb >> 16) & 0xFF),
                   (uint8)((argb >>  8) & 0xFF),
                   (uint8)( argb        & 0xFF));
}

std::string
ColorRgb8::toString() const
{
  return format("RGB(%d, %d, %d)", (int)c[0], (int)c[1], (int)c[2]);
}

} // namespace Thea
