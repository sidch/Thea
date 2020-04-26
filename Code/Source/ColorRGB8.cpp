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

#include "ColorRGB8.hpp"
#include "ColorRGB.hpp"
#include "ColorRGBA.hpp"
#include "Math.hpp"

namespace Thea {

ColorRGB8::ColorRGB8(ColorRGB const & src)
{
  c[0] = (uint8)Math::clamp((Real)Math::round(src.r() * 255), (Real)0, (Real)255);
  c[1] = (uint8)Math::clamp((Real)Math::round(src.g() * 255), (Real)0, (Real)255);
  c[2] = (uint8)Math::clamp((Real)Math::round(src.b() * 255), (Real)0, (Real)255);
}

ColorRGB8::ColorRGB8(ColorRGBA const & src)
{
  *this = ColorRGB8(src.rgb());
}

ColorRGB8
ColorRGB8::fromARGB(uint32 argb)
{
  return ColorRGB8((uint8)((argb >> 16) & 0xFF),
                   (uint8)((argb >>  8) & 0xFF),
                   (uint8)( argb        & 0xFF));
}

std::string
ColorRGB8::toString() const
{
  return format("RGB(%d, %d, %d)", (int)c[0], (int)c[1], (int)c[2]);
}

} // namespace Thea
