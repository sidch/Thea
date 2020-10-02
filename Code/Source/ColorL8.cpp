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

 @file ColorL8.cpp

 @author Morgan McGuire, http://graphics.cs.williams.edu

 @created 2007-01-30
 @edited  2007-01-30
 */

#include "ColorL8.hpp"
#include "ColorL.hpp"
#include "ColorRgba.hpp"
#include "Math.hpp"

namespace Thea {

ColorL8::ColorL8(ColorL const & c)
: val(Math::clamp((uint8)Math::round(c.value() * 255), (uint8)0, (uint8)255))
{}

ColorL::ColorL(ColorRgba const & other)
{
  Real lum = 0.299f * other.r() + 0.587f * other.g() + 0.114f * other.b();
  val = Math::clamp((uint8)Math::round(lum * 255), (uint8)0, (uint8)255);
}

std::string
ColorL8::toString() const
{
  return format("L8(%d)", (int)val);
}

} // namespace Thea
