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

 @file ColorL.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu

 @created 2007-01-30
 @edited  2009-03-27
 */

#include "ColorL.hpp"
#include "ColorL8.hpp"
#include "ColorRgba.hpp"

namespace Thea {

ColorL::ColorL(ColorL8 const & other)
: val(other.value() / 255.0f)
{}

ColorL::ColorL(ColorRgba const & other)
: val(0.299f * other.r() + 0.587f * other.g() + 0.114f * other.b())
{}

std::string
ColorL::toString() const
{
  return format("L(%g)", val);
}

} // namespace Thea
