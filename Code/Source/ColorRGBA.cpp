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

 @file ColorRGBA.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu
 @cite Portions by Laura Wollstadt, graphics3d.com
 @cite Portions based on Dave Eberly's Magic Software Library at http://www.magic-software.com


 @created 2002-06-25
 @edited  2009-11-10
*/

#include "ColorRGBA.hpp"
#include "ColorRGBA8.hpp"

namespace Thea {

ColorRGBA const &
ColorRGBA::zero()
{
  static ColorRGBA const col(0, 0, 0, 0);
  return col;
}

ColorRGBA::ColorRGBA(ColorRGBA8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
  c[3] = src.a() / 255.0f;
}

ColorRGBA
ColorRGBA::fromARGB(uint32 x)
{
  return ColorRGBA((float)((x >> 16) & 0xFF),
                   (float)((x >>  8) & 0xFF),
                   (float)( x        & 0xFF),
                   (float)((x >> 24) & 0xFF)) / 255.0;
}

std::string
ColorRGBA::toString() const
{
  return format("RGBA(%g, %g, %g, %g)", c[0], c[1], c[2], c[3]);
}

} // namespace Thea
