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

 @file ColorRgba.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu
 @cite Portions by Laura Wollstadt, graphics3d.com
 @cite Portions based on Dave Eberly's Magic Software Library at http://www.magic-software.com


 @created 2002-06-25
 @edited  2009-11-10
*/

#include "ColorRgba.hpp"
#include "ColorRgb8.hpp"
#include "ColorRgba8.hpp"

namespace Thea {

ColorRgba const & ColorRgba::red()    { static ColorRgba const col(ColorRgb::red(),    1); return col; }
ColorRgba const & ColorRgba::green()  { static ColorRgba const col(ColorRgb::green(),  1); return col; }
ColorRgba const & ColorRgba::blue()   { static ColorRgba const col(ColorRgb::blue(),   1); return col; }
ColorRgba const & ColorRgba::purple() { static ColorRgba const col(ColorRgb::purple(), 1); return col; }
ColorRgba const & ColorRgba::cyan()   { static ColorRgba const col(ColorRgb::cyan(),   1); return col; }
ColorRgba const & ColorRgba::yellow() { static ColorRgba const col(ColorRgb::yellow(), 1); return col; }
ColorRgba const & ColorRgba::brown()  { static ColorRgba const col(ColorRgb::brown(),  1); return col; }
ColorRgba const & ColorRgba::orange() { static ColorRgba const col(ColorRgb::orange(), 1); return col; }
ColorRgba const & ColorRgba::black()  { static ColorRgba const col(ColorRgb::black(),  1); return col; }
ColorRgba const & ColorRgba::gray()   { static ColorRgba const col(ColorRgb::gray(),   1); return col; }
ColorRgba const & ColorRgba::white()  { static ColorRgba const col(ColorRgb::white(),  1); return col; }

ColorRgba const &
ColorRgba::zero()
{
  static ColorRgba const col(0, 0, 0, 0);
  return col;
}

ColorRgba::ColorRgba(ColorRgba8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
  c[3] = src.a() / 255.0f;
}

ColorRgba::ColorRgba(ColorRgb8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
  c[3] = 1;
}

ColorRgba
ColorRgba::fromARGB(uint32 x)
{
  return ColorRgba((float)((x >> 16) & 0xFF),
                   (float)((x >>  8) & 0xFF),
                   (float)( x        & 0xFF),
                   (float)((x >> 24) & 0xFF)) / 255.0;
}

std::string
ColorRgba::toString() const
{
  return format("RGBA(%g, %g, %g, %g)", c[0], c[1], c[2], c[3]);
}

} // namespace Thea
