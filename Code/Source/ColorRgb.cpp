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

 @file ColorRgb.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu

 @created 2001-06-02
 @edited  2010-01-28
*/

#include "ColorRgb.hpp"
#include "ColorRgb8.hpp"
#include "ColorRgba.hpp"
#include "Crypto.hpp"

namespace Thea {

ColorRgb
ColorRgb::ansiMap(int i)
{
  static ColorRgb const map[] =
  {
    ColorRgb::black(),
    ColorRgb::red()    * 0.75f,
    ColorRgb::green()  * 0.75f,
    ColorRgb::yellow() * 0.75f,
    ColorRgb::blue()   * 0.75f,
    ColorRgb::purple() * 0.75f,
    ColorRgb::cyan()   * 0.75f,
    ColorRgb::white()  * 0.75f,
    ColorRgb::white()  * 0.90f,
    ColorRgb::red(),
    ColorRgb::green(),
    ColorRgb::yellow(),
    ColorRgb::blue(),
    ColorRgb::purple(),
    ColorRgb::cyan(),
    ColorRgb::white()
  };

  return map[i % 16];
}

ColorRgb
ColorRgb::pastelMap(int i)
{
  uint32 x = Crypto::crc32(&i, sizeof(uint32));

  // Create fairly bright, saturated colors
  Vector3 v(((x >> 22) & 1023) / 1023.0f,
            (((x >> 11) & 2047) / 2047.0f) * 0.5f + 0.25f,
            ((x & 2047) / 2047.0f) * 0.75f + 0.25f);
  return ColorRgb::fromHSV(v);
}

ColorRgb const &
ColorRgb::red()
{
  static ColorRgb const col(1.0f, 0.0f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::green()
{
  static ColorRgb const col(0.0f, 1.0f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::blue()
{
  static ColorRgb const col(0.0f, 0.0f, 1.0f);
  return col;
}

ColorRgb const &
ColorRgb::purple()
{
  static ColorRgb const col(0.7f, 0.0f, 1.0f);
  return col;
}

ColorRgb const &
ColorRgb::cyan()
{
  static ColorRgb const col(0.0f, 0.7f, 1.0f);
  return col;
}

ColorRgb const &
ColorRgb::yellow()
{
  static ColorRgb const col(1.0f, 1.0f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::brown()
{
  static ColorRgb const col(0.5f, 0.5f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::orange()
{
  static ColorRgb const col(1.0f, 0.5f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::black()
{
  static ColorRgb const col(0.0f, 0.0f, 0.0f);
  return col;
}

ColorRgb const &
ColorRgb::gray()
{
  static ColorRgb const col(0.7f, 0.7f, 0.7f);
  return col;
}

ColorRgb const &
ColorRgb::white()
{
  static ColorRgb const col(1, 1, 1);
  return col;
}

ColorRgb const &
ColorRgb::zero()
{
  static ColorRgb const col(0, 0, 0);
  return col;
}

ColorRgb const & ColorRgb::wheelRandom()
{
  static ColorRgb const color_array[8] =
  {
    ColorRgb::blue(),   ColorRgb::red(),    ColorRgb::green(),
    ColorRgb::orange(), ColorRgb::yellow(),
    ColorRgb::cyan(),   ColorRgb::purple(), ColorRgb::brown()
  };

  return color_array[Random::common().integer(0, 7)];
}

ColorRgb::ColorRgb(ColorRgba const & src)
{
  *this = src.rgb();
}

ColorRgb::ColorRgb(ColorRgb8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
}

ColorRgb
ColorRgb::fromARGB(uint32 x)
{
  return ColorRgb((Real)((x >> 16) & 0xFF), (Real)((x >> 8) & 0xFF), (Real)(x & 0xFF)) / 255.0f;
}

ColorRgb ColorRgb::random()
{
  return ColorRgb(Random::common().uniform01(), Random::common().uniform01(), Random::common().uniform01()).normalized();
}

ColorRgb
ColorRgb::fromHSV(Vector3 const & hsv)
{
  debugAssertM((hsv[0] <= 1.0f && hsv[0] >= 0.0f)
            && (hsv[1] <= 1.0f && hsv[1] >= 0.0f)
            && (hsv[2] <= 1.0f && hsv[2] >= 0.0f), "ColorRgb: H, S, V must be in [0, 1]");

  int const i = std::min(5, (int)std::floor(6.0 * hsv[0]));
  Real const f = 6.0f * hsv[0] - i;
  Real const m = hsv[2] * (1.0f - (hsv[1]));
  Real const n = hsv[2] * (1.0f - (hsv[1] * f));
  Real const k = hsv[2] * (1.0f - (hsv[1] * (1 - f)));

  switch (i)
  {
    case 0: return ColorRgb(hsv[2], k, m);
    case 1: return ColorRgb(n, hsv[2], m);
    case 2: return ColorRgb(m, hsv[2], k);
    case 3: return ColorRgb(m, n, hsv[2]);
    case 4: return ColorRgb(k, m, hsv[2]);
    case 5: return ColorRgb(hsv[2], m, n);

    default: debugAssertM(false, "ColorRgb: Fell through switch when attempting conversion from HSV");
  }

  return ColorRgb::black();
}

Vector3
ColorRgb::toHSV() const
{
  debugAssertM((c[0] <= 1.0f && c[0] >= 0.0f)
            && (c[1] <= 1.0f && c[1] >= 0.0f)
            && (c[2] <= 1.0f && c[2] >= 0.0f), "ColorRgb: R, G, B must be in [0, 1]");

  Vector3 hsv = Vector3::Zero();
  hsv[2] = std::max(std::max(c[0], c[1]), c[2]);

  if (Math::fuzzyEq(hsv[2], (Real)0))
    return hsv;

  Real const x =  std::min(std::min(c[0], c[1]), c[2]);
  hsv[1] = (hsv[2] - x) / hsv[2];

  if (Math::fuzzyEq(hsv[1], (Real)0))
    return hsv;

  Vector3 rgbN;
  rgbN[0] = (hsv[2] - c[0]) / (hsv[2] - x);
  rgbN[1] = (hsv[2] - c[1]) / (hsv[2] - x);
  rgbN[2] = (hsv[2] - c[2]) / (hsv[2] - x);

  if (c[0] == hsv[2])  // Note: from the max we know that it exactly equals one of the three
  {
    hsv[0] = (c[1] == x) ? 5.0f + rgbN[2] : 1.0f - rgbN[1];
  }
  else if (c[1] == hsv[2])
  {
    hsv[0] = (c[2] == x) ? 1.0f + rgbN[0] : 3.0f - rgbN[2];
  }
  else
  {
    hsv[0] = (c[0] == x) ? 3.0f + rgbN[1] : 5.0f - rgbN[0];
  }

  hsv[0] /= 6.0f;

  return hsv;
}

namespace ColorRgbInternal {

// http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
// (solution there is for val in [-1, 1])
Real
interpolate(Real val, Real y0, Real x0, Real y1, Real x1)
{
  return (val - x0) * (y1 - y0) / (x1 - x0) + y0;
}

Real
base(Real val)
{
  if      (val <= 0.125) return 0;
  else if (val <= 0.375) return interpolate(val, 0.0, 0.125, 1.0, 0.375);
  else if (val <= 0.625) return 1.0;
  else if (val <= 0.875) return interpolate(val, 1.0, 0.625, 0.0, 0.875);
  else                   return 0.0;
}

} // namespace ColorRgbInternal

ColorRgb
ColorRgb::jetColorMap(Real val)
{
  val = Math::clamp(val, 0, 1);
  return ColorRgb(ColorRgbInternal::base(val - 0.25),
                  ColorRgbInternal::base(val),
                  ColorRgbInternal::base(val + 0.25));
}

std::string
ColorRgb::toString() const
{
  return format("RGB(%g, %g, %g)", c[0], c[1], c[2]);
}

ColorRgb
ColorRgb::rainbowColorMap(Real hue)
{
  return fromHSV(Vector3(hue, 1, 1));
}

} // namespace Thea
