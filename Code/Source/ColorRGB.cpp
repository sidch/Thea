//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

/*
 ORIGINAL HEADER

 @file ColorRGB.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu

 @created 2001-06-02
 @edited  2010-01-28
*/

#include "ColorRGB.hpp"
#include "ColorRGB8.hpp"
#include "ColorRGBA.hpp"
#include "Crypto.hpp"

namespace Thea {

ColorRGB
ColorRGB::ansiMap(int i)
{
  static ColorRGB const map[] =
  {
    ColorRGB::black(),
    ColorRGB::red()    * 0.75f,
    ColorRGB::green()  * 0.75f,
    ColorRGB::yellow() * 0.75f,
    ColorRGB::blue()   * 0.75f,
    ColorRGB::purple() * 0.75f,
    ColorRGB::cyan()   * 0.75f,
    ColorRGB::white()  * 0.75f,
    ColorRGB::white()  * 0.90f,
    ColorRGB::red(),
    ColorRGB::green(),
    ColorRGB::yellow(),
    ColorRGB::blue(),
    ColorRGB::purple(),
    ColorRGB::cyan(),
    ColorRGB::white()
  };

  return map[i % 16];
}

ColorRGB
ColorRGB::pastelMap(int i)
{
  uint32 x = Crypto::crc32(&i, sizeof(uint32));

  // Create fairly bright, saturated colors
  Vector3 v(((x >> 22) & 1023) / 1023.0f,
            (((x >> 11) & 2047) / 2047.0f) * 0.5f + 0.25f,
            ((x & 2047) / 2047.0f) * 0.75f + 0.25f);
  return ColorRGB::fromHSV(v);
}

ColorRGB const &
ColorRGB::red()
{
  static ColorRGB const col(1.0f, 0.0f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::green()
{
  static ColorRGB const col(0.0f, 1.0f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::blue()
{
  static ColorRGB const col(0.0f, 0.0f, 1.0f);
  return col;
}

ColorRGB const &
ColorRGB::purple()
{
  static ColorRGB const col(0.7f, 0.0f, 1.0f);
  return col;
}

ColorRGB const &
ColorRGB::cyan()
{
  static ColorRGB const col(0.0f, 0.7f, 1.0f);
  return col;
}

ColorRGB const &
ColorRGB::yellow()
{
  static ColorRGB const col(1.0f, 1.0f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::brown()
{
  static ColorRGB const col(0.5f, 0.5f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::orange()
{
  static ColorRGB const col(1.0f, 0.5f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::black()
{
  static ColorRGB const col(0.0f, 0.0f, 0.0f);
  return col;
}

ColorRGB const &
ColorRGB::gray()
{
  static ColorRGB const col(0.7f, 0.7f, 0.7f);
  return col;
}

ColorRGB const &
ColorRGB::white()
{
  static ColorRGB const col(1, 1, 1);
  return col;
}

ColorRGB const &
ColorRGB::zero()
{
  static ColorRGB const col(0, 0, 0);
  return col;
}

ColorRGB const & ColorRGB::wheelRandom()
{
  static ColorRGB const color_array[8] =
  {
    ColorRGB::blue(),   ColorRGB::red(),    ColorRGB::green(),
    ColorRGB::orange(), ColorRGB::yellow(),
    ColorRGB::cyan(),   ColorRGB::purple(), ColorRGB::brown()
  };

  return color_array[Random::common().integer(0, 7)];
}

ColorRGB::ColorRGB(ColorRGBA const & src)
{
  *this = src.rgb();
}

ColorRGB::ColorRGB(ColorRGB8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
}

ColorRGB
ColorRGB::fromARGB(uint32 x)
{
  return ColorRGB((Real)((x >> 16) & 0xFF), (Real)((x >> 8) & 0xFF), (Real)(x & 0xFF)) / 255.0f;
}

ColorRGB ColorRGB::random()
{
  return ColorRGB(Random::common().uniform01(), Random::common().uniform01(), Random::common().uniform01()).unit();
}

ColorRGB
ColorRGB::fromHSV(Vector3 const & hsv)
{
  debugAssertM((hsv[0] <= 1.0f && hsv[0] >= 0.0f)
            && (hsv[1] <= 1.0f && hsv[1] >= 0.0f)
            && (hsv[2] <= 1.0f && hsv[2] >= 0.0f), "ColorRGB: H, S, V must be in [0, 1]");

  int const i = std::min(5, (int)std::floor(6.0 * hsv[0]));
  Real const f = 6.0f * hsv[0] - i;
  Real const m = hsv[2] * (1.0f - (hsv[1]));
  Real const n = hsv[2] * (1.0f - (hsv[1] * f));
  Real const k = hsv[2] * (1.0f - (hsv[1] * (1 - f)));

  switch (i)
  {
    case 0: return ColorRGB(hsv[2], k, m);
    case 1: return ColorRGB(n, hsv[2], m);
    case 2: return ColorRGB(m, hsv[2], k);
    case 3: return ColorRGB(m, n, hsv[2]);
    case 4: return ColorRGB(k, m, hsv[2]);
    case 5: return ColorRGB(hsv[2], m, n);

    default: debugAssertM(false, "ColorRGB: Fell through switch when attempting conversion from HSV");
  }

  return ColorRGB::black();
}

Vector3
ColorRGB::toHSV() const
{
  debugAssertM((c[0] <= 1.0f && c[0] >= 0.0f)
            && (c[1] <= 1.0f && c[1] >= 0.0f)
            && (c[2] <= 1.0f && c[2] >= 0.0f), "ColorRGB: R, G, B must be in [0, 1]");

  Vector3 hsv = Vector3::zero();
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

namespace ColorRGBInternal {

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

} // namespace ColorRGBInternal

ColorRGB
ColorRGB::jetColorMap(Real val)
{
  val = Math::clamp(val, 0, 1);
  return ColorRGB(ColorRGBInternal::base(val - 0.25),
                  ColorRGBInternal::base(val),
                  ColorRGBInternal::base(val + 0.25));
}

std::string
ColorRGB::toString() const
{
  return format("RGB(%g, %g, %g)", c[0], c[1], c[2]);
}

ColorRGB
ColorRGB::rainbowColorMap(Real hue)
{
  return fromHSV(Vector3(hue, 1, 1));
}

} // namespace Thea
