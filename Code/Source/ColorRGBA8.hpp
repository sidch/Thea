//============================================================================
//
// This file is part of the Thea project.
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

 @file Color4uint8.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2003-04-07
 @edited  2010-03-24

 Copyright 2000-2010, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_ColorRGBA8_hpp__
#define __Thea_ColorRGBA8_hpp__

#include "Common.hpp"
#include "ColorL.hpp"
#include "ColorRGB8.hpp"
#include "Math.hpp"

namespace Thea {

// Forward declarations
class ColorRGBA;

/**
 * A color with four byte-sized channels: red, green, blue and alpha, each in [0, 255]. Derived from the G3D library:
 * http://g3d.sourceforge.net
 */
THEA_BEGIN_PACKED_CLASS(1)
class THEA_API ColorRGBA8
{
  private:
    uint8 c[4];  ///< Four components: red, green, blue and alpha.

  public:
    /** Default constructor. Does not initialize fields. */
    ColorRGBA8() {}

    /** Construct from red, green, blue and alpha components. */
    ColorRGBA8(uint8 r_, uint8 g_, uint8 b_, uint8 a_) { c[0] = r_; c[1] = g_; c[2] = b_; c[3] = a_; }

    /** Construct from an RGB color and an alpha component. */
    ColorRGBA8(ColorRGB8 const & rgb_, uint8 a_ = 255) { c[0] = rgb_.r(); c[1] = rgb_.g(); c[2] = rgb_.b(); c[3] = a_; }

    /** Construct a color from three components in an array. */
    explicit ColorRGBA8(uint8 v[4]) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Copy constructor. */
    ColorRGBA8(ColorRGBA8 const & other) { c[0] = other.c[0]; c[1] = other.c[1]; c[2] = other.c[2]; c[3] = other.c[3]; }

    /** Construct from a color with floating-point channels, with automatic scaling from [0, 1] to [0, 255]. */
    ColorRGBA8(ColorRGBA const & src);

    /** Initialize from an HTML-style color (e.g. 0xFF0000 == RED) */
    static ColorRGBA8 fromARGB(uint32 argb);

    /** The value of the red channel. */
    uint8 r() const { return c[0]; }

    /** A reference to the red channel. */
    uint8 & r() { return c[0]; }

    /** The value of the green channel. */
    uint8 g() const { return c[1]; }

    /** A reference to the green channel. */
    uint8 & g() { return c[1]; }

    /** The value of the blue channel. */
    uint8 b() const { return c[2]; }

    /** A reference to the blue channel. */
    uint8 & b() { return c[2]; }

    /** The value of the alpha channel. */
    uint8 a() const { return c[3]; }

    /** A reference to the alpha channel. */
    uint8 & a() { return c[3]; }

    /** Get the red, green and blue channels as a ColorRGB8. */
    ColorRGB8 rgb() const { return ColorRGB8(c[0], c[1], c[2]); }

    /** Set all channels simultaneously. */
    void set(uint8 r_, uint8 g_, uint8 b_, uint8 a_)
    {
      c[0] = r_;
      c[1] = g_;
      c[2] = b_;
      c[3] = a_;
    }

    /** Addition. Upper-bounds channels to 255. */
    ColorRGBA8 operator+(ColorRGBA8 const & rhs) const
    {
      return ColorRGBA8((uint8)std::min(255, (int)c[0] + (int)rhs.c[0]),
                        (uint8)std::min(255, (int)c[1] + (int)rhs.c[1]),
                        (uint8)std::min(255, (int)c[2] + (int)rhs.c[2]),
                        (uint8)std::min(255, (int)c[3] + (int)rhs.c[3]));
    }

    /** Subtraction. Lower-bounds channels to 0. */
    ColorRGBA8 operator-(ColorRGBA8 const & rhs) const
    {
      return ColorRGBA8((uint8)std::max(0, (int)c[0] - (int)rhs.c[0]),
                        (uint8)std::max(0, (int)c[1] - (int)rhs.c[1]),
                        (uint8)std::max(0, (int)c[2] - (int)rhs.c[2]),
                        (uint8)std::max(0, (int)c[3] - (int)rhs.c[3]));
    }

    /** Multiplication by a scalar. Channels are rounded to the nearest byte values. */
    ColorRGBA8 operator*(Real s) const
    {
      return ColorRGBA8((uint8)Math::clamp((Real)Math::round(c[0] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[1] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[2] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[3] * s), (Real)0, (Real)255));
    }

    /** Division by a scalar. Channels are rounded to the nearest byte values. */
    ColorRGBA8 operator/(Real s) const
    {
      return ColorRGBA8((uint8)Math::clamp((Real)Math::round(c[0] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[1] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[2] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[3] / s), (Real)0, (Real)255));
    }

    /** Add and assign. Upper-bounds channels to 255. */
    ColorRGBA8 & operator+=(ColorRGBA8 const & rhs)
    {
      *this = *this + rhs;
      return *this;
    }

    /** Subtract and assign. Lower-bounds channels to 0. */
    ColorRGBA8 & operator-=(ColorRGBA8 const & rhs)
    {
      *this = *this - rhs;
      return *this;
    }

    /** Multiply by a scalar and assign. Channels are rounded to the nearest byte values. */
    ColorRGBA8 & operator*=(Real s)
    {
      *this = *this * s;
      return *this;
    }

    /** Divide by a scalar and assign. Channels are rounded to the nearest byte values. */
    ColorRGBA8 & operator/=(Real s)
    {
      *this = *this / s;
      return *this;
    }

    /** Swap the red and blue channels. */
    ColorRGBA8 bgra() const
    {
      return ColorRGBA8(c[2], c[1], c[0], c[3]);
    }

    /**
     * Get the color packed into a 32-bit unsigned integer. The most significant byte is alpha, the next byte is red, then
     * green, then blue. Note that the actual memory ordering of the bytes depends upon the endianness of the system.
     */
    uint32 asUInt32() const
    {
      return ((uint32)c[3] << 24) | ((uint32)c[0] << 16) | ((uint32)c[1] << 8) | (uint32)c[2];
    }

    /** Get a string representation of the color. */
    std::string toString() const;

} THEA_END_PACKED_CLASS(1)  // class ColorRGBA8

/** Multiply by a scalar. */
inline ColorRGBA8
operator*(Real s, ColorRGBA8 const & c)
{
  return c * s;
}

/** Multiply by a one-channel color. */
inline ColorRGBA8
operator*(ColorL & s, ColorRGBA8 const & c)
{
  return c * s.value();
}

/** Multiply by a one-channel color. */
inline ColorRGBA8
operator*(ColorRGBA8 const & c, ColorL & s)
{
  return c * s.value();
}

} // namespace Thea

#endif
