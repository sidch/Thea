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

 @file Color4uint8.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2003-04-07
 @edited  2010-03-24

 Copyright 2000-2010, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_ColorRgba8_hpp__
#define __Thea_ColorRgba8_hpp__

#include "Common.hpp"
#include "ColorL.hpp"
#include "ColorRgb8.hpp"
#include "Math.hpp"

namespace Thea {

// Forward declarations
class ColorRgba;

/**
 * A color with four byte-sized channels: red, green, blue and alpha, each in [0, 255]. Derived from the G3D library:
 * http://g3d.sourceforge.net
 */
THEA_BEGIN_PACKED_CLASS(1)
class THEA_API ColorRgba8
{
  private:
    uint8 c[4];  ///< Four components: red, green, blue and alpha.

  public:
    /** Default constructor. Does not initialize fields. */
    ColorRgba8() {}

    /** Construct from red, green, blue and alpha components. */
    ColorRgba8(uint8 r_, uint8 g_, uint8 b_, uint8 a_) { c[0] = r_; c[1] = g_; c[2] = b_; c[3] = a_; }

    /** Construct from an RGB color and an alpha component. */
    ColorRgba8(ColorRgb8 const & rgb_, uint8 a_ = 255) { c[0] = rgb_.r(); c[1] = rgb_.g(); c[2] = rgb_.b(); c[3] = a_; }

    /** Construct a color from four components in an array. */
    explicit ColorRgba8(uint8 const * v) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Construct from a color with floating-point channels, with automatic scaling from [0, 1] to [0, 255]. */
    ColorRgba8(ColorRgba const & src);

    /** Initialize from an HTML-style color (e.g. 0xFF0000 == RED) */
    static ColorRgba8 fromARGB(uint32 argb);

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

    /** Get the red, green and blue channels as a ColorRgb8. */
    ColorRgb8 rgb() const { return ColorRgb8(c[0], c[1], c[2]); }

    /** Get the address of the array storing color channel values in RGBA order. */
    uint8 const * data() const { return c; }

    /** Get the address of the array storing color channel values in RGBA order. */
    uint8 * data() { return c; }

    /** Array-style channel access. */
    template <typename IntegerT> uint8 const & operator[](IntegerT channel) const
    {
      debugAssertM(channel >= 0 && channel <= 3, "ColorRgba8: Channel must be 0, 1, 2 or 3");
      return c[channel];
    }

    /** Array-style channel access. */
    template <typename IntegerT> uint8 & operator[](IntegerT channel)
    {
      debugAssertM(channel >= 0 && channel <= 3, "ColorRgba8: Channel must be 0, 1, 2 or 3");
      return c[channel];
    }

    /** Set all channels simultaneously. */
    void set(uint8 r_, uint8 g_, uint8 b_, uint8 a_)
    {
      c[0] = r_;
      c[1] = g_;
      c[2] = b_;
      c[3] = a_;
    }

    /** Addition. Upper-bounds channels to 255. */
    ColorRgba8 operator+(ColorRgba8 const & rhs) const
    {
      return ColorRgba8((uint8)std::min(255, (int)c[0] + (int)rhs.c[0]),
                        (uint8)std::min(255, (int)c[1] + (int)rhs.c[1]),
                        (uint8)std::min(255, (int)c[2] + (int)rhs.c[2]),
                        (uint8)std::min(255, (int)c[3] + (int)rhs.c[3]));
    }

    /** Subtraction. Lower-bounds channels to 0. */
    ColorRgba8 operator-(ColorRgba8 const & rhs) const
    {
      return ColorRgba8((uint8)std::max(0, (int)c[0] - (int)rhs.c[0]),
                        (uint8)std::max(0, (int)c[1] - (int)rhs.c[1]),
                        (uint8)std::max(0, (int)c[2] - (int)rhs.c[2]),
                        (uint8)std::max(0, (int)c[3] - (int)rhs.c[3]));
    }

    /** Multiplication by a scalar. Channels are rounded to the nearest byte values. */
    ColorRgba8 operator*(Real s) const
    {
      return ColorRgba8((uint8)Math::clamp((Real)Math::round(c[0] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[1] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[2] * s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[3] * s), (Real)0, (Real)255));
    }

    /** Division by a scalar. Channels are rounded to the nearest byte values. */
    ColorRgba8 operator/(Real s) const
    {
      return ColorRgba8((uint8)Math::clamp((Real)Math::round(c[0] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[1] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[2] / s), (Real)0, (Real)255),
                        (uint8)Math::clamp((Real)Math::round(c[3] / s), (Real)0, (Real)255));
    }

    /** Add and assign. Upper-bounds channels to 255. */
    ColorRgba8 & operator+=(ColorRgba8 const & rhs)
    {
      *this = *this + rhs;
      return *this;
    }

    /** Subtract and assign. Lower-bounds channels to 0. */
    ColorRgba8 & operator-=(ColorRgba8 const & rhs)
    {
      *this = *this - rhs;
      return *this;
    }

    /** Multiply by a scalar and assign. Channels are rounded to the nearest byte values. */
    ColorRgba8 & operator*=(Real s)
    {
      *this = *this * s;
      return *this;
    }

    /** Divide by a scalar and assign. Channels are rounded to the nearest byte values. */
    ColorRgba8 & operator/=(Real s)
    {
      *this = *this / s;
      return *this;
    }

    /** Swap the red and blue channels. */
    ColorRgba8 bgra() const
    {
      return ColorRgba8(c[2], c[1], c[0], c[3]);
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

} THEA_END_PACKED_CLASS(1)  // class ColorRgba8

/** Multiply by a scalar. */
inline ColorRgba8
operator*(Real s, ColorRgba8 const & c)
{
  return c * s;
}

/** Multiply by a one-channel color. */
inline ColorRgba8
operator*(ColorL & s, ColorRgba8 const & c)
{
  return c * s.value();
}

/** Multiply by a one-channel color. */
inline ColorRgba8
operator*(ColorRgba8 const & c, ColorL & s)
{
  return c * s.value();
}

} // namespace Thea

#endif
