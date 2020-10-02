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

  @file Color1uint8uint8.h

  @maintainer Morgan McGuire, graphics3d.com

  @created 2007-01-30
  @edited  2007-01-30

  Copyright 2000-2007, Morgan McGuire.
  All rights reserved.
 */

#ifndef __Thea_ColorL8_hpp__
#define __Thea_ColorL8_hpp__

#include "Common.hpp"
#include "Math.hpp"

namespace Thea {

// Forward declarations
class ColorL;
class ColorRgba;

/**
 * Monochrome color represented as a single byte value in [0, 255], with automatic scaling by 255 when switching between integer
 * (ColorL8) and floating point (ColorL) formats. Derived from the G3D library: http://g3d.sourceforge.net
 */
THEA_BEGIN_PACKED_CLASS(1)
class THEA_API ColorL8
{
  private:
    uint8 val;  ///< Luminance value.

  public:
    /** Default constructor, initializes color to 0. */
    ColorL8() : val(0) {}

    /** Initializing constructor. */
    explicit ColorL8(uint8 v) : val(v) {}

    /** Initialize from a 32-bit RGBA color. For conversion from a consistent source type. */
    ColorL8(ColorRgba const & other);

    /** Initialize from a floating-point color, automatically multiplying by 255. */
    ColorL8(ColorL const & other);

    /** The value of the color. */
    uint8 value() const { return val; }

    /** The value of the color. */
    uint8 & value() { return val; }

    /** Get the address of the color value. */
    uint8 const * data() const { return &val; }

    /** Get the address of the color value. */
    uint8 * data() { return &val; }

    /** Array-style channel access. */
    template <typename IntegerT> uint8 const & operator[](IntegerT channel) const
    {
      debugAssertM(channel == 0, "ColorL8: Channel must be 0");
      return val;
    }

    /** Array-style channel access. */
    template <typename IntegerT> uint8 & operator[](IntegerT channel)
    {
      debugAssertM(channel == 0, "ColorL8: Channel must be 0");
      return val;
    }

    /** Addition operator. Upper-bounds result to 255. */
    ColorL8 operator+(ColorL8 const & other) const
    {
      return ColorL8((uint8)std::min(255, (int)val + (int)other.val));
    }

    /** Subtraction operator. Lower-bounds result to 0. */
    ColorL8 operator-(ColorL8 const & other) const
    {
      return ColorL8((uint8)std::max(0, (int)val - (int)other.val));
    }

    /** Multiply by a scalar. The result is rounded to the nearest byte value. */
    ColorL8 operator*(Real other) const
    {
      return ColorL8((uint8)Math::clamp((Real)Math::round(val * other), (Real)0, (Real)255));
    }

    /** Divide by a scalar. The result is rounded to the nearest byte value. */
    ColorL8 operator/(Real other) const
    {
      return ColorL8((uint8)Math::clamp((Real)Math::round(val / other), (Real)0, (Real)255));
    }

    /** Add-and-assign. Upper-bounds result to 255. */
    ColorL8 & operator+=(ColorL8 const & other)
    {
      *this = *this + other;
      return *this;
    }

    /** Subtract-and-assign. Lower-bounds result to 0. */
    ColorL8 & operator-=(ColorL8 const & other)
    {
      *this = *this - other;
      return *this;
    }

    /** Multiply by a scalar and assign. The result is rounded to the nearest byte value. */
    ColorL8 & operator*=(Real other)
    {
      *this = *this * other;
      return *this;
    }

    /** Divide by a scalar and assign. The result is rounded to the nearest byte value. */
    ColorL8 & operator/=(Real other)
    {
      *this = *this / other;
      return *this;
    }

    /** Get the maximum of two color values. */
    ColorL8 max(ColorL8 const & other) const
    {
      return ColorL8(std::max(val, other.val));
    }

    /** Get the minimum of two color values. */
    ColorL8 min(ColorL8 const & other) const
    {
      return ColorL8(std::min(val, other.val));
    }

    /** Color with zero luminance (black). */
    static ColorL8 const & zero() { static ColorL8 const col(0); return col; }

    /** Get a string representation of the color. */
    std::string toString() const;

} THEA_END_PACKED_CLASS(1)  // class ColorL8

/** Multiply by a scalar. The result is rounded to the nearest byte value. */
inline ColorL8
operator*(Real lhs, ColorL8 const & rhs)
{
  return rhs * lhs;
}

} // namespace Thea

#endif
