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
class ColorRGBA;

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
    ColorL8(ColorRGBA const & other);

    /** Copy constructor. */
    ColorL8(ColorL8 const & other) : val(other.val) {}

    /** Initialize from a floating-point color, automatically multiplying by 255. */
    ColorL8(ColorL const & other);

    /** The value of the color. */
    uint8 value() const { return val; }

    /** The value of the color. */
    uint8 & value() { return val; }

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
