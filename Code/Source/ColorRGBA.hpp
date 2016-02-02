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

 @file Color4.h

 Color class

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu
 @cite Portions based on Dave Eberly's Magic Software Library
      at <A HREF="http://www.magic-software.com">http://www.magic-software.com</A>

 @created 2002-06-25
 @edited  2009-11-15

 Copyright 2000-2009, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_ColorRGBA_hpp__
#define __Thea_ColorRGBA_hpp__

#include "Common.hpp"
#include "ColorL.hpp"
#include "ColorRGB.hpp"
#include "Math.hpp"
#include "Vector4.hpp"

namespace Thea {

// Forward declarations
class ColorRGBA8;

/**
 * A color with four floating-point channels: red, green and blue, each in [0, 1]. Derived from the G3D library:
 * http://g3d.sourceforge.net
 */
class THEA_API ColorRGBA
{
  private:
    Real c[4];  ///< Four components: red, gree, blue and alpha.

  public:
    /** Default constructor. Does not initialize fields. */
    ColorRGBA() {}

    /** Construct from red, green, blue and alpha components. */
    ColorRGBA(Real r_, Real g_, Real b_, Real a_) { c[0] = r_; c[1] = g_; c[2] = b_; c[3] = a_; }

    /** Construct from an RGB color and an alpha component. */
    ColorRGBA(ColorRGB const & rgb_, Real a_ = 1) { c[0] = rgb_.r(); c[1] = rgb_.g(); c[2] = rgb_.b(); c[3] = a_; }

    /** Construct a color from a 4-vector. */
    explicit ColorRGBA(Vector4 const & v) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Construct a color from three components in an array. */
    explicit ColorRGBA(Real v[4]) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Copy constructor. */
    ColorRGBA(ColorRGBA const & other) { c[0] = other.c[0]; c[1] = other.c[1]; c[2] = other.c[2]; c[3] = other.c[3]; }

    /** Construct from a color with byte channels, with automatic scaling from [0, 255] to [0, 1]. */
    ColorRGBA(ColorRGBA8 const & src);

    /** Initialize from an HTML-style color (e.g. 0xFFFF0000 == RED) */
    static ColorRGBA fromARGB(uint32 argb);

    /** The value of the red channel. */
    Real r() const { return c[0]; }

    /** A reference to the red channel. */
    Real & r() { return c[0]; }

    /** The value of the green channel. */
    Real g() const { return c[1]; }

    /** A reference to the green channel. */
    Real & g() { return c[1]; }

    /** The value of the blue channel. */
    Real b() const { return c[2]; }

    /** A reference to the blue channel. */
    Real & b() { return c[2]; }

    /** The value of the alpha channel. */
    Real a() const { return c[3]; }

    /** A reference to the alpha channel. */
    Real & a() { return c[3]; }

    /** Get the red, green and blue channels as a ColorRGB. */
    ColorRGB rgb() const { return ColorRGB(c[0], c[1], c[2]); }

    /** Set all channels simultaneously. */
    void set(Real r_, Real g_, Real b_, Real a_)
    {
      c[0] = r_;
      c[1] = g_;
      c[2] = b_;
      c[3] = a_;
    }

    /** Addition. */
    ColorRGBA operator+(ColorRGBA const & rhs) const
    {
      return ColorRGBA(c[0] + rhs.c[0], c[1] + rhs.c[1], c[2] + rhs.c[2], c[3] + rhs.c[3]);
    }

    /** Subtraction. */
    ColorRGBA operator-(ColorRGBA const & rhs) const
    {
      return ColorRGBA(c[0] - rhs.c[0], c[1] - rhs.c[1], c[2] - rhs.c[2], c[3] - rhs.c[3]);
    }

    /** Multiplication by a scalar. */
    ColorRGBA operator*(Real s) const
    {
      return ColorRGBA(c[0] * s, c[1] * s, c[2] * s, c[3] * s);
    }

    /** Component-wise multiplication by another color. */
    ColorRGBA operator*(ColorRGBA const & rhs) const
    {
      return ColorRGBA(c[0] * rhs.c[0], c[1] * rhs.c[1], c[2] * rhs.c[2], c[3] * rhs.c[3]);
    }

    /** Division by a scalar. */
    ColorRGBA operator/(Real s) const
    {
      return ColorRGBA(c[0] / s, c[1] / s, c[2] / s, c[3] / s);
    }

    /** Component-wise division by another color. */
    ColorRGBA operator/(ColorRGBA const & rhs) const
    {
      return ColorRGBA(c[0] / rhs.c[0], c[1] / rhs.c[1], c[2] / rhs.c[2], c[3] / rhs.c[3]);
    }

    /** Negation of the color. */
    ColorRGBA operator-() const
    {
      return ColorRGBA(-c[0], -c[1], -c[2], -c[3]);
    }

    /** Add and assign. */
    ColorRGBA & operator+=(ColorRGBA const & rhs)
    {
      c[0] += rhs.c[0];
      c[1] += rhs.c[1];
      c[2] += rhs.c[2];
      c[3] += rhs.c[3];
      return *this;
    }

    /** Subtract and assign. */
    ColorRGBA & operator-=(ColorRGBA const & rhs)
    {
      c[0] -= rhs.c[0];
      c[1] -= rhs.c[1];
      c[2] -= rhs.c[2];
      c[3] -= rhs.c[3];
      return *this;
    }

    /** Multiply component-wise and assign. */
    ColorRGBA & operator*=(ColorRGBA const & rhs)
    {
      c[0] *= rhs.c[0];
      c[1] *= rhs.c[1];
      c[2] *= rhs.c[2];
      c[3] *= rhs.c[3];
      return *this;
    }

    /** Divide component-wise and assign. */
    ColorRGBA & operator/=(ColorRGBA const & rhs)
    {
      c[0] /= rhs.c[0];
      c[1] /= rhs.c[1];
      c[2] /= rhs.c[2];
      c[3] /= rhs.c[3];
      return *this;
    }

    /** Multiply by a scalar and assign. */
    ColorRGBA & operator*=(Real s)
    {
      c[0] *= s;
      c[1] *= s;
      c[2] *= s;
      c[3] *= s;
      return *this;
    }

    /** Divide by a scalar and assign. */
    ColorRGBA & operator/=(Real s)
    {
      c[0] /= s;
      c[1] /= s;
      c[2] /= s;
      c[3] /= s;
      return *this;
    }

    /** Check if two colors are approximately equal. */
    bool fuzzyEq(ColorRGBA const & other) const
    {
      return Math::fuzzyEq((*this - other).squaredLength(), (Real)0);
    }

    /** Check if two colors are not approximately equal. */
    bool fuzzyNe(ColorRGBA const & other) const
    {
      return Math::fuzzyNe((*this - other).squaredLength(), (Real)0);
    }

    /** Swap the red and blue channels. */
    ColorRGBA bgra() const
    {
      return ColorRGBA(c[2], c[1], c[0], c[3]);
    }

    /** Get a string representation of the color. */
    std::string toString() const;

    /** Color with all channels zero (transparent black). */
    static ColorRGBA const & zero();

  private:
    /** Get the square of the magnitude of the color. */
    Real squaredLength() const { return c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + c[3] * c[3]; }

}; // class ColorRGBA

/** Multiply by a scalar. */
inline ColorRGBA
operator*(Real s, ColorRGBA const & c)
{
  return c * s;
}

/** Multiply by a one-channel color. */
inline ColorRGBA
operator*(ColorL & s, ColorRGBA const & c)
{
  return c * s.value();
}

/** Multiply by a one-channel color. */
inline ColorRGBA
operator*(ColorRGBA const & c, ColorL & s)
{
  return c * s.value();
}

} // namespace Thea

#endif
