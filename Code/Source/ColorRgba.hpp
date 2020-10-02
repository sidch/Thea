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

#ifndef __Thea_ColorRgba_hpp__
#define __Thea_ColorRgba_hpp__

#include "Common.hpp"
#include "ColorL.hpp"
#include "ColorRgb.hpp"
#include "Math.hpp"
#include "MatVec.hpp"

namespace Thea {

// Forward declarations
class ColorRgb8;
class ColorRgba8;

/**
 * A color with four floating-point channels: red, green and blue, each in [0, 1]. Derived from the G3D library:
 * http://g3d.sourceforge.net
 */
class THEA_API ColorRgba
{
  private:
    Real c[4];  ///< Four components: red, gree, blue and alpha.

  public:
    /** Default constructor. Does not initialize fields. */
    ColorRgba() {}

    /** Construct from red, green, blue and alpha components. */
    ColorRgba(Real r_, Real g_, Real b_, Real a_ = 1) { c[0] = r_; c[1] = g_; c[2] = b_; c[3] = a_; }

    /** Construct from an RGB color and an alpha component. */
    ColorRgba(ColorRgb const & rgb_, Real a_ = 1) { c[0] = rgb_.r(); c[1] = rgb_.g(); c[2] = rgb_.b(); c[3] = a_; }

    /** Construct a color from a 4-vector. */
    explicit ColorRgba(Vector4 const & v) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Construct a color from four components in an array. */
    explicit ColorRgba(Real const * v) { c[0] = v[0]; c[1] = v[1]; c[2] = v[2]; c[3] = v[3]; }

    /** Construct from a color with 4 byte channels, with automatic scaling from [0, 255] to [0, 1]. */
    ColorRgba(ColorRgba8 const & src);

    /** Construct from a color with 3 byte channels, with automatic scaling from [0, 255] to [0, 1]. */
    ColorRgba(ColorRgb8 const & src);

    /** Initialize from an HTML-style color (e.g. 0xFFFF0000 == RED) */
    static ColorRgba fromARGB(uint32 argb);

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

    /** Get the red, green and blue channels as a ColorRgb. */
    ColorRgb rgb() const { return ColorRgb(c[0], c[1], c[2]); }

    /** Get the address of the array storing color channel values in RGBA order. */
    Real const * data() const { return c; }

    /** Get the address of the array storing color channel values in RGBA order. */
    Real * data() { return c; }

    /** Array-style channel access. */
    template <typename IntegerT> Real const & operator[](IntegerT channel) const
    {
      debugAssertM(channel >= 0 && channel <= 3, "ColorRgba: Channel must be 0, 1, 2 or 3");
      return c[channel];
    }

    /** Array-style channel access. */
    template <typename IntegerT> Real & operator[](IntegerT channel)
    {
      debugAssertM(channel >= 0 && channel <= 3, "ColorRgba: Channel must be 0, 1, 2 or 3");
      return c[channel];
    }

    /** Set all channels simultaneously. */
    void set(Real r_, Real g_, Real b_, Real a_)
    {
      c[0] = r_;
      c[1] = g_;
      c[2] = b_;
      c[3] = a_;
    }

    /** Addition. */
    ColorRgba operator+(ColorRgba const & rhs) const
    {
      return ColorRgba(c[0] + rhs.c[0], c[1] + rhs.c[1], c[2] + rhs.c[2], c[3] + rhs.c[3]);
    }

    /** Subtraction. */
    ColorRgba operator-(ColorRgba const & rhs) const
    {
      return ColorRgba(c[0] - rhs.c[0], c[1] - rhs.c[1], c[2] - rhs.c[2], c[3] - rhs.c[3]);
    }

    /** Multiplication by a scalar. */
    ColorRgba operator*(Real s) const
    {
      return ColorRgba(c[0] * s, c[1] * s, c[2] * s, c[3] * s);
    }

    /** Component-wise multiplication by another color. */
    ColorRgba operator*(ColorRgba const & rhs) const
    {
      return ColorRgba(c[0] * rhs.c[0], c[1] * rhs.c[1], c[2] * rhs.c[2], c[3] * rhs.c[3]);
    }

    /** Division by a scalar. */
    ColorRgba operator/(Real s) const
    {
      return ColorRgba(c[0] / s, c[1] / s, c[2] / s, c[3] / s);
    }

    /** Component-wise division by another color. */
    ColorRgba operator/(ColorRgba const & rhs) const
    {
      return ColorRgba(c[0] / rhs.c[0], c[1] / rhs.c[1], c[2] / rhs.c[2], c[3] / rhs.c[3]);
    }

    /** Negation of the color. */
    ColorRgba operator-() const
    {
      return ColorRgba(-c[0], -c[1], -c[2], -c[3]);
    }

    /** Add and assign. */
    ColorRgba & operator+=(ColorRgba const & rhs)
    {
      c[0] += rhs.c[0];
      c[1] += rhs.c[1];
      c[2] += rhs.c[2];
      c[3] += rhs.c[3];
      return *this;
    }

    /** Subtract and assign. */
    ColorRgba & operator-=(ColorRgba const & rhs)
    {
      c[0] -= rhs.c[0];
      c[1] -= rhs.c[1];
      c[2] -= rhs.c[2];
      c[3] -= rhs.c[3];
      return *this;
    }

    /** Multiply component-wise and assign. */
    ColorRgba & operator*=(ColorRgba const & rhs)
    {
      c[0] *= rhs.c[0];
      c[1] *= rhs.c[1];
      c[2] *= rhs.c[2];
      c[3] *= rhs.c[3];
      return *this;
    }

    /** Divide component-wise and assign. */
    ColorRgba & operator/=(ColorRgba const & rhs)
    {
      c[0] /= rhs.c[0];
      c[1] /= rhs.c[1];
      c[2] /= rhs.c[2];
      c[3] /= rhs.c[3];
      return *this;
    }

    /** Multiply by a scalar and assign. */
    ColorRgba & operator*=(Real s)
    {
      c[0] *= s;
      c[1] *= s;
      c[2] *= s;
      c[3] *= s;
      return *this;
    }

    /** Divide by a scalar and assign. */
    ColorRgba & operator/=(Real s)
    {
      c[0] /= s;
      c[1] /= s;
      c[2] /= s;
      c[3] /= s;
      return *this;
    }

    /** Check if two colors are approximately equal. */
    bool fuzzyEq(ColorRgba const & other) const
    {
      return Math::fuzzyEq((*this - other).squaredNorm(), (Real)0);
    }

    /** Check if two colors are not approximately equal. */
    bool fuzzyNe(ColorRgba const & other) const
    {
      return Math::fuzzyNe((*this - other).squaredNorm(), (Real)0);
    }

    /** Swap the red and blue channels. */
    ColorRgba bgra() const
    {
      return ColorRgba(c[2], c[1], c[0], c[3]);
    }

    /** Get a string representation of the color. */
    std::string toString() const;

    // Named colors
    static ColorRgba const & red();     ///< Opaque red color.
    static ColorRgba const & green();   ///< Opaque green color.
    static ColorRgba const & blue();    ///< Opaque blue color.
    static ColorRgba const & purple();  ///< Opaque purple color.
    static ColorRgba const & cyan();    ///< Opaque cyan color.
    static ColorRgba const & yellow();  ///< Opaque yellow color.
    static ColorRgba const & brown();   ///< Opaque brown color.
    static ColorRgba const & orange();  ///< Opaque orange color.
    static ColorRgba const & black();   ///< Opaque black color.
    static ColorRgba const & gray();    ///< Opaque gray color.
    static ColorRgba const & white();   ///< Opaque white color.
    static ColorRgba const & zero();    ///< Color with all channels zero (transparent black).

  private:
    /** Get the square of the magnitude of the color. */
    Real squaredNorm() const { return c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + c[3] * c[3]; }

}; // class ColorRgba

/** Multiply by a scalar. */
inline ColorRgba
operator*(Real s, ColorRgba const & c)
{
  return c * s;
}

/** Multiply by a one-channel color. */
inline ColorRgba
operator*(ColorL & s, ColorRgba const & c)
{
  return c * s.value();
}

/** Multiply by a one-channel color. */
inline ColorRgba
operator*(ColorRgba const & c, ColorL & s)
{
  return c * s.value();
}

} // namespace Thea

#endif
