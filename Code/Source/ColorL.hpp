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

 @file ColorL.h

 Monochrome Color class

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu
 @created 2007-01-31
 @edited  2009-03-20

 Copyright 2000-2009, Morgan McGuire.
 All rights reserved.
 */

#ifndef __Thea_ColorL_hpp__
#define __Thea_ColorL_hpp__

#include "Common.hpp"
#include "Math.hpp"

namespace Thea {

// Forward declarations
class ColorL8;
class ColorRgba;

/**
 * Monochrome luminance value in [0, 1], with automatic scaling by 255 when switching between integer (ColorL8) and floating
 * point (ColorL) formats. Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API ColorL
{
  private:
    Real val;  ///< Luminance value.

  public:
    /** Default constructor, initializes color to 0. */
    ColorL() : val(0) {}

    /** Initializing constructor. */
    explicit ColorL(Real v) : val(v) {}

    /** Initialize from an integer color, automatically dividing by 255. */
    ColorL(ColorL8 const & other);

    /** Initialize from a 32-bit RGBA color. For conversion from a consistent source type. */
    ColorL(ColorRgba const & other);

    /** The value of the color. */
    Real value() const { return val; }

    /** The value of the color. */
    Real & value() { return val; }

    /** Get the address of the color value. */
    Real const * data() const { return &val; }

    /** Get the address of the color value. */
    Real * data() { return &val; }

    /** Array-style channel access. */
    template <typename IntegerT> Real const & operator[](IntegerT channel) const
    {
      debugAssertM(channel == 0, "ColorL: Channel must be 0");
      return val;
    }

    /** Array-style channel access. */
    template <typename IntegerT> Real & operator[](IntegerT channel)
    {
      debugAssertM(channel == 0, "ColorL: Channel must be 0");
      return val;
    }

    /** Negation operator. */
    ColorL operator-() const
    {
      return ColorL(-val);
    }

    /** Addition operator. */
    ColorL operator+(ColorL const & other) const
    {
      return ColorL(val + other.val);
    }

    /** Subtraction operator. */
    ColorL operator-(ColorL const & other) const
    {
      return ColorL(val - other.val);
    }

    /** Multiplication operator. */
    ColorL operator*(ColorL const & other) const
    {
      return ColorL(val * other.val);
    }

    /** Multiply by a scalar. */
    ColorL operator*(Real other) const
    {
      return ColorL(val * other);
    }

    /** Division operator. */
    ColorL operator/(ColorL const & other) const
    {
      return ColorL(val / other.val);
    }

    /** Divide by a scalar. */
    ColorL operator/(Real other) const
    {
      return ColorL(val / other);
    }

    /** Add-and-assign. */
    ColorL & operator+=(ColorL const & other)
    {
      val += other.val;
      return *this;
    }

    /** Subtract-and-assign. */
    ColorL & operator-=(ColorL const & other)
    {
      val -= other.val;
      return *this;
    }

    /** Multiply-and-assign. */
    ColorL & operator*=(ColorL const & other)
    {
      val *= other.val;
      return *this;
    }

    /** Divide-and-assign. */
    ColorL & operator/=(ColorL const & other)
    {
      val /= other.val;
      return *this;
    }

    /** Multiply by a scalar and assign. */
    ColorL & operator*=(Real other)
    {
      val *= other;
      return *this;
    }

    /** Divide by a scalar and assign. */
    ColorL & operator/=(Real other)
    {
      val /= other;
      return *this;
    }

    /** Get the maximum of two color values. */
    ColorL max(ColorL const & other) const
    {
      return ColorL(std::max(val, other.val));
    }

    /** Get the minimum of two color values. */
    ColorL min(ColorL const & other) const
    {
      return ColorL(std::min(val, other.val));
    }

    /** Color with zero luminance (black). */
    static ColorL const & zero() { static ColorL const col(0); return col; }

    /** Get a string representation of the color. */
    std::string toString() const;

}; // class ColorL

/** Multiply by a scalar. */
inline ColorL
operator*(Real lhs, ColorL const & rhs)
{
  return rhs * lhs;
}

} // namespace Thea

#endif
