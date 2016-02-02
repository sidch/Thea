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
class ColorRGBA;

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

    /** Copy constructor. */
    ColorL(ColorL const & other) : val(other.val) {}

    /** Initialize from an integer color, automatically dividing by 255. */
    ColorL(ColorL8 const & other);

    /** Initialize from a 32-bit RGBA color. For conversion from a consistent source type. */
    ColorL(ColorRGBA const & other);

    /** The value of the color. */
    Real value() const { return val; }

    /** The value of the color. */
    Real & value() { return val; }

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
