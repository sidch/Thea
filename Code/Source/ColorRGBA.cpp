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

 @file ColorRGBA.cpp

 Color class.

 @author Morgan McGuire, http://graphics.cs.williams.edu
 @cite Portions by Laura Wollstadt, graphics3d.com
 @cite Portions based on Dave Eberly's Magic Software Library at http://www.magic-software.com


 @created 2002-06-25
 @edited  2009-11-10
*/

#include "ColorRGBA.hpp"
#include "ColorRGBA8.hpp"

namespace Thea {

ColorRGBA const &
ColorRGBA::zero()
{
  static ColorRGBA const col(0, 0, 0, 0);
  return col;
}

ColorRGBA::ColorRGBA(ColorRGBA8 const & src)
{
  c[0] = src.r() / 255.0f;
  c[1] = src.g() / 255.0f;
  c[2] = src.b() / 255.0f;
  c[3] = src.a() / 255.0f;
}

ColorRGBA
ColorRGBA::fromARGB(uint32 x)
{
  return ColorRGBA((float)((x >> 16) & 0xFF),
                   (float)((x >>  8) & 0xFF),
                   (float)( x        & 0xFF),
                   (float)((x >> 24) & 0xFF)) / 255.0;
}

std::string
ColorRGBA::toString() const
{
  return format("RGBA(%g, %g, %g, %g)", c[0], c[1], c[2], c[3]);
}

} // namespace Thea
