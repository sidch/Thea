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

 @file ImageFormat.cpp

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2003-05-23
 @edited  2010-03-30
*/

#include "TextureFormat.hpp"
#include "../Plugins/GL/glew.h"

// Let's hope the EXT versions exist at least... else we're on really really old GL
#ifndef GL_BGR
#  define GL_BGR GL_BGR_EXT
#endif

#ifndef GL_BGRA
#  define GL_BGRA GL_BGRA_EXT
#endif

namespace Thea {
namespace Graphics {

TextureFormat::TextureFormat(
  int             _numComponents,
  bool            _compressed,
  int             _glFormat,
  int             _glBaseFormat,
  int             _luminanceBits,
  int             _alphaBits,
  int             _redBits,
  int             _greenBits,
  int             _blueBits,
  int             _depthBits,
  int             _stencilBits,
  int             _glBitsPerPixel,
  int             _cpuBitsPerPixel,
  int             _glDataFormat,
  bool            _opaque,
  bool            _floatingPoint,
  Code            _code,
  ColorSpace      _colorSpace,
  BayerPattern    _bayerPattern)
:
  numComponents(_numComponents),
  compressed(_compressed),
  code(_code),
  colorSpace(_colorSpace),
  bayerPattern(_bayerPattern),
  openGLFormat(_glFormat),
  openGLBaseFormat(_glBaseFormat),
  luminanceBits(_luminanceBits),
  alphaBits(_alphaBits),
  redBits(_redBits),
  greenBits(_greenBits),
  blueBits(_blueBits),
  stencilBits(_stencilBits),
  depthBits(_depthBits),
  cpuBitsPerPixel(_cpuBitsPerPixel),
  openGLBitsPerPixel(_glBitsPerPixel),
  openGLDataFormat(_glDataFormat),
  opaque(_opaque),
  floatingPoint(_floatingPoint)
{
  debugAssertM(_cpuBitsPerPixel <= _glBitsPerPixel, "TextureFormat: Too many packed bits");
}

bool
TextureFormat::isDepth() const
{
  return depthBits > 0;
}

TextureFormat const *
TextureFormat::depth(int depth_bits)
{
  switch (depth_bits)
  {
    case 16:
      return DEPTH16();

    case 24:
      return DEPTH24();

    case 32:
      return DEPTH32();

    default:
      alwaysAssertM(false, "TextureFormat: Number of depth bits must be 16, 24, or 32");
      return DEPTH32();
  }
}

TextureFormat const *
TextureFormat::stencil(int bits)
{
  switch (bits)
  {
    case 1:
      return STENCIL1();

    case 4:
      return STENCIL4();

    case 8:
      return STENCIL8();

    case 16:
      return STENCIL16();

    default:
      alwaysAssertM(false, "TextureFormat: Number of stencil bits must be 1, 4, 8 or 16");
      return STENCIL16();
  }
}

namespace TextureFormatInternal {

static std::string const nameArray[] =
{
  "L8",
  "L16",
  "L16F",
  "L32F",

  "A8",
  "A16",
  "A16F",
  "A32F",

  "LA4",
  "LA8",
  "LA16",
  "LA16F",
  "LA32F",

  "RGB5",
  "RGB5A1",
  "RGB8",
  "RGB10",
  "RGB10A2",
  "RGB16",
  "RGB16F",
  "RGB32F",
  "R11G11B10F",
  "RGB9E10F",

  "RGB8I",
  "RGB8UI",

  "RGBA8UI",

  "R8",

  "RG8",
  "RG8I",
  "RG8UI",

  "RG16F",

  "RGBA8",
  "RGBA16",
  "RGBA16F",
  "RGBA32F",

  "RGBA32UI",

  "BGR8",
  "BGRA8",
  "BGR16",
  "BGRA16",
  "BGR32F",
  "BGRA32F",
  "ARGB8",

  "BAYER_RGGB8",
  "BAYER_GRBG8",
  "BAYER_GBRG8",
  "BAYER_BGGR8",
  "BAYER_RGGB32F",
  "BAYER_GRBG32F",
  "BAYER_GBRG32F",
  "BAYER_BGGR32F",

  "HSV8",
  "HSV32F",

  "YUV420_PLANAR",
  "YUV422",
  "YUV444",

  "RGB_DXT1",
  "RGBA_DXT1",
  "RGBA_DXT3",
  "RGBA_DXT5",

  "SRGB8",
  "SRGBA8",

  "SL8",
  "SLA8",

  "SRGB_DXT1",
  "SRGBA_DXT1",
  "SRGBA_DXT3",
  "SRGBA_DXT5",

  "DEPTH16",
  "DEPTH24",
  "DEPTH32",
  "DEPTH32F",

  "STENCIL1",
  "STENCIL4",
  "STENCIL8",
  "STENCIL16",

  "DEPTH24_STENCIL8",
  ""
};

} // namespace TextureFormatInternal

std::string const &
TextureFormat::name() const
{
  debugAssertM(code < Code::NUM, "TextureFormat: Invalid code");
  return TextureFormatInternal::nameArray[code];
}

TextureFormat const *
TextureFormat::fromString(std::string const & s)
{
  for (int i = 0; i < Code::NUM; ++i)
  {
    if (s == TextureFormatInternal::nameArray[i])
    {
      return fromCode(Code(i));
    }
  }

  return NULL;
}

TextureFormat const *
TextureFormat::fromCode(Code code)
{
  switch (code)
  {
    case Code::L8:
      return L8();

    case Code::L16:
      return L16();

    case Code::L16F:
      return L16F();

    case Code::L32F:
      return L32F();

    case Code::A8:
      return A8();

    case Code::A16:
      return A16();

    case Code::A16F:
      return A16F();

    case Code::A32F:
      return A32F();

    case Code::LA4:
      return LA4();

    case Code::LA8:
      return LA8();

    case Code::LA16:
      return LA16();

    case Code::LA16F:
      return LA16F();
      break;

    case Code::LA32F:
      return LA32F();

    case Code::RGB5:
      return RGB5();

    case Code::RGB5A1:
      return RGB5A1();

    case Code::RGB8:
      return RGB8();

    case Code::RGB10:
      return RGB10();

    case Code::RGB10A2:
      return RGB10A2();

    case Code::RGB16:
      return RGB16();

    case Code::RGB32F:
      return RGB32F();

    case Code::R11G11B10F:
      return R11G11B10F();

    case Code::RGB9E5F:
      return RGB9E5F();

    case Code::RGB8I:
      return RGB8I();

    case Code::RGB8UI:
      return RGB8UI();

    case Code::ARGB8:
      return NULL;

    case Code::BGR8:
      return BGR8();

    case Code::BGRA8:
      return BGRA8();

    case Code::BGR16:
      return BGR16();

    case Code::BGRA16:
      return BGRA16();

    case Code::BGR32F:
      return BGR32F();

    case Code::BGRA32F:
      return BGRA32F();

    case Code::R8:
      return R8();

    case Code::RG8:
      return RG8();

    case Code::RG8I:
      return RG8I();

    case Code::RG8UI:
      return RG8UI();

    case Code::RG16F:
      return RG16F();

    case Code::RGBA8:
      return RGBA8();

    case Code::RGBA16:
      return RGBA16();

    case Code::RGBA16F:
      return RGBA16F();

    case Code::RGBA32F:
      return RGBA32F();

    case Code::RGBA32UI:
      return RGBA32UI();

    case Code::BAYER_RGGB8:

      // TODO
    case Code::BAYER_GRBG8:

      // TODO
    case Code::BAYER_GBRG8:

      // TODO
    case Code::BAYER_BGGR8:

      // TODO
    case Code::BAYER_RGGB32F:

      // TODO
    case Code::BAYER_GRBG32F:

      // TODO
    case Code::BAYER_GBRG32F:

      // TODO
    case Code::BAYER_BGGR32F:

      // TODO
    case Code::HSV8:

      // TODO
    case Code::HSV32F:
      // TODO
      return NULL;
      break;

    case Code::RGB_DXT1:
      return RGB_DXT1();
      break;

    case Code::RGBA_DXT1:
      return RGBA_DXT1();
      break;

    case Code::RGBA_DXT3:
      return RGBA_DXT3();
      break;

    case Code::RGBA_DXT5:
      return RGBA_DXT5();
      break;

    case Code::SRGB8:
      return SRGB8();
      break;

    case Code::SRGBA8:
      return SRGBA8();
      break;

    case Code::SL8:
      return SL8();
      break;

    case Code::SLA8:
      return SLA8();
      break;

    case Code::SRGB_DXT1:
      return SRGB_DXT1();
      break;

    case Code::SRGBA_DXT1:
      return SRGBA_DXT1();
      break;

    case Code::SRGBA_DXT3:
      return SRGBA_DXT3();
      break;

    case Code::SRGBA_DXT5:
      return SRGBA_DXT5();
      break;

    case Code::DEPTH16:
      return DEPTH16();
      break;

    case Code::DEPTH24:
      return DEPTH24();
      break;

    case Code::DEPTH32:
      return DEPTH32();
      break;

    case Code::DEPTH32F:
      return DEPTH32F();
      break;

    case Code::STENCIL1:
      return STENCIL1();
      break;

    case Code::STENCIL4:
      return STENCIL4();
      break;

    case Code::STENCIL8:
      return STENCIL8();
      break;

    case Code::STENCIL16:
      return STENCIL16();
      break;

    case Code::DEPTH24_STENCIL8:
      return DEPTH24_STENCIL8();
      break;

    case Code::YUV420_PLANAR:
      return YUV420_PLANAR();
      break;

    case Code::YUV422:
      return YUV422();
      break;

    case Code::YUV444:
      return YUV444();
      break;

    default:
      return NULL;
  }
}

TextureFormat const *
TextureFormat::fromImageType(AbstractImage::Type type, bool is_depth)
{
  if (is_depth)
  {
    switch (type)
    {
      case AbstractImage::Type::LUMINANCE_16U : return DEPTH16();
      case AbstractImage::Type::LUMINANCE_32U : return DEPTH32();
      case AbstractImage::Type::LUMINANCE_32F : return DEPTH32F();
      default: throw Error("TextureFormat: No supported depth texture format corresponds to the specified image format");
    }
  }

  enum { COLOR_ORDER_RGB, COLOR_ORDER_BGR } color_order;
  color_order = (AbstractImage::Channel::RED < AbstractImage::Channel::BLUE ? COLOR_ORDER_RGB : COLOR_ORDER_BGR);

  switch (type)
  {
    case AbstractImage::Type::LUMINANCE_8U  : return L8();
    case AbstractImage::Type::LUMINANCE_16U : return L16();
    case AbstractImage::Type::LUMINANCE_32F : return L32F();
    case AbstractImage::Type::RGB_8U        : return color_order == COLOR_ORDER_RGB ? RGB8() : BGR8();
    case AbstractImage::Type::RGBA_8U       : return color_order == COLOR_ORDER_RGB ? RGBA8() : BGRA8();
    case AbstractImage::Type::RGB_16U       : return color_order == COLOR_ORDER_RGB ? RGB16() : BGR16();
    case AbstractImage::Type::RGBA_16U      : return color_order == COLOR_ORDER_RGB ? RGBA16() : BGRA16();
    case AbstractImage::Type::RGB_32F       : return color_order == COLOR_ORDER_RGB ? RGB32F() : BGR32F();
    case AbstractImage::Type::RGBA_32F      : return color_order == COLOR_ORDER_RGB ? RGBA32F() : BGRA32F();
    default: throw Error("TextureFormat: No supported texture format corresponds to the specified image format");
  }
}

// Helper variables for defining texture formats

// Is floating point format
static bool const FLOAT_FORMAT  = true;
static bool const INT_FORMAT    = false;

// Is opaque format (no alpha)
static bool const OPAQUE_FORMAT = true;
static bool const CLEAR_FORMAT  = false;

// Is compressed format (not raw component data)
static bool const COMP_FORMAT   = true;
static bool const UNCOMP_FORMAT = false;

#define THEA_TEXTURE_FORMAT(enumname, cmpnts, cmprssd, glf, glbf, lb, ab, rb, gb, bb, db, sb, gbpp, cbpp, gldf, opq, fp, code, cs) \
  TextureFormat const * TextureFormat::enumname() { \
    static TextureFormat const format(cmpnts, cmprssd, glf, glbf, lb, ab, rb, gb, bb, db, sb, gbpp, cbpp, gldf, opq, fp, TextureFormat::Code::code, TextureFormat::ColorSpace::cs); \
    return &format; }

THEA_TEXTURE_FORMAT(L8,         1, UNCOMP_FORMAT,   GL_LUMINANCE8,                           GL_LUMINANCE,          8,  0,  0,  0,  0,  0,  0, 8, 8, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, L8, NONE);

THEA_TEXTURE_FORMAT(L16,        1, UNCOMP_FORMAT,   GL_LUMINANCE16,                          GL_LUMINANCE,         16,  0,  0,  0,  0,  0,  0, 16, 16, GL_UNSIGNED_SHORT, OPAQUE_FORMAT, INT_FORMAT, L16, NONE);

THEA_TEXTURE_FORMAT(L16F,       1, UNCOMP_FORMAT,   GL_LUMINANCE16F_ARB,                     GL_LUMINANCE,         16,  0,  0,  0,  0,  0,  0, 16, 16, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, L16F, NONE);

THEA_TEXTURE_FORMAT(L32F,       1, UNCOMP_FORMAT,   GL_LUMINANCE32F_ARB,                     GL_LUMINANCE,         32,  0,  0,  0,  0,  0,  0, 32, 32, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, L32F, NONE);

THEA_TEXTURE_FORMAT(A8,         1, UNCOMP_FORMAT,   GL_ALPHA8,                               GL_ALPHA,              0,  8,  0,  0,  0,  0,  0,  8, 8, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, A8, NONE);

THEA_TEXTURE_FORMAT(A16,        1, UNCOMP_FORMAT,   GL_ALPHA16,                              GL_ALPHA,              0, 16,  0,  0,  0,  0,  0, 16, 16, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, A16, NONE);

THEA_TEXTURE_FORMAT(A16F,       1, UNCOMP_FORMAT,   GL_ALPHA16F_ARB,                         GL_ALPHA,              0, 16,  0,  0,  0,  0,  0, 16, 16, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, A16F, NONE);

THEA_TEXTURE_FORMAT(A32F,       1, UNCOMP_FORMAT,   GL_ALPHA32F_ARB,                         GL_ALPHA,              0, 32,  0,  0,  0,  0,  0, 32, 32, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, A32F, NONE);

THEA_TEXTURE_FORMAT(LA4,        2, UNCOMP_FORMAT,   GL_LUMINANCE4_ALPHA4,                    GL_LUMINANCE_ALPHA,    4,  4,  0,  0,  0,  0,  0,  8, 8, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, LA4, NONE);

THEA_TEXTURE_FORMAT(LA8,        2, UNCOMP_FORMAT,   GL_LUMINANCE8_ALPHA8,                    GL_LUMINANCE_ALPHA,    8,  8,  0,  0,  0,  0,  0, 16, 16, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, LA8, NONE);

THEA_TEXTURE_FORMAT(LA16,       2, UNCOMP_FORMAT,   GL_LUMINANCE16_ALPHA16,                  GL_LUMINANCE_ALPHA,   16, 16,  0,  0,  0,  0,  0, 16 * 2, 16 * 2, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, LA16, NONE);

THEA_TEXTURE_FORMAT(LA16F,      2, UNCOMP_FORMAT,   GL_LUMINANCE_ALPHA16F_ARB,               GL_LUMINANCE_ALPHA,   16, 16,  0,  0,  0,  0,  0, 16 * 2, 16 * 2, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, LA16F, NONE);

THEA_TEXTURE_FORMAT(LA32F,      2, UNCOMP_FORMAT,   GL_LUMINANCE_ALPHA32F_ARB,               GL_LUMINANCE_ALPHA,   32, 32,  0,  0,  0,  0,  0, 32 * 2, 32 * 2, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, LA32F, NONE);

THEA_TEXTURE_FORMAT(R8,         1, UNCOMP_FORMAT,   GL_R8,                                   GL_RED,                0,  0,  8,  0,  0,  0,  0,  8, 8, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, R8, RGB);

THEA_TEXTURE_FORMAT(RG8,        2, UNCOMP_FORMAT,   GL_RG8,                                  GL_RG,                 0,  0,  8,  8,  0,  0,  0, 16, 16, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RG8, RGB);

// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
THEA_TEXTURE_FORMAT(RG8I,       2, UNCOMP_FORMAT,   GL_RG8I,                                 GL_RG_INTEGER,         0,  0,  8,  8,  0,  0,  0, 16, 16, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RG8I, RGB);

THEA_TEXTURE_FORMAT(RG8UI,      2, UNCOMP_FORMAT,   GL_RG8UI,                                GL_RG_INTEGER,         0,  0,  8,  8,  0,  0,  0, 16, 16, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RG8UI, RGB);

THEA_TEXTURE_FORMAT(RG16F,      2, UNCOMP_FORMAT,   GL_RG16F,                                GL_RG,                 0,  0, 16, 16,  0,  0,  0, 32, 32, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, RG16F, RGB);

THEA_TEXTURE_FORMAT(RGB5,       3, UNCOMP_FORMAT,   GL_RGB5,                                 GL_RGBA,               0,  0,  5,  5,  5,  0,  0, 16, 16, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB5, RGB);

THEA_TEXTURE_FORMAT(RGB5A1,     4, UNCOMP_FORMAT,   GL_RGB5_A1,                              GL_RGBA,               0,  1,  5,  5,  5,  0,  0, 16, 16, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB5A1, RGB);

THEA_TEXTURE_FORMAT(RGB8,       3, UNCOMP_FORMAT,   GL_RGB8,                                 GL_RGB,                0,  0,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB8, RGB);

THEA_TEXTURE_FORMAT(RGB10,      3, UNCOMP_FORMAT,   GL_RGB10,                                GL_RGB,                0,  0, 10, 10, 10,  0,  0, 32, 10 * 3, GL_UNSIGNED_SHORT, OPAQUE_FORMAT, INT_FORMAT, RGB10, RGB);

THEA_TEXTURE_FORMAT(RGB10A2,    4, UNCOMP_FORMAT,   GL_RGB10_A2,                             GL_RGBA,               0,  2, 10, 10, 10,  0,  0, 32, 32, GL_UNSIGNED_INT_10_10_10_2, OPAQUE_FORMAT, INT_FORMAT, RGB10A2, RGB);

THEA_TEXTURE_FORMAT(RGB16,      3, UNCOMP_FORMAT,   GL_RGB16,                                GL_RGB,                0,  0, 16, 16, 16,  0,  0, 16 * 3, 16 * 3, GL_UNSIGNED_SHORT, OPAQUE_FORMAT, INT_FORMAT, RGB16, RGB);

THEA_TEXTURE_FORMAT(RGB16F,     3, UNCOMP_FORMAT,   GL_RGB16F_ARB,                           GL_RGB,                0,  0, 16, 16, 16,  0,  0, 16 * 3, 16 * 3, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB16F, RGB);

THEA_TEXTURE_FORMAT(RGB32F,     3, UNCOMP_FORMAT,   GL_RGB32F_ARB,                           GL_RGB,                0,  0, 32, 32, 32,  0,  0, 32 * 3, 32 * 3, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB32F, RGB);

THEA_TEXTURE_FORMAT(RGBA8,      4, UNCOMP_FORMAT,   GL_RGBA8,                                GL_RGBA,               0,  8,  8,  8,  8,  0,  0, 32, 32, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, RGBA8, RGB);

THEA_TEXTURE_FORMAT(RGBA16,     4, UNCOMP_FORMAT,   GL_RGBA16,                               GL_RGBA,               0, 16, 16, 16, 16,  0,  0, 16 * 4, 16 * 4, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, RGBA16, RGB);

THEA_TEXTURE_FORMAT(RGBA16F,    4, UNCOMP_FORMAT,   GL_RGBA16F_ARB,                          GL_RGBA,               0, 16, 16, 16, 16,  0,  0, 16 * 4, 16 * 4, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, RGBA16F, RGB);

THEA_TEXTURE_FORMAT(RGBA32F,    4, UNCOMP_FORMAT,   GL_RGBA32F_ARB,                          GL_RGBA,               0, 32, 32, 32, 32,  0,  0, 32 * 4, 32 * 4, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, RGBA32F, RGB);

// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
THEA_TEXTURE_FORMAT(RGBA32UI,   4, UNCOMP_FORMAT,   GL_RGBA32UI,                             GL_RGBA_INTEGER,       0, 32, 32, 32, 32,  0,  0, 32 * 4, 32 * 4, GL_UNSIGNED_INT, CLEAR_FORMAT, INT_FORMAT, RGBA32UI, RGB);

// Unsigned
THEA_TEXTURE_FORMAT(R11G11B10F, 3, UNCOMP_FORMAT,   GL_R11F_G11F_B10F_EXT,                   GL_RGB,                0,  0, 11, 11, 10,  0,  0, 32, 32, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, R11G11B10F, RGB);

// Unsigned
THEA_TEXTURE_FORMAT(RGB9E5F,    3, UNCOMP_FORMAT,   GL_RGB9_E5_EXT,                          GL_RGB,                0,  0, 14, 14, 14,  0,  0, 32, 32, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB9E5F, RGB);

// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
THEA_TEXTURE_FORMAT(RGB8I,      3, UNCOMP_FORMAT,   GL_RGB8I_EXT,                            GL_RGB_INTEGER,        0,  0,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB8I, RGB);

THEA_TEXTURE_FORMAT(RGB8UI,     3, UNCOMP_FORMAT,   GL_RGB8UI_EXT,                           GL_RGB_INTEGER,        0,  0,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB8UI, RGB);

THEA_TEXTURE_FORMAT(RGBA8UI,    4, UNCOMP_FORMAT,   GL_RGBA8UI_EXT,                          GL_RGBA_INTEGER,       0,  0,  8,  8,  8,  8,  0, 32, 32, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGBA8UI, RGB);


THEA_TEXTURE_FORMAT(BGR8,       3, UNCOMP_FORMAT,   GL_RGB8,                                 GL_BGR,                0,  0,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, BGR8, RGB);

THEA_TEXTURE_FORMAT(BGRA8,      4, UNCOMP_FORMAT,   GL_RGBA8,                                GL_BGRA,               0,  8,  8,  8,  8,  0,  0, 32, 32, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, BGRA8, RGB);

THEA_TEXTURE_FORMAT(BGR16,      3, UNCOMP_FORMAT,   GL_RGB16,                                GL_BGR,                0,  0, 16, 16, 16,  0,  0, 16 * 3, 16 * 3, GL_UNSIGNED_SHORT, OPAQUE_FORMAT, INT_FORMAT, BGR16, RGB);

THEA_TEXTURE_FORMAT(BGRA16,     4, UNCOMP_FORMAT,   GL_RGBA16,                               GL_BGRA,               0, 16, 16, 16, 16,  0,  0, 16 * 4, 16 * 4, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, BGRA16, RGB);

THEA_TEXTURE_FORMAT(BGR32F,     3, UNCOMP_FORMAT,   GL_RGB32F_ARB,                           GL_BGR,                0,  0, 32, 32, 32,  0,  0, 32 * 3, 32 * 3, GL_FLOAT, OPAQUE_FORMAT, FLOAT_FORMAT, BGR32F, RGB);

THEA_TEXTURE_FORMAT(BGRA32F,    4, UNCOMP_FORMAT,   GL_RGBA32F_ARB,                          GL_BGRA,               0, 32, 32, 32, 32,  0,  0, 32 * 4, 32 * 4, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, BGRA32F, RGB);


THEA_TEXTURE_FORMAT(RGB_DXT1,   3, COMP_FORMAT,     GL_COMPRESSED_RGB_S3TC_DXT1_EXT,         GL_RGB,                0,  0,  0,  0,  0,  0,  0, 64, 64, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, RGB_DXT1, RGB);

THEA_TEXTURE_FORMAT(RGBA_DXT1,  4, COMP_FORMAT,     GL_COMPRESSED_RGBA_S3TC_DXT1_EXT,        GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 64, 64, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, RGBA_DXT1, RGB);

THEA_TEXTURE_FORMAT(RGBA_DXT3,  4, COMP_FORMAT,     GL_COMPRESSED_RGBA_S3TC_DXT3_EXT,        GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 128, 128, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, RGBA_DXT3, RGB);

THEA_TEXTURE_FORMAT(RGBA_DXT5,  4, COMP_FORMAT,     GL_COMPRESSED_RGBA_S3TC_DXT5_EXT,        GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 128, 128, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, RGBA_DXT5, RGB);

THEA_TEXTURE_FORMAT(SRGB8,      3, UNCOMP_FORMAT,   GL_SRGB8,                                GL_RGB,                0,  0,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, SRGB8, SRGB);

THEA_TEXTURE_FORMAT(SRGBA8,     4, UNCOMP_FORMAT,   GL_SRGB8_ALPHA8,                         GL_RGBA,               0,  8,  8,  8,  8,  0,  0, 32, 24, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, SRGBA8, SRGB);

THEA_TEXTURE_FORMAT(SL8,        1, UNCOMP_FORMAT,   GL_SLUMINANCE8,                          GL_LUMINANCE,          8,  0,  0,  0,  0,  0,  0, 8, 8, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, SL8, SRGB);

THEA_TEXTURE_FORMAT(SLA8,       2, UNCOMP_FORMAT,   GL_SLUMINANCE8_ALPHA8,                   GL_LUMINANCE_ALPHA,    8,  8,  0,  0,  0,  0,  0, 16, 16, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, SLA8, SRGB);

THEA_TEXTURE_FORMAT(SRGB_DXT1,  3, COMP_FORMAT,     GL_COMPRESSED_SRGB_S3TC_DXT1_EXT,        GL_RGB,                0,  0,  0,  0,  0,  0,  0, 64, 64, GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, SRGB_DXT1, SRGB);

THEA_TEXTURE_FORMAT(SRGBA_DXT1, 4, COMP_FORMAT,     GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT,  GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 64, 64, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, SRGBA_DXT1, SRGB);

THEA_TEXTURE_FORMAT(SRGBA_DXT3, 4, COMP_FORMAT,     GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT,  GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 128, 128,  GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, SRGBA_DXT3, SRGB);

THEA_TEXTURE_FORMAT(SRGBA_DXT5, 4, COMP_FORMAT,     GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT,  GL_RGBA,               0,  0,  0,  0,  0,  0,  0, 128, 128,  GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, SRGBA_DXT5, SRGB);

THEA_TEXTURE_FORMAT(DEPTH16,    1, UNCOMP_FORMAT,   GL_DEPTH_COMPONENT16_ARB,                GL_DEPTH_COMPONENT,    0,  0,  0,  0,  0, 16, 0, 16, 16, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, DEPTH16, NONE);

THEA_TEXTURE_FORMAT(DEPTH24,    1, UNCOMP_FORMAT,   GL_DEPTH_COMPONENT24_ARB,                GL_DEPTH_COMPONENT,    0,  0,  0,  0,  0, 24, 0, 32, 24, GL_UNSIGNED_INT, CLEAR_FORMAT, INT_FORMAT, DEPTH24, NONE);

THEA_TEXTURE_FORMAT(DEPTH32,    1, UNCOMP_FORMAT,   GL_DEPTH_COMPONENT32_ARB,                GL_DEPTH_COMPONENT,    0,  0,  0,  0,  0, 32, 0, 32, 32, GL_UNSIGNED_INT, CLEAR_FORMAT, INT_FORMAT, DEPTH32, NONE);

THEA_TEXTURE_FORMAT(DEPTH32F,   1, UNCOMP_FORMAT,   GL_DEPTH_COMPONENT32_ARB,                GL_DEPTH_COMPONENT,    0,  0,  0,  0,  0, 32, 0, 32, 32, GL_FLOAT, CLEAR_FORMAT, FLOAT_FORMAT, DEPTH32F, NONE);

// These formats are for use with Renderbuffers only!
THEA_TEXTURE_FORMAT(STENCIL1,   1, UNCOMP_FORMAT,   GL_STENCIL_INDEX1_EXT,                   GL_STENCIL_INDEX,      0,  0,  0,  0,  0,  0,  1,  1,  1, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, STENCIL1, NONE);

THEA_TEXTURE_FORMAT(STENCIL4,   1, UNCOMP_FORMAT,   GL_STENCIL_INDEX4_EXT,                   GL_STENCIL_INDEX,      0,  0,  0,  0,  0,  0,  4,  4,  4, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, STENCIL4, NONE);

THEA_TEXTURE_FORMAT(STENCIL8,   1, UNCOMP_FORMAT,   GL_STENCIL_INDEX8_EXT,                   GL_STENCIL_INDEX,      0,  0,  0,  0,  0,  0,  8,  8,  8, GL_UNSIGNED_BYTE, CLEAR_FORMAT, INT_FORMAT, STENCIL8, NONE);

THEA_TEXTURE_FORMAT(STENCIL16,  1, UNCOMP_FORMAT,   GL_STENCIL_INDEX16_EXT,                  GL_STENCIL_INDEX,      0,  0,  0,  0,  0,  0, 16, 16, 16, GL_UNSIGNED_SHORT, CLEAR_FORMAT, INT_FORMAT, STENCIL16, NONE);

THEA_TEXTURE_FORMAT(DEPTH24_STENCIL8,   2, UNCOMP_FORMAT,   GL_DEPTH24_STENCIL8_EXT,         GL_DEPTH_STENCIL_EXT,  0,  0,  0,  0,  0, 24,  8, 32, 32, GL_UNSIGNED_INT, CLEAR_FORMAT, INT_FORMAT, DEPTH24_STENCIL8, NONE);

THEA_TEXTURE_FORMAT(YUV420_PLANAR,  3, UNCOMP_FORMAT,   GL_NONE,    GL_NONE,   0,  0,  0,  0,  0,  0,  0, 12, 12,  GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, YUV420_PLANAR, YUV);

THEA_TEXTURE_FORMAT(YUV422,         3, UNCOMP_FORMAT,   GL_NONE,    GL_NONE,   0,  0,  0,  0,  0,  0,  0, 16, 16,  GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, YUV422, YUV);

THEA_TEXTURE_FORMAT(YUV444,         3, UNCOMP_FORMAT,   GL_NONE,    GL_NONE,   0,  0,  0,  0,  0,  0,  0, 24, 24,  GL_UNSIGNED_BYTE, OPAQUE_FORMAT, INT_FORMAT, YUV444, YUV);

} // namespace Graphics
} // namespace Thea
