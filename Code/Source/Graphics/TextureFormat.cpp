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
  int             num_components_,
  bool            is_compressed_,
  bool            is_opaque_,
  bool            is_floating_point_,
  Code            code_,
  ColorSpace      color_space_,
  BayerPattern    bayer_pattern_,
  int             luminance_bits_,
  int             red_bits_,
  int             green_bits_,
  int             blue_bits_,
  int             alpha_bits_,
  int             stencil_bits_,
  int             depth_bits_,
  int             cpu_bits_per_pixel_,
  int             gl_format_,
  int             gl_base_format_,
  int             gl_bits_per_pixel_,
  int             gl_data_format_
)
: num_components(num_components_),
  is_compressed(is_compressed_),
  is_opaque(is_opaque_),
  is_floating_point(is_floating_point_),
  code_val(code_),
  color_space(color_space_),
  bayer_pattern(bayer_pattern_),
  luminance_bits(luminance_bits_),
  red_bits(red_bits_),
  green_bits(green_bits_),
  blue_bits(blue_bits_),
  alpha_bits(alpha_bits_),
  stencil_bits(stencil_bits_),
  depth_bits(depth_bits_),
  cpu_bits_per_pixel(cpu_bits_per_pixel_),
  gl_format(gl_format_),
  gl_base_format(gl_base_format_),
  gl_bits_per_pixel(gl_bits_per_pixel_),
  gl_data_format(gl_data_format_)
{
  debugAssertM(cpu_bits_per_pixel_ <= gl_bits_per_pixel_, "TextureFormat: Too many packed bits");
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

char const *
TextureFormat::getName() const
{
  debugAssertM(code_val < Code::NUM, "TextureFormat: Invalid code");
  return TextureFormatInternal::nameArray[code_val].c_str();
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

TextureFormat const *
TextureFormat::fromCode(Code code_)
{
  switch (code_)
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
      return nullptr;

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
      return nullptr;
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
      return nullptr;
  }
}

TextureFormat const *
TextureFormat::fromImageType(IImage::Type type, bool is_depth)
{
  if (is_depth)
  {
    switch (type)
    {
      case IImage::Type::LUMINANCE_16U : return DEPTH16();
      case IImage::Type::LUMINANCE_32U : return DEPTH32();
      case IImage::Type::LUMINANCE_32F : return DEPTH32F();
      default: throw Error("TextureFormat: No supported depth texture format corresponds to the specified image format");
    }
  }

  enum { COLOR_ORDER_RGB, COLOR_ORDER_BGR } color_order;
  color_order = (IImage::Channel::RED < IImage::Channel::BLUE ? COLOR_ORDER_RGB : COLOR_ORDER_BGR);

  switch (type)
  {
    case IImage::Type::LUMINANCE_8U  : return L8();
    case IImage::Type::LUMINANCE_16U : return L16();
    case IImage::Type::LUMINANCE_32F : return L32F();
    case IImage::Type::RGB_8U        : return color_order == COLOR_ORDER_RGB ? RGB8() : BGR8();
    case IImage::Type::RGBA_8U       : return color_order == COLOR_ORDER_RGB ? RGBA8() : BGRA8();
    case IImage::Type::RGB_16U       : return color_order == COLOR_ORDER_RGB ? RGB16() : BGR16();
    case IImage::Type::RGBA_16U      : return color_order == COLOR_ORDER_RGB ? RGBA16() : BGRA16();
    case IImage::Type::RGB_32F       : return color_order == COLOR_ORDER_RGB ? RGB32F() : BGR32F();
    case IImage::Type::RGBA_32F      : return color_order == COLOR_ORDER_RGB ? RGBA32F() : BGRA32F();
    default: throw Error("TextureFormat: No supported texture format corresponds to the specified image format");
  }
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

  return nullptr;
}

// Helper variables for defining texture formats

// Is floating point format
static int8 const FLOAT_FORMAT  = 1;
static int8 const INT_FORMAT    = 0;

// Is opaque format (no alpha)
static int8 const OPAQUE_FORMAT = 1;
static int8 const CLEAR_FORMAT  = 0;

// Is compressed format (not raw component data)
static int8 const COMP_FORMAT   = 1;
static int8 const UNCOMP_FORMAT = 0;

#define THEA_TEXTURE_FORMAT(enumname, cmpnts, cmprssd, opq, fp, code, cs, lb, rb, gb, bb, ab, sb, db, cbpp, glf, glbf, gbpp, gldf) \
  TextureFormat const * TextureFormat::enumname() { \
    static TextureFormat const format(cmpnts, cmprssd, opq, fp, ITextureFormat::Code::code, ITextureFormat::ColorSpace::cs, BayerPattern::NONE, lb, rb, gb, bb, ab, sb, db, cbpp, glf, glbf, gbpp, gldf); \
    return &format; }

THEA_TEXTURE_FORMAT(L8,               1, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   L8,               NONE,    8,  0,  0,  0,  0,  0,  0, 8,      GL_LUMINANCE8,                          GL_LUMINANCE,            8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(L16,              1, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   L16,              NONE,   16,  0,  0,  0,  0,  0,  0, 16,     GL_LUMINANCE16,                         GL_LUMINANCE,           16, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(L16F,             1, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, L16F,             NONE,   16,  0,  0,  0,  0,  0,  0, 16,     GL_LUMINANCE16F_ARB,                    GL_LUMINANCE,           16, GL_FLOAT);
THEA_TEXTURE_FORMAT(L32F,             1, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, L32F,             NONE,   32,  0,  0,  0,  0,  0,  0, 32,     GL_LUMINANCE32F_ARB,                    GL_LUMINANCE,           32, GL_FLOAT);
THEA_TEXTURE_FORMAT(A8,               1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   A8,               NONE,    0,  0,  0,  0,  8,  0,  0, 8,      GL_ALPHA8,                              GL_ALPHA,                8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(A16,              1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   A16,              NONE,    0,  0,  0,  0, 16,  0,  0, 16,     GL_ALPHA16,                             GL_ALPHA,               16, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(A16F,             1, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, A16F,             NONE,    0,  0,  0,  0, 16,  0,  0, 16,     GL_ALPHA16F_ARB,                        GL_ALPHA,               16, GL_FLOAT);
THEA_TEXTURE_FORMAT(A32F,             1, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, A32F,             NONE,    0,  0,  0,  0, 32,  0,  0, 32,     GL_ALPHA32F_ARB,                        GL_ALPHA,               32, GL_FLOAT);
THEA_TEXTURE_FORMAT(LA4,              2, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   LA4,              NONE,    4,  0,  0,  0,  4,  0,  0, 8,      GL_LUMINANCE4_ALPHA4,                   GL_LUMINANCE_ALPHA,      8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(LA8,              2, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   LA8,              NONE,    8,  0,  0,  0,  8,  0,  0, 16,     GL_LUMINANCE8_ALPHA8,                   GL_LUMINANCE_ALPHA,     16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(LA16,             2, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   LA16,             NONE,   16,  0,  0,  0, 16,  0,  0, 16 * 2, GL_LUMINANCE16_ALPHA16,                 GL_LUMINANCE_ALPHA, 16 * 2, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(LA16F,            2, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, LA16F,            NONE,   16,  0,  0,  0, 16,  0,  0, 16 * 2, GL_LUMINANCE_ALPHA16F_ARB,              GL_LUMINANCE_ALPHA, 16 * 2, GL_FLOAT);
THEA_TEXTURE_FORMAT(LA32F,            2, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, LA32F,            NONE,   32,  0,  0,  0, 32,  0,  0, 32 * 2, GL_LUMINANCE_ALPHA32F_ARB,              GL_LUMINANCE_ALPHA, 32 * 2, GL_FLOAT);
THEA_TEXTURE_FORMAT(R8,               1, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   R8,               RGB,     0,  8,  0,  0,  0,  0,  0, 8,      GL_R8,                                  GL_RED,                  8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RG8,              2, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RG8,              RGB,     0,  8,  8,  0,  0,  0,  0, 16,     GL_RG8,                                 GL_RG,                  16, GL_UNSIGNED_BYTE);
// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
THEA_TEXTURE_FORMAT(RG8I,             2, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RG8I,             RGB,     0,  8,  8,  0,  0,  0,  0, 16,     GL_RG8I,                                GL_RG_INTEGER,          16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RG8UI,            2, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RG8UI,            RGB,     0,  8,  8,  0,  0,  0,  0, 16,     GL_RG8UI,                               GL_RG_INTEGER,          16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RG16F,            2, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, RG16F,            RGB,     0, 16, 16,  0,  0,  0,  0, 32,     GL_RG16F,                               GL_RG,                  32, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGB5,             3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB5,             RGB,     0,  5,  5,  5,  0,  0,  0, 16,     GL_RGB5,                                GL_RGBA,                16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGB5A1,           4, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB5A1,           RGB,     0,  5,  5,  5,  1,  0,  0, 16,     GL_RGB5_A1,                             GL_RGBA,                16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGB8,             3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB8,             RGB,     0,  8,  8,  8,  0,  0,  0, 24,     GL_RGB8,                                GL_RGB,                 32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGB10,            3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB10,            RGB,     0, 10, 10, 10,  0,  0,  0, 10 * 3, GL_RGB10,                               GL_RGB,                 32, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(RGB10A2,          4, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB10A2,          RGB,     0, 10, 10, 10,  2,  0,  0, 32,     GL_RGB10_A2,                            GL_RGBA,                32, GL_UNSIGNED_INT_10_10_10_2);
THEA_TEXTURE_FORMAT(RGB16,            3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB16,            RGB,     0, 16, 16, 16,  0,  0,  0, 16 * 3, GL_RGB16,                               GL_RGB,             16 * 3, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(RGB16F,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB16F,           RGB,     0, 16, 16, 16,  0,  0,  0, 16 * 3, GL_RGB16F_ARB,                          GL_RGB,             16 * 3, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGB32F,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB32F,           RGB,     0, 32, 32, 32,  0,  0,  0, 32 * 3, GL_RGB32F_ARB,                          GL_RGB,             32 * 3, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGBA8,            4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   RGBA8,            RGB,     0,  8,  8,  8,  8,  0,  0, 32,     GL_RGBA8,                               GL_RGBA,                32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGBA16,           4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   RGBA16,           RGB,     0, 16, 16, 16, 16,  0,  0, 16 * 4, GL_RGBA16,                              GL_RGBA,            16 * 4, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(RGBA16F,          4, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, RGBA16F,          RGB,     0, 16, 16, 16, 16,  0,  0, 16 * 4, GL_RGBA16F_ARB,                         GL_RGBA,            16 * 4, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGBA32F,          4, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, RGBA32F,          RGB,     0, 32, 32, 32, 32,  0,  0, 32 * 4, GL_RGBA32F_ARB,                         GL_RGBA,            32 * 4, GL_FLOAT);
// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
THEA_TEXTURE_FORMAT(RGBA32UI,         4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   RGBA32UI,         RGB,     0, 32, 32, 32, 32,  0,  0, 32 * 4, GL_RGBA32UI,                            GL_RGBA_INTEGER,    32 * 4, GL_UNSIGNED_INT);
THEA_TEXTURE_FORMAT(R11G11B10F,       3, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, R11G11B10F,       RGB,     0, 11, 11, 10,  0,  0,  0, 32,     GL_R11F_G11F_B10F_EXT,                  GL_RGB,                 32, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGB9E5F,          3, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, RGB9E5F,          RGB,     0, 14, 14, 14,  0,  0,  0, 32,     GL_RGB9_E5_EXT,                         GL_RGB,                 32, GL_FLOAT);
// The base format for integer formats must be *_INTEGER even though the spec doesn't state this
// FIXME: Should openGlBitsPerPixel() really be 32 and not 24 for RGB8, BGR8 etc?
THEA_TEXTURE_FORMAT(RGB8I,            3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB8I,            RGB,     0,  8,  8,  8,  0,  0,  0, 24,     GL_RGB8I_EXT,                           GL_RGB_INTEGER,         32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGB8UI,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGB8UI,           RGB,     0,  8,  8,  8,  0,  0,  0, 24,     GL_RGB8UI_EXT,                          GL_RGB_INTEGER,         32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGBA8UI,          4, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   RGBA8UI,          RGB,     0,  8,  8,  8,  0,  0,  8, 32,     GL_RGBA8UI_EXT,                         GL_RGBA_INTEGER,        32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(BGR8,             3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   BGR8,             RGB,     0,  8,  8,  8,  0,  0,  0, 24,     GL_RGB8,                                GL_BGR,                 32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(BGRA8,            4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   BGRA8,            RGB,     0,  8,  8,  8,  8,  0,  0, 32,     GL_RGBA8,                               GL_BGRA,                32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(BGR16,            3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   BGR16,            RGB,     0, 16, 16, 16,  0,  0,  0, 16 * 3, GL_RGB16,                               GL_BGR,             16 * 3, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(BGRA16,           4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   BGRA16,           RGB,     0, 16, 16, 16, 16,  0,  0, 16 * 4, GL_RGBA16,                              GL_BGRA,            16 * 4, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(BGR32F,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, FLOAT_FORMAT, BGR32F,           RGB,     0, 32, 32, 32,  0,  0,  0, 32 * 3, GL_RGB32F_ARB,                          GL_BGR,             32 * 3, GL_FLOAT);
THEA_TEXTURE_FORMAT(BGRA32F,          4, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, BGRA32F,          RGB,     0, 32, 32, 32, 32,  0,  0, 32 * 4, GL_RGBA32F_ARB,                         GL_BGRA,            32 * 4, GL_FLOAT);
THEA_TEXTURE_FORMAT(RGB_DXT1,         3, COMP_FORMAT,   OPAQUE_FORMAT, INT_FORMAT,   RGB_DXT1,         RGB,     0,  0,  0,  0,  0,  0,  0, 64,     GL_COMPRESSED_RGB_S3TC_DXT1_EXT,        GL_RGB,                 64, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGBA_DXT1,        4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   RGBA_DXT1,        RGB,     0,  0,  0,  0,  0,  0,  0, 64,     GL_COMPRESSED_RGBA_S3TC_DXT1_EXT,       GL_RGBA,                64, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGBA_DXT3,        4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   RGBA_DXT3,        RGB,     0,  0,  0,  0,  0,  0,  0, 128,    GL_COMPRESSED_RGBA_S3TC_DXT3_EXT,       GL_RGBA,               128, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(RGBA_DXT5,        4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   RGBA_DXT5,        RGB,     0,  0,  0,  0,  0,  0,  0, 128,    GL_COMPRESSED_RGBA_S3TC_DXT5_EXT,       GL_RGBA,               128, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGB8,            3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   SRGB8,            SRGB,    0,  8,  8,  8,  0,  0,  0, 24,     GL_SRGB8,                               GL_RGB,                 32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGBA8,           4, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   SRGBA8,           SRGB,    0,  8,  8,  8,  8,  0,  0, 24,     GL_SRGB8_ALPHA8,                        GL_RGBA,                32, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SL8,              1, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   SL8,              SRGB,    8,  0,  0,  0,  0,  0,  0, 8,      GL_SLUMINANCE8,                         GL_LUMINANCE,            8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SLA8,             2, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   SLA8,             SRGB,    8,  0,  0,  0,  8,  0,  0, 16,     GL_SLUMINANCE8_ALPHA8,                  GL_LUMINANCE_ALPHA,     16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGB_DXT1,        3, COMP_FORMAT,   OPAQUE_FORMAT, INT_FORMAT,   SRGB_DXT1,        SRGB,    0,  0,  0,  0,  0,  0,  0, 64,     GL_COMPRESSED_SRGB_S3TC_DXT1_EXT,       GL_RGB,                 64, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGBA_DXT1,       4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   SRGBA_DXT1,       SRGB,    0,  0,  0,  0,  0,  0,  0, 64,     GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT, GL_RGBA,                64, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGBA_DXT3,       4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   SRGBA_DXT3,       SRGB,    0,  0,  0,  0,  0,  0,  0, 128,    GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT, GL_RGBA,               128, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(SRGBA_DXT5,       4, COMP_FORMAT,   CLEAR_FORMAT,  INT_FORMAT,   SRGBA_DXT5,       SRGB,    0,  0,  0,  0,  0,  0,  0, 128,    GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT, GL_RGBA,               128, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(DEPTH16,          1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   DEPTH16,          NONE,    0,  0,  0,  0,  0, 0, 16, 16,      GL_DEPTH_COMPONENT16_ARB,               GL_DEPTH_COMPONENT,     16, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(DEPTH24,          1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   DEPTH24,          NONE,    0,  0,  0,  0,  0, 0, 24, 24,      GL_DEPTH_COMPONENT24_ARB,               GL_DEPTH_COMPONENT,     32, GL_UNSIGNED_INT);
THEA_TEXTURE_FORMAT(DEPTH32,          1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   DEPTH32,          NONE,    0,  0,  0,  0,  0, 0, 32, 32,      GL_DEPTH_COMPONENT32_ARB,               GL_DEPTH_COMPONENT,     32, GL_UNSIGNED_INT);
THEA_TEXTURE_FORMAT(DEPTH32F,         1, UNCOMP_FORMAT, CLEAR_FORMAT,  FLOAT_FORMAT, DEPTH32F,         NONE,    0,  0,  0,  0,  0, 0, 32, 32,      GL_DEPTH_COMPONENT32_ARB,               GL_DEPTH_COMPONENT,     32, GL_FLOAT);
// These formats are for use with Renderbuffers only!
THEA_TEXTURE_FORMAT(STENCIL1,         1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   STENCIL1,         NONE,    0,  0,  0,  0,  0,  1,  0,  1,     GL_STENCIL_INDEX1_EXT,                  GL_STENCIL_INDEX,        1, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(STENCIL4,         1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   STENCIL4,         NONE,    0,  0,  0,  0,  0,  4,  0,  4,     GL_STENCIL_INDEX4_EXT,                  GL_STENCIL_INDEX,        4, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(STENCIL8,         1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   STENCIL8,         NONE,    0,  0,  0,  0,  0,  8,  0,  8,     GL_STENCIL_INDEX8_EXT,                  GL_STENCIL_INDEX,        8, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(STENCIL16,        1, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   STENCIL16,        NONE,    0,  0,  0,  0,  0, 16,  0, 16,     GL_STENCIL_INDEX16_EXT,                 GL_STENCIL_INDEX,       16, GL_UNSIGNED_SHORT);
THEA_TEXTURE_FORMAT(DEPTH24_STENCIL8, 2, UNCOMP_FORMAT, CLEAR_FORMAT,  INT_FORMAT,   DEPTH24_STENCIL8, NONE,    0,  0,  0,  0,  0,  8, 24, 32,     GL_DEPTH24_STENCIL8_EXT,                GL_DEPTH_STENCIL_EXT,   32, GL_UNSIGNED_INT);
THEA_TEXTURE_FORMAT(YUV420_PLANAR,    3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   YUV420_PLANAR,    YUV,     0,  0,  0,  0,  0,  0,  0, 12,     GL_NONE,                                GL_NONE,                12, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(YUV422,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   YUV422,           YUV,     0,  0,  0,  0,  0,  0,  0, 16,     GL_NONE,                                GL_NONE,                16, GL_UNSIGNED_BYTE);
THEA_TEXTURE_FORMAT(YUV444,           3, UNCOMP_FORMAT, OPAQUE_FORMAT, INT_FORMAT,   YUV444,           YUV,     0,  0,  0,  0,  0,  0,  0, 24,     GL_NONE,                                GL_NONE,                24, GL_UNSIGNED_BYTE);

} // namespace Graphics
} // namespace Thea
