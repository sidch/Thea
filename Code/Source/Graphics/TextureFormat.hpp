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

  @file ImageFormat.h

  @maintainer Morgan McGuire, http://graphics.cs.williams.edu

  @created 2003-05-23
  @edited  2010-05-01
*/

#ifndef __Thea_Graphics_TextureFormat_hpp__
#define __Thea_Graphics_TextureFormat_hpp__

#include "../Common.hpp"
#include "../Image.hpp"

namespace Thea {
namespace Graphics {

/**
 * Information about common GPU texture formats. Don't construct these; use the methods provided. For most formats, the number
 * indicates the number of bits per channel and a suffix of "F" indicates floating point. This does not hold for the YUV and DXT
 * formats.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API TextureFormat
{
  public:

    /** Code identifying the format (enum class). */
    struct THEA_API Code
    {
      /** Supported values. */
      enum Value
      {
        NONE = -1,
        L8,
        L16,
        L16F,
        L32F,

        A8,
        A16,
        A16F,
        A32F,

        LA4,
        LA8,
        LA16,
        LA16F,
        LA32F,

        RGB5,
        RGB5A1,
        RGB8,
        RGB10,
        RGB10A2,
        RGB16,
        RGB16F,
        RGB32F,
        R11G11B10F,
        RGB9E5F,

        RGB8I,
        RGB8UI,

        RGBA8UI,

        R8,

        RG8,
        RG8I,
        RG8UI,

        RG16F,

        RGBA8,
        RGBA16,
        RGBA16F,
        RGBA32F,

        RGBA32UI,

        BGR8,
        BGRA8,
        BGR16,
        BGRA16,
        BGR32F,
        BGRA32F,
        ARGB8,

        BAYER_RGGB8,
        BAYER_GRBG8,
        BAYER_GBRG8,
        BAYER_BGGR8,
        BAYER_RGGB32F,
        BAYER_GRBG32F,
        BAYER_GBRG32F,
        BAYER_BGGR32F,

        HSV8,
        HSV32F,

        YUV420_PLANAR,
        YUV422,
        YUV444,

        RGB_DXT1,
        RGBA_DXT1,
        RGBA_DXT3,
        RGBA_DXT5,

        SRGB8,
        SRGBA8,

        SL8,
        SLA8,

        SRGB_DXT1,
        SRGBA_DXT1,
        SRGBA_DXT3,
        SRGBA_DXT5,

        DEPTH16,
        DEPTH24,
        DEPTH32,
        DEPTH32F,

        STENCIL1,
        STENCIL4,
        STENCIL8,
        STENCIL16,

        DEPTH24_STENCIL8,

        NUM
      };

      THEA_ENUM_CLASS_BODY(Code)
    };

    /** Color space of the format (enum class). */
    struct THEA_API ColorSpace
    {
      /** Supported values. */
      enum Value
      {
        NONE,
        RGB,
        HSV,
        YUV,
        SRGB
      };

      THEA_ENUM_CLASS_BODY(ColorSpace)
    };

    /** Bayer pattern of the format (enum class). */
    struct THEA_API BayerPattern
    {
      /** Supported values. */
      enum Value
      {
        NONE,
        RGGB,
        GRBG,
        GBRG,
        BGGR
      };

      THEA_ENUM_CLASS_BODY(BayerPattern)
    };

    /**  Number of channels (1 for a depth texture). */
    int                 numComponents;

    /** Is the format compressed? */
    bool                compressed;

    /** Code identifying the format. */
    Code                code;

    /** Color space of the format. */
    ColorSpace          colorSpace;

    /** If this is a Bayer format, what is the pattern? */
    BayerPattern        bayerPattern;

    /** The OpenGL format equivalent to this one, such as GL_RGB8. Zero if there is no equivalent. This is actually a GLenum. */
    int                 openGLFormat;

    /** The OpenGL base format equivalent to this one (such as GL_RGB, GL_ALPHA). Zero if there is no equivalent.  */
    int                 openGLBaseFormat;

    /** Number of bits assigned to luminance. */
    int                 luminanceBits;

    /** Number of bits per pixel storage for alpha values; Zero for compressed textures and non-RGB. */
    int                 alphaBits;

    /** Number of bits per pixel storage for red values; Zero for compressed textures and non-RGB. */
    int                 redBits;

    /** Number of bits per pixel storage for green values; Zero for compressed textures and non-RGB. */
    int                 greenBits;

    /** Number of bits per pixel storage for blue values; Zero for compressed textures  and non-RGB. */
    int                 blueBits;

    /** Number of bits per pixel */
    int                 stencilBits;

    /** Number of depth bits (for depth textures, such as shadow maps) */
    int                 depthBits;

    /** Amount of CPU memory per pixel when packed into an array, discounting any end-of-row padding. */
    int                 cpuBitsPerPixel;

    /**
     * Amount of GPU memory per pixel on most graphics cards, for formats supported by OpenGL. This is only an estimate -- the
     * actual amount of memory may be different on your actual card.
     *
     * This may be greater than the sum of the per-channel bits because graphics cards need to pad to the nearest 1, 2, or 4
     * bytes.
     */
    int                 openGLBitsPerPixel;

    /** The OpenGL bytes (type) format of the data buffer used with this texture format, such as GL_UNSIGNED_BYTE. */
    int                 openGLDataFormat;

    /** True if there is no alpha channel for this texture. */
    bool                opaque;

    /** True if the format has floating-point channels. */
    bool                floatingPoint;

    /** Check if this is a depth texture format. */
    bool isDepth() const;

    /** Human readable name of this format.*/
    std::string const & name() const;

    /** Get an image format given its name. Takes the same values that name() returns. */
    static TextureFormat const * fromString(std::string const & s);

  private:
    TextureFormat(
      int             numComponents,
      bool            compressed,
      int             glFormat,
      int             glBaseFormat,
      int             luminanceBits,
      int             alphaBits,
      int             redBits,
      int             greenBits,
      int             blueBits,
      int             depthBits,
      int             stencilBits,
      int             openGLBitsPerPixel,
      int             cpuBitsPerPixel,
      int             glDataFormat,
      bool            opaque,
      bool            floatingPoint,
      Code            code,
      ColorSpace      colorSpace,
      BayerPattern    bayerPattern = BayerPattern::NONE);

  public:
    /** L8 format. */
    static TextureFormat const * L8();

    /** L16 format. */
    static TextureFormat const * L16();

    /** L16F format. */
    static TextureFormat const * L16F();

    /** L32F format. */
    static TextureFormat const * L32F();

    /** A8 format. */
    static TextureFormat const * A8();

    /** A16 format. */
    static TextureFormat const * A16();

    /** A16F format. */
    static TextureFormat const * A16F();

    /** A32F format. */
    static TextureFormat const * A32F();

    /** LA4 format. */
    static TextureFormat const * LA4();

    /** LA8 format. */
    static TextureFormat const * LA8();

    /** LA16 format. */
    static TextureFormat const * LA16();

    /** LA16F format. */
    static TextureFormat const * LA16F();

    /** LA32F format. */
    static TextureFormat const * LA32F();

    /** R8 format. */
    static TextureFormat const * R8();

    /** RG8 format. */
    static TextureFormat const * RG8();

    /** RG8I format. */
    static TextureFormat const * RG8I();

    /** RG8UI format. */
    static TextureFormat const * RG8UI();

    /** RG16F format. */
    static TextureFormat const * RG16F();

    /** RGB5 format. */
    static TextureFormat const * RGB5();

    /** RGB5A1 format. */
    static TextureFormat const * RGB5A1();

    /** RGB8 format. */
    static TextureFormat const * RGB8();

    /** RGB10 format. */
    static TextureFormat const * RGB10();

    /** RGB10A2 format. */
    static TextureFormat const * RGB10A2();

    /** RGB16 format. */
    static TextureFormat const * RGB16();

    /** RGB16F format. */
    static TextureFormat const * RGB16F();

    /** RGB32F format. */
    static TextureFormat const * RGB32F();

    /** RGBA8 format. */
    static TextureFormat const * RGBA8();

    /** RGBA16 format. */
    static TextureFormat const * RGBA16();

    /** RGBA16F format. */
    static TextureFormat const * RGBA16F();

    /** RGBA32F format. */
    static TextureFormat const * RGBA32F();

    /** RGBA32UI format. */
    static TextureFormat const * RGBA32UI();

    /** R11G11B10F format. */
    static TextureFormat const * R11G11B10F();

    /** RGB9E5F format. */
    static TextureFormat const * RGB9E5F();

    /** RGB8I format. */
    static TextureFormat const * RGB8I();

    /** RGB8UI format. */
    static TextureFormat const * RGB8UI();

    /** RGBA8UI format. */
    static TextureFormat const * RGBA8UI();

    /** BGR8 format. */
    static TextureFormat const * BGR8();

    /** BGRA8 format. */
    static TextureFormat const * BGRA8();

    /** BGR16 format. */
    static TextureFormat const * BGR16();

    /** BGRA16 format. */
    static TextureFormat const * BGRA16();

    /** BGR32F format. */
    static TextureFormat const * BGR32F();

    /** BGRA32F format. */
    static TextureFormat const * BGRA32F();

    /** RGB_DXT1 format. */
    static TextureFormat const * RGB_DXT1();

    /** RGBA_DXT1 format. */
    static TextureFormat const * RGBA_DXT1();

    /** RGBA_DXT3 format. */
    static TextureFormat const * RGBA_DXT3();

    /** RGBA_DXT5 format. */
    static TextureFormat const * RGBA_DXT5();

    /** SRGB8 format. */
    static TextureFormat const * SRGB8();

    /** SRGBA8 format. */
    static TextureFormat const * SRGBA8();

    /** SL8 format. */
    static TextureFormat const * SL8();

    /** SLA8 format. */
    static TextureFormat const * SLA8();

    /** SRGB_DXT1 format. */
    static TextureFormat const * SRGB_DXT1();

    /** SRGBA_DXT1 format. */
    static TextureFormat const * SRGBA_DXT1();

    /** SRGBA_DXT3 format. */
    static TextureFormat const * SRGBA_DXT3();

    /** SRGBA_DXT5 format. */
    static TextureFormat const * SRGBA_DXT5();

    /** DEPTH16 format. */
    static TextureFormat const * DEPTH16();

    /** DEPTH24 format. */
    static TextureFormat const * DEPTH24();

    /** DEPTH32 format. */
    static TextureFormat const * DEPTH32();

    /** DEPTH32F format. */
    static TextureFormat const * DEPTH32F();

    /** STENCIL1 format. */
    static TextureFormat const * STENCIL1();

    /** STENCIL4 format. */
    static TextureFormat const * STENCIL4();

    /** STENCIL8 format. */
    static TextureFormat const * STENCIL8();

    /** STENCIL16 format. */
    static TextureFormat const * STENCIL16();

    /** DEPTH24_STENCIL8 format. */
    static TextureFormat const * DEPTH24_STENCIL8();

    /** YUV420_PLANAR format. */
    static TextureFormat const * YUV420_PLANAR();

    /** YUV422 format. */
    static TextureFormat const * YUV422();

    /** YUV444 format. */
    static TextureFormat const * YUV444();

    /**
     * Indicates the format should be automatically detected. The usual choice is either RGBA8 or RGB8 depending on the presence
     * of an alpha channel in the input.
     */
    static TextureFormat const * AUTO()
    {
      return NULL;
    }

    /**
     * Get a depth image format. Returns DEPTH16, DEPTH24, or DEPTH32 according to the bits specified. You can use
     * "glGetInteger(GL_DEPTH_BITS)" to match the screen's format.
     */
    static TextureFormat const * depth(int depth_bits = 24);

    /** Returns STENCIL1, STENCIL4, STENCIL8 or STENCIL16 according to the bits
      specified. You can use "glGetInteger(GL_STENCIL_BITS)" to match
     the screen's format.*/
    static TextureFormat const * stencil(int bits = 8);

    /** Returns the matching TextureFormat* identified by the Code.  May return NULL
      if this format's code is reserved but not yet implemented by Thea. */
    static TextureFormat const * fromCode(TextureFormat::Code code);

    /** Internal convenience function to convert an image type to a texture storage format. */
    static TextureFormat const * fromImageType(AbstractImage::Type type, bool is_depth = false);

}; // class TextureFormat

} // namespace Graphics
} // namespace Thea

#endif
