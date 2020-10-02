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

  @file ImageFormat.h

  @maintainer Morgan McGuire, http://graphics.cs.williams.edu

  @created 2003-05-23
  @edited  2010-05-01
*/

#ifndef __Thea_Graphics_TextureFormat_hpp__
#define __Thea_Graphics_TextureFormat_hpp__

#include "../Common.hpp"
#include "../IImage.hpp"
#include "../NamedObject.hpp"

namespace Thea {
namespace Graphics {

/**
 * Interface for accessing common GPU texture formats.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API ITextureFormat : public virtual INamedObject
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

    /** Destructor. */
    virtual ~ITextureFormat() = 0;

    /**  Number of channels (1 for a depth texture). */
    virtual int32 THEA_ICALL numComponents() const = 0;

    /** Is the format compressed? */
    virtual int8 THEA_ICALL isCompressed() const = 0;

    /** Is this a depth texture format? */
    virtual int8 THEA_ICALL isDepth() const = 0;

    /** True if there is no alpha channel for this texture. */
    virtual int8 THEA_ICALL isOpaque() const = 0;

    /** True if the format has floating-point channels. */
    virtual int8 THEA_ICALL isFloatingPoint() const = 0;

    /** Code identifying the format, as a value from the Code enum. */
    virtual int32 THEA_ICALL code() const = 0;

    /** Color space of the format, as a value from the ColorSpace enum. */
    virtual int32 THEA_ICALL colorSpace() const = 0;

    /** If this is a Bayer format, what is the pattern? Returns a value from the BayerPattern enum. */
    virtual int32 THEA_ICALL bayerPattern() const = 0;

    /** Number of bits assigned to luminance. */
    virtual int32 THEA_ICALL luminanceBits() const = 0;

    /** Number of bits per pixel for the red channel. Zero for compressed textures and non-RGB. */
    virtual int32 THEA_ICALL redBits() const = 0;

    /** Number of bits per pixel for the green channel. Zero for compressed textures and non-RGB. */
    virtual int32 THEA_ICALL greenBits() const = 0;

    /** Number of bits per pixel for the blue channel. Zero for compressed textures and non-RGB. */
    virtual int32 THEA_ICALL blueBits() const = 0;

    /** Number of bits per pixel for the alpha channel. Zero for compressed textures and non-RGB. */
    virtual int32 THEA_ICALL alphaBits() const = 0;

    /** Number of bits per pixel for stencil buffers */
    virtual int32 THEA_ICALL stencilBits() const = 0;

    /** Number of depth bits (for depth textures, such as shadow maps) */
    virtual int32 THEA_ICALL depthBits() const = 0;

    /** Amount of CPU memory per pixel when packed into an array, discounting any end-of-row padding. */
    virtual int32 THEA_ICALL cpuBitsPerPixel() const = 0;

    /** The OpenGL format equivalent to this one, such as GL_RGB8. Zero if there is no equivalent. This is actually a GLenum. */
    virtual int32 THEA_ICALL openGlFormat() const = 0;

    /** The OpenGL base format equivalent to this one (such as GL_RGB, GL_ALPHA). Zero if there is no equivalent.  */
    virtual int32 THEA_ICALL openGlBaseFormat() const = 0;

    /**
     * Amount of GPU memory per pixel on most graphics cards, for formats supported by OpenGL. This is only an estimate -- the
     * actual amount of memory may be different on your actual card.
     *
     * This may be greater than the sum of the per-channel bits because graphics cards need to pad to the nearest 1, 2, or 4
     * bytes.
     */
    virtual int32 THEA_ICALL openGlBitsPerPixel() const = 0;

    /** The OpenGL bytes (type) format of the data buffer used with this texture format, such as GL_UNSIGNED_BYTE. */
    virtual int32 THEA_ICALL openGlDataFormat() const = 0;

}; // class ITextureFormat

inline ITextureFormat::~ITextureFormat() {}

/**
 * Encapsulates information about common GPU texture formats. Don't construct these; use the static methods provided. For most
 * formats, the number indicates the number of bits per channel and a suffix of "F" indicates floating point. This does not hold
 * for the YUV and DXT formats.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API TextureFormat : public virtual ITextureFormat
{
  private:
    TextureFormat(
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
    );

    int32  num_components;
    int8   is_compressed;
    int8   is_opaque;
    int8   is_floating_point;
    int32  code_val;
    int32  color_space;
    int32  bayer_pattern;
    int32  luminance_bits;
    int32  red_bits;
    int32  green_bits;
    int32  blue_bits;
    int32  alpha_bits;
    int32  stencil_bits;
    int32  depth_bits;
    int32  cpu_bits_per_pixel;
    int32  gl_format;
    int32  gl_base_format;
    int32  gl_bits_per_pixel;
    int32  gl_data_format;

  public:
    int32 THEA_ICALL numComponents() const       { return num_components;     }
    int8 THEA_ICALL isCompressed() const         { return is_compressed;      }
    int8 THEA_ICALL isDepth() const              { return depth_bits > 0;     }
    int8 THEA_ICALL isOpaque() const             { return is_opaque;          }
    int8 THEA_ICALL isFloatingPoint() const      { return is_floating_point;  }
    int32 THEA_ICALL code() const                { return code_val;           }
    int32 THEA_ICALL colorSpace() const          { return color_space;        }
    int32 THEA_ICALL bayerPattern() const        { return bayer_pattern;      }
    int32 THEA_ICALL luminanceBits() const       { return luminance_bits;     }
    int32 THEA_ICALL redBits() const             { return red_bits;           }
    int32 THEA_ICALL greenBits() const           { return green_bits;         }
    int32 THEA_ICALL blueBits() const            { return blue_bits;          }
    int32 THEA_ICALL alphaBits() const           { return alpha_bits;         }
    int32 THEA_ICALL stencilBits() const         { return stencil_bits;       }
    int32 THEA_ICALL depthBits() const           { return depth_bits;         }
    int32 THEA_ICALL cpuBitsPerPixel() const     { return cpu_bits_per_pixel; }
    int32 THEA_ICALL openGlFormat() const        { return gl_format;          }
    int32 THEA_ICALL openGlBaseFormat() const    { return gl_base_format;     }
    int32 THEA_ICALL openGlBitsPerPixel() const  { return gl_bits_per_pixel;  }
    int32 THEA_ICALL openGlDataFormat() const    { return gl_data_format;     }

    char const * THEA_ICALL getName() const;

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
    static TextureFormat const * AUTO() { return nullptr; }

    /**
     * Get a depth image format. Returns DEPTH16, DEPTH24, or DEPTH32 according to the bits specified. You can use
     * "glGetInteger(GL_DEPTH_BITS)" to match the screen's format.
     */
    static TextureFormat const * depth(int depth_bits = 24);

    /** Returns STENCIL1, STENCIL4, STENCIL8 or STENCIL16 according to the bits
      specified. You can use "glGetInteger(GL_STENCIL_BITS)" to match
     the screen's format.*/
    static TextureFormat const * stencil(int bits = 8);

    /** Returns the matching TextureFormat* identified by the Code.  May return nullptr
      if this format's code is reserved but not yet implemented by Thea. */
    static TextureFormat const * fromCode(TextureFormat::Code code);

    /** Internal convenience function to convert an image type to a texture storage format. */
    static TextureFormat const * fromImageType(IImage::Type type, bool is_depth = false);

    /** Get an image format given its name. Takes the same values that name() returns. */
    static TextureFormat const * fromString(std::string const & s);

}; // class TextureFormat

} // namespace Graphics
} // namespace Thea

#endif
