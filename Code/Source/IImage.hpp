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
// First version: 2009
//
//============================================================================

#ifndef __Thea_IImage_hpp__
#define __Thea_IImage_hpp__

#include "Common.hpp"

namespace Thea {

/**
 * Interface for 2D and 3D images, useful for passing images across shared library boundaries. The image is assumed to be
 * depth-major, i.e. the depth dimension is the most significant, followed by height, followed by width. A 3D image is therefore
 * a set of 2D images stacked along the depth dimension.
 */
class THEA_API IImage
{
  public:
    THEA_DECL_SMART_POINTERS(IImage)

    /**
     * Indices for the color channels in a single pixel. E.g. for a 24-bit RGB bitmap (8 bits per channel), the individual
     * components for the pixel at address <code>unsigned char * p</code> may be accessed as <code>p[Channel::RED]</code>,
     * <code>p[Channel::GREEN]</code> and <code>p[Channel::BLUE]</code>.
     */
    struct THEA_API Channel
    {
      static int32 const RED;    ///< Index of red channel
      static int32 const GREEN;  ///< Index of green channel
      static int32 const BLUE;   ///< Index of blue channel
      static int32 const ALPHA;  ///< Index of alpha channel
    };

    /**
     * Different image types (enum class plus extra functions). Note that the channel ordering can be OS dependent, hence access
     * is best achieved using the implementation-specific RED, GREEN, BLUE and ALPHA indices.
     */
    struct THEA_API Type
    {
      /** Supported values. */
      enum Value
      {
        LUMINANCE_1U,   ///< 1-bit luminance image.
        LUMINANCE_2U,   ///< 2-bit luminance image.
        LUMINANCE_4U,   ///< 4-bit luminance image.
        LUMINANCE_8U,   ///< 8-bit luminance image.
        LUMINANCE_16,   ///< 16-bit signed luminance image.
        LUMINANCE_16U,  ///< 16-bit unsigned luminance image.
        LUMINANCE_32,   ///< 32-bit signed luminance image.
        LUMINANCE_32U,  ///< 32-bit unsigned luminance image.
        LUMINANCE_32F,  ///< 32-bit floating-point luminance image.
        LUMINANCE_64F,  ///< 64-bit floating-point luminance image.
        RGB_8U,         ///< RGB image, 8 bits per channel.
        RGBA_8U,        ///< RGBA image, 8 bits per channel.
        RGB_16U,        ///< RGB image, 16 bits per channel.
        RGBA_16U,       ///< RGBA image, 16 bits per channel.
        RGB_32F,        ///< RGB image, 32-bit floating point per channel.
        RGBA_32F,       ///< RGBA image, 32-bit floating point per channel.
        COMPLEX_64F,    ///< Image of complex numbers with 64-bit floating point real and imaginary parts.
        UNKNOWN,        ///< Unknown format.
      };

      THEA_ENUM_CLASS_BODY(Type)

      /**
       * Get the number of channels per pixel. For example: 1 for luminance, 3 for RGB, 4 for RGBA. A complex number corresponds
       * to a single channel (not two). Returns -1 if the image is of unknown type.
       */
      int numChannels() const;

      /** Check if the channels hold complex (as opposed to real) values. */
      bool isComplex() const;

      /** Check if the channels hold floating-point values (at any precision). */
      bool isFloatingPoint() const;

      /** Get the number of bits assigned to each pixel. */
      int getBitsPerPixel() const;

      /** Get the number of bits assigned to each channel. Returns -1 if the channels don't all have the same number of bits. */
      int getBitsPerChannel() const;

      /**
       * Get the number of bits assigned to a particular channel. For luminance images, the single channel is assumed to
       * correspond to the enum value Channel::ALPHA. If the image doesn't contain the specific channel (e.g. luminance images
       * don't have red, green or blue channels) a value of zero is returned.
       */
      int getBitsInChannel(int channel) const;

      /** Check if the image pixels start at byte addresses. For example, RGB_8U does, but LUMINANCE_1U does not. */
      bool hasByteAlignedPixels() const;

      /** Check if all pixel channels are aligned to byte addresses. */
      bool hasByteAlignedChannels() const;

    }; // struct Type

    /** Image resampling filters (enum class). */
    struct Filter
    {
      /** Supported values (copied from FreeImage). */
      enum Value
      {
        BOX,          ///< Box, pulse, Fourier window, 1st order (constant) B-Spline.
        BILINEAR,     ///< Bilinear filter.
        BSPLINE,      ///< 4th order (cubic) B-Spline.
        BICUBIC,      ///< Mitchell and Netravali's two-parameter cubic filter.
        CATMULL_ROM,  ///< Catmull-Rom spline, Overhauser spline.
        LANCZOS3,     ///< Lanczos-windowed sinc filter.
        AUTO,         ///< Automatically choose an appropriate filter.
      };

      THEA_ENUM_CLASS_BODY(Filter);

    }; // struct Filter

    /** Destructor. */
    virtual ~IImage() {}

    /**
     * Check if the image has been allocated non-zero memory space (hence has valid type and dimensions) or not. An image
     * created by the default constructor is invalid and must be further initialized using Image::read() or a similar function.
     */
    virtual int8 THEA_ICALL isValid() const = 0;

    /**
     * Destroy all image data, resetting the image to an invalid state.
     *
     * @return True on success, false on error.
     *
     * @see isValid()
     */
    virtual int8 THEA_ICALL clear() = 0;

    /**
     * Resize the image, changing its type and dimensions. All previous image data is discarded. \a type should be one of the
     * Type values.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL resize(int64 type, int64 width, int64 height, int64 depth = 1) = 0;

    /** Get the width of the image in pixels. */
    virtual int64 THEA_ICALL getWidth() const = 0;

    /** Get the height of the image in pixels. */
    virtual int64 THEA_ICALL getHeight() const = 0;

    /** Get the depth of the image in pixels. */
    virtual int64 THEA_ICALL getDepth() const { return 1; }

    /** Get the type of the image pixels, corresponding to one of the values of the Type enum class. */
    virtual int32 THEA_ICALL getType() const = 0;

    /** Get a pointer to the image data. */
    virtual void const * THEA_ICALL getData() const = 0;

    /** Get a pointer to the image data. */
    virtual void * THEA_ICALL getData() = 0;

    /** Get a pointer to the beginning of a specified row of pixels, optionally for a specific depth slice. */
    virtual void const * THEA_ICALL getScanLine(int64 row, int64 z = 0) const = 0;

    /** Get a pointer to the beginning of a specified row of pixels, optionally for a specific depth slice. */
    virtual void * THEA_ICALL getScanLine(int64 row, int64 z = 0) = 0;

    /**
     * Get the number of bytes consumed by a row of pixels. Rows may be aligned to 32-bit (or other) boundaries for performance
     * reasons, so this is <b>not</b> necessarily equal to the number of pixels in a row times the size of a pixel.
     */
    virtual int64 THEA_ICALL getScanWidth() const = 0;

    /**
     * Get the byte alignment of a pixel row. Each row is padded to take up a number of bytes which is a multiple of the
     * alignment value.
     */
    virtual int32 THEA_ICALL getRowAlignment() const = 0;

}; // class IImage

} // namespace Thea

#endif
