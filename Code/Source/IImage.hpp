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
 * a set of 2D images stacked along the depth dimension. Images are stored <b>bottom-up</b>: scanline 0 is the visually lowest
 * one.
 */
class THEA_API IImage
{
  public:
    THEA_DECL_SMART_POINTERS(IImage)

    /**
     * Different image types (enum class plus extra functions). The channel ordering within a pixel is exactly as in the type
     * name, e.g. RGB_8U has red as channel 0, green as channel 1, and blue as channel 2.
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
        LA_8U,          ///< Luminance + alpha image, 8 bits per channel.
        LA_16U,         ///< Luminance + alpha image, 16 bits per channel.
        LA_32F,         ///< Floating-point luminance + alpha image, 32 bits per channel.
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
       * Get the number of bits assigned to a particular channel. If the image doesn't contain the specific channel (e.g.
       * luminance images don't have red, green or blue channels) a value of zero is returned.
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

    /**
     * Get a pointer to the beginning of a specified row of pixels, optionally for a specific depth slice. Scanline 0 is the
     * <b>bottom/lowest</b> one.
     */
    virtual void const * THEA_ICALL getScanLine(int64 row, int64 z = 0) const = 0;

    /**
     * Get a pointer to the beginning of a specified row of pixels, optionally for a specific depth slice. Scanline 0 is the
     * <b>bottom/lowest</b> one.
     */
    virtual void * THEA_ICALL getScanLine(int64 row, int64 z = 0) = 0;

    /**
     * Get the distance in bytes between the start of one row of pixels and the start of the next row. Rows may be aligned to
     * 32-bit (or other) boundaries for performance reasons, so this is <b>not</b> necessarily equal to the number of pixels in
     * a row times the size of a pixel. Also note that the last row of pixels may not occupy a full stride.
     */
    virtual int64 THEA_ICALL getStrideBytes() const = 0;

}; // class IImage

} // namespace Thea

#endif
