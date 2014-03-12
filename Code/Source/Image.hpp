//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Image_hpp__
#define __Thea_Image_hpp__

#include "Common.hpp"
#include "IOStream.hpp"
#include "Serializable.hpp"

// Forward declaration
class fipImage;

namespace Thea {

/** Abstract base class for images, useful for passing images across shared library boundaries. */
class THEA_API AbstractImage
{
  public:
    /**
     * Indices for the color channels in a single pixel. E.g. for a 24-bit RGB bitmap (8 bits per channel), the individual
     * components for the pixel at address <code>unsigned char * p</code> may be accessed as <code>p[Channel::RED]</code>,
     * <code>p[Channel::GREEN]</code> and <code>p[Channel::BLUE]</code>.
     */
    struct THEA_API Channel
    {
      static int const RED;    ///< Index of red channel
      static int const GREEN;  ///< Index of green channel
      static int const BLUE;   ///< Index of blue channel
      static int const ALPHA;  ///< Index of alpha channel
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
    virtual ~AbstractImage() {}

    /**
     * Check if the image has been allocated non-zero memory space (hence has valid type and dimensions) or not. An image
     * created by the default constructor is invalid and must be further initialized using deserialize() or a similar function.
     */
    virtual bool isValid() const = 0;

    /**
     * Destroy all image data, resetting the image to an invalid state.
     *
     * @see isValid()
     */
    virtual void clear() = 0;

    /** Resize the image, changing its type and dimensions. All previous image data is discarded. */
    virtual void resize(Type type, int width, int height) = 0;

    /** Get the width of the image in pixels. */
    virtual int getWidth() const = 0;

    /** Get the height of the image in pixels. */
    virtual int getHeight() const = 0;

    /** Get the type of the image pixels. */
    virtual Type getType() const = 0;

    /** Get a pointer to the image data. */
    virtual void const * getData() const = 0;

    /** Get a pointer to the image data. */
    virtual void * getData() = 0;

    /** Get a pointer to the beginning of a specified row of pixels. */
    virtual void const * getScanLine(int row) const = 0;

    /** Get a pointer to the beginning of a specified row of pixels. */
    virtual void * getScanLine(int row) = 0;

    /**
     * Get the number of bytes consumed by a row of pixels. Rows may be aligned to 32-bit (or other) boundaries for performance
     * reasons, so this is <b>not</b> necessarily equal to the number of pixels in a row times the size of a pixel.
     */
    virtual int getScanWidth() const = 0;

    /** Get the byte alignment of a pixel row. */
    virtual int getRowAlignment() const = 0;

}; // class AbstractImage

/** 2D images. */
class THEA_API Image : public AbstractImage, public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(Image, shared_ptr, weak_ptr)

    /** Construct an empty image with no initial data. */
    Image();

    /** Construct an uninitialized image of the specified type and pixel dimensions, which must have valid non-zero values. */
    Image(Type type_, int width, int height);

    /**
     * Construct an image by deserializing it from an input stream.
     *
     * @see deserialize()
     */
    Image(BinaryInputStream & input, Codec const & codec = Codec_AUTO());

    /**
     * Construct an image by loading it from a file.
     *
     * @see load()
     */
    Image(std::string const & filename, Codec const & codec = Codec_AUTO());

    /* Copy constructor. */
    Image(Image const & src);

    /** Destructor. */
    ~Image();

    /* Assignment operator. */
    Image & operator=(Image const & src);

    bool isValid() const;
    void clear();
    void resize(Type type, int width, int height);
    int getWidth() const;
    int getHeight() const;
    Type getType() const { return type; }

    /**
     * Get the number of channels per pixel. For example: 1 for luminance, 3 for RGB, 4 for RGBA. A complex number corresponds
     * to a single channel (not two). Returns -1 if the image is of unknown type.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    int numChannels() const { return num_channels; }

    /**
     * Check if the channels hold complex (as opposed to real) values.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    bool isComplex() const { return is_complex; }

    /**
     * Check if the channels hold floating-point values (at any precision).
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    bool isFloatingPoint() const { return is_floating_point; }

    /**
     * Get the number of bits assigned to each pixel.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    int getBitsPerPixel() const { return bits_per_pixel; }

    /**
     * Get the number of bits assigned to each channel. Returns -1 if the channels don't all have the same number of bits.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    int getBitsPerChannel() const { return bits_per_channel; }

    /**
     * Get the number of bits assigned to a particular channel. For luminance images, the single channel is assumed to
     * correspond to the enum value Channel::ALPHA. If the image doesn't contain the specific channel (e.g. luminance images
     * don't have red, green or blue channels) a value of zero is returned.
     */
    int getBitsInChannel(int channel) const { return type.getBitsInChannel(channel); }

    /**
     * Check if the image pixels start at byte addresses. For example, RGB_8U does, but LUMINANCE_1U does not.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    bool hasByteAlignedPixels() const { return has_byte_aligned_pixels; }

    /**
     * Check if all pixel channels are aligned to byte addresses.
     *
     * This returns a cached value and is likely to be faster than calling the equivalent function in Image::Type.
     */
    bool hasByteAlignedChannels() const { return has_byte_aligned_channels; }

    void const * getData() const;
    void * getData();
    void const * getScanLine(int row) const;
    void * getScanLine(int row);
    int getScanWidth() const;
    int getRowAlignment() const;

    /**
     * Get the value of a channel of a particular pixel, normalized to the range [0, 1]. The following caveats should be
     * noted:
     *   - Signed integer channels are scaled to the range [-1, 1). E.g. Type::LUMINANCE_16 natively stores values in the range
     *     [-32768, 32767], which is mapped to [1, 1) by dividing by 32768.
     *   - Floating point channels are assumed to be pre-normalized and no further normalization is done.
     *   - For complex channels, the function returns the magnitude of the complex number.
     *   - For single-channel images, the luminance is extracted with Channel::ALPHA.
     *
     * This is a relatively slow way to iterate over pixel values and is provided only for convenience.
     *
     * This function <b>does not support channels smaller than a byte</b> (e.g. Type::LUMINANCE_1U), and returns a value of
     * zero for such images.
     */
    double getNormalizedValue(void const * pixel, int channel) const;

    /** Invert the pixel values of the image. */
    bool invert();

    /**
     * Convert this image to a different format. Currently only source/target format combinations supported by FreeImage are
     * available.
     */
    bool convert(Type dst_type);

    /**
     * Convert an image from one format to another. Currently only source/target format combinations supported by FreeImage are
     * available.
     */
    bool convert(Type dst_type, Image & dst) const;

    /** Rescale the image to a new width and height. */
    bool rescale(int new_width, int new_height, Filter filter = Filter::AUTO);

    /**
     * {@inheritDoc}
     *
     * The serialized image is prefixed with a header indicating its encoded size.
     */
    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const;

    /**
     * {@inheritDoc}
     *
     * The image <b>must</b> have been serialized using the layout specified in serialize().
     */
    void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO());

    /**
     * Save the image to an image file. Unlike serialize(), the file will <b>not</b> have a prefixed header. An exception will
     * be thrown if the image cannot be saved.
     */
    void save(std::string const & filename, Codec const & codec = Codec_AUTO()) const;

    /**
     * Load the image from an image file. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception
     * will be thrown if the image cannot be loaded.
     */
    void load(std::string const & filename, Codec const & codec = Codec_AUTO());

    /** <b>[Internal use only]</b> Get the wrapped FreeImage bitmap. */
    fipImage const * _getFreeImage() const { return fip_img; }

    /** <b>[Internal use only]</b> Get the wrapped FreeImage bitmap. */
    fipImage * _getFreeImage() { return fip_img; }

    /** <b>[Internal use only]</b> Set the type of the image. */
    void _setType(Type type_);

  private:
    /** Automatically detect the type of the encoded image and deserialize it appropriately. */
    void deserialize_AUTO(BinaryInputStream & input, bool read_prefixed_info);

    /** Cache properties related to the image type. */
    void cacheTypeProperties();

    fipImage * fip_img;
    Type type;

    // Cached type properties for fast access
    int   num_channels;
    bool  is_complex;
    bool  is_floating_point;
    int   bits_per_pixel;
    int   bits_per_channel;
    bool  has_byte_aligned_pixels;
    bool  has_byte_aligned_channels;
};

/** Abstract base class for all image codecs. */
class THEA_API ImageCodec : public Codec
{
  public:
    /**
     * Serialize an image to a binary output stream. Optionally prefixes extra information about the image block such as its
     * size and type (which may have not been specified in the encoding format itself).
     *
     * @return The number of bytes written (this will include the information field if it was written)
     */
    virtual long serializeImage(Image const & image, BinaryOutputStream & output, bool prefix_info) const = 0;

    /**
     * Deserialize an image from a binary input stream. If the <code>read_prefixed_info</code> parameter is true, extra
     * information about the image block (such as its size and type) will be read first from the input stream. Else, the entire
     * input will be treated as the image block (the size() function of the stream must return the correct value in this case).
     *
     * @see serializeMeshGroup
     */
    virtual void deserializeImage(Image & image, BinaryInputStream & input, bool read_prefixed_info) const = 0;
};

#define THEA_DEF_IMAGE_CODEC_BODY(name, desc)                                                                                 \
    public:                                                                                                                   \
      char const * getName() const { static char const * my_name = desc; return my_name; }                                    \
      long serializeImage(Image const & image, BinaryOutputStream & output, bool prefix_info) const;                          \
      void deserializeImage(Image & image, BinaryInputStream & input, bool read_prefixed_info) const;

#define THEA_DEF_IMAGE_CODEC(name, desc)                                                                                      \
  class THEA_API name : public ImageCodec                                                                                     \
  {                                                                                                                           \
    THEA_DEF_IMAGE_CODEC_BODY(name, desc)                                                                                     \
  };

// TODO: Add options to all the ones that support them
THEA_DEF_IMAGE_CODEC(CodecBMP,     "Windows or OS/2 Bitmap File (*.BMP)")
THEA_DEF_IMAGE_CODEC(CodecCUT,     "Dr. Halo (*.CUT)")
THEA_DEF_IMAGE_CODEC(CodecDDS,     "DirectDraw Surface (*.DDS)")
THEA_DEF_IMAGE_CODEC(CodecEXR,     "ILM OpenEXR (*.EXR)")
THEA_DEF_IMAGE_CODEC(CodecFAXG3,   "Raw Fax format CCITT G3 (*.G3)")
THEA_DEF_IMAGE_CODEC(CodecGIF,     "Graphics Interchange Format (*.GIF)")
THEA_DEF_IMAGE_CODEC(CodecHDR,     "High Dynamic Range (*.HDR)")
THEA_DEF_IMAGE_CODEC(CodecICO,     "Windows Icon (*.ICO)")
THEA_DEF_IMAGE_CODEC(CodecIFF,     "Amiga IFF (*.IFF, *.LBM)")
THEA_DEF_IMAGE_CODEC(CodecJ2K,     "JPEG-2000 codestream (*.J2K, *.J2C)")
THEA_DEF_IMAGE_CODEC(CodecJNG,     "JPEG Network Graphics (*.JNG)")
THEA_DEF_IMAGE_CODEC(CodecJP2,     "JPEG-2000 File Format (*.JP2)")
THEA_DEF_IMAGE_CODEC(CodecKOALA,   "Commodore 64 Koala format (*.KOA)")
THEA_DEF_IMAGE_CODEC(CodecMNG,     "Multiple Network Graphics (*.MNG)")
THEA_DEF_IMAGE_CODEC(CodecPBM,     "Portable Bitmap (ASCII) (*.PBM)")
THEA_DEF_IMAGE_CODEC(CodecPBMRAW,  "Portable Bitmap (BINARY) (*.PBM)")
THEA_DEF_IMAGE_CODEC(CodecPCD,     "Kodak PhotoCD (*.PCD)")
THEA_DEF_IMAGE_CODEC(CodecPCX,     "Zsoft Paintbrush PCX bitmap format (*.PCX)")
THEA_DEF_IMAGE_CODEC(CodecPFM,     "Portable Floatmap (*.PFM)")
THEA_DEF_IMAGE_CODEC(CodecPGM,     "Portable Graymap (ASCII) (*.PGM)")
THEA_DEF_IMAGE_CODEC(CodecPGMRAW,  "Portable Graymap (BINARY) (*.PGM)")
THEA_DEF_IMAGE_CODEC(CodecPNG,     "Portable Network Graphics (*.PNG)")
THEA_DEF_IMAGE_CODEC(CodecPPM,     "Portable Pixelmap (ASCII) (*.PPM)")
THEA_DEF_IMAGE_CODEC(CodecPPMRAW,  "Portable Pixelmap (BINARY) (*.PPM)")
THEA_DEF_IMAGE_CODEC(CodecPSD,     "Adobe Photoshop (*.PSD)")
THEA_DEF_IMAGE_CODEC(CodecRAS,     "Sun Rasterfile (*.RAS)")
THEA_DEF_IMAGE_CODEC(CodecSGI,     "Silicon Graphics SGI image format (*.SGI)")
THEA_DEF_IMAGE_CODEC(CodecTARGA,   "Truevision Targa files (*.TGA, *.TARGA)")
THEA_DEF_IMAGE_CODEC(CodecTIFF,    "Tagged Image File Format (*.TIF, *.TIFF)")
THEA_DEF_IMAGE_CODEC(CodecWBMP,    "Wireless Bitmap (*.WBMP)")
THEA_DEF_IMAGE_CODEC(CodecXBM,     "X11 Bitmap Format (*.XBM)")
THEA_DEF_IMAGE_CODEC(CodecXPM,     "X11 Pixmap Format (*.XPM)")

/** JPEG image codec. */
class THEA_API CodecJPEG : public ImageCodec
{
  public:
    /** %Options for JPEG encoding. */
    struct THEA_API Options
    {
      int quality;
      bool progressive;

      /** Constructor. */
      Options(int quality_, bool progressive_) : quality(quality_), progressive(progressive_) {}

      /** Copy constructor. */
      Options(Options const & src) : quality(src.quality), progressive(src.progressive) {}

      /** The set of default options (quality 75, non-progressive encoding). */
      static Options const & defaults() { static Options const def(75, false); return def; }
    };

    /** Default constructor. Sets default options (quality 75, non-progressive encoding). */
    CodecJPEG() : options(Options::defaults()) {}

    /** Constructor to set encoding options. */
    CodecJPEG(Options const & options_) : options(options_) {}

    THEA_DEF_IMAGE_CODEC_BODY(CodecJPEG, "Independent JPEG Group (*.JPG, *.JIF, *.JPEG, *.JPE)")

  private:
    Options options;
};

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Image)

#endif
