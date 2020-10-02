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

#ifndef __Thea_Image_hpp__
#define __Thea_Image_hpp__

#include "Common.hpp"
#include "IImage.hpp"
#include "AlignedAllocator.hpp"
#include "Array.hpp"
#include "Iostream.hpp"
#include "Serializable.hpp"

// Forward declaration
class fipImage;

namespace Thea {

/** A 2D image. */
class THEA_API Image : public IImage, public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Image)

    /** Construct an empty image with no initial data. */
    Image();

    /** Construct an uninitialized image of the specified type and pixel dimensions, which must have valid non-zero values. */
    Image(Type type_, int64 width_, int64 height_, int64 depth_ = 1);

    /**
     * Construct an image by deserializing it from an input stream.
     *
     * @see read()
     */
    Image(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);

    /**
     * Construct an image by loading it from a file.
     *
     * @see load()
     */
    Image(std::string const & path, Codec const & codec = CodecAuto());

    /* Copy constructor. */
    Image(Image const & src);

    /** Destructor. */
    ~Image();

    /* Assignment operator. */
    Image & operator=(Image const & src);

    // Abstract interface functions
    int8 THEA_ICALL isValid() const;
    int8 THEA_ICALL clear();
    int8 THEA_ICALL resize(int64 type, int64 width_, int64 height_, int64 depth_ = 1);
    int64 THEA_ICALL getWidth() const { return width; }
    int64 THEA_ICALL getHeight() const { return height; }
    int64 THEA_ICALL getDepth() const { return depth; }
    int32 THEA_ICALL getType() const { return (int)type; }

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

    // Abstract interface functions
    void const * THEA_ICALL getData() const;
    void * THEA_ICALL getData();
    void const * THEA_ICALL getScanLine(int64 row, int64 z = 0) const;
    void * THEA_ICALL getScanLine(int64 row, int64 z = 0);
    int64 THEA_ICALL getScanWidth() const;
    int32 THEA_ICALL getRowAlignment() const;

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
    bool rescale(int64 new_width, int64 new_height, int64 new_depth = 1, Filter filter = Filter::AUTO);

    void THEA_ICALL read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);
    void THEA_ICALL write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const;

    /**
     * Load the image from an image file. The file should <b>not</b> have a prefixed Codec::BlockHeader. An exception will be
     * thrown if the image cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec = CodecAuto());

    /**
     * Save the image to an image file. The file will <b>not</b> have a prefixed Codec::BlockHeader. An exception will be thrown
     * if the image cannot be saved.
     */
    void save(std::string const & path, Codec const & codec = CodecAuto()) const;

    /** <b>[Internal use only]</b> Get the wrapped FreeImage bitmap. */
    fipImage const * _getFreeImage() const { return fip_img; }

    /** <b>[Internal use only]</b> Get the wrapped FreeImage bitmap. */
    fipImage * _getFreeImage() { return fip_img; }

    /** <b>[Internal use only]</b> Set the type of the image. */
    void _setType(Type type_);

  private:
    /** Automatically detect the type of the encoded image and read it appropriately. */
    void readAuto(BinaryInputStream & input, bool read_block_header);

    /** Cache properties related to the image type. */
    void cacheTypeProperties();

    // Scanline alignment when allocating custom arrays
    static size_t const ROW_ALIGNMENT = 8;  // would prefer 16 for SSE compatibility, but OpenGL supports a max of 8

    // Image parameters
    Type type;
    int64 width;
    int64 height;
    int64 depth;

    // Image data
    fipImage * fip_img;
    Array< uint8, AlignedAllocator<uint8, ROW_ALIGNMENT> > data;  // pixel buffer when fipImage won't work, e.g. 3D images

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
     * Read an image from a binary input stream. If \a read_block_header is true, extra information about the image block (such
     * as its size and type) will be read first from the input stream. Else, the entire input will be treated as the image block
     * (the size() function of the stream must return the correct value in this case).
     *
     * @see writeImage
     */
    virtual void readImage(Image & image, BinaryInputStream & input, bool read_block_header) const = 0;

    /**
     * Write an image to a binary output stream. Optionally prefixes extra information about the image block such as its size
     * and type (which may have not been specified in the encoding format itself).
     *
     * @see readImage
     */
    virtual void writeImage(Image const & image, BinaryOutputStream & output, bool write_block_header) const = 0;
};

#define THEA_DEF_IMAGE_CODEC_BODY(name, magic, desc)                                                                          \
    public:                                                                                                                   \
      char const * getName() const { static char const * my_name = desc; return my_name; }                                    \
      MagicString const & getMagic() const { static MagicString const magic_ = toMagic(magic); return magic_; }               \
      void readImage(Image & image, BinaryInputStream & input, bool read_block_header) const;                                 \
      void writeImage(Image const & image, BinaryOutputStream & output, bool write_block_header) const;

#define THEA_DEF_IMAGE_CODEC(name, magic, desc)                                                                               \
  class THEA_API name : public ImageCodec                                                                                     \
  {                                                                                                                           \
    public: name() {}                                                                                                         \
    THEA_DEF_IMAGE_CODEC_BODY(name, magic, desc)                                                                              \
  };

// TODO: Add options to all the ones that support them

// 2D formats
THEA_DEF_IMAGE_CODEC(CodecBMP,     "BMP",     "Windows or OS/2 Bitmap File (*.BMP)")
THEA_DEF_IMAGE_CODEC(CodecCUT,     "CUT",     "Dr. Halo (*.CUT)")
THEA_DEF_IMAGE_CODEC(CodecDDS,     "DDS",     "DirectDraw Surface (*.DDS)")
THEA_DEF_IMAGE_CODEC(CodecEXR,     "EXR",     "ILM OpenEXR (*.EXR)")
THEA_DEF_IMAGE_CODEC(CodecFAXG3,   "FAXG3",   "Raw Fax format CCITT G3 (*.G3)")
THEA_DEF_IMAGE_CODEC(CodecGIF,     "GIF",     "Graphics Interchange Format (*.GIF)")
THEA_DEF_IMAGE_CODEC(CodecHDR,     "HDR",     "High Dynamic Range (*.HDR)")
THEA_DEF_IMAGE_CODEC(CodecICO,     "ICO",     "Windows Icon (*.ICO)")
THEA_DEF_IMAGE_CODEC(CodecIFF,     "IFF",     "Amiga IFF (*.IFF, *.LBM)")
THEA_DEF_IMAGE_CODEC(CodecJ2K,     "J2K",     "JPEG-2000 codestream (*.J2K, *.J2C)")
THEA_DEF_IMAGE_CODEC(CodecJNG,     "JNG",     "JPEG Network Graphics (*.JNG)")
THEA_DEF_IMAGE_CODEC(CodecJP2,     "JP2",     "JPEG-2000 File Format (*.JP2)")
THEA_DEF_IMAGE_CODEC(CodecKOALA,   "KOALA",   "Commodore 64 Koala format (*.KOA)")
THEA_DEF_IMAGE_CODEC(CodecMNG,     "MNG",     "Multiple Network Graphics (*.MNG)")
THEA_DEF_IMAGE_CODEC(CodecPBM,     "PBM",     "Portable Bitmap (ASCII) (*.PBM)")
THEA_DEF_IMAGE_CODEC(CodecPBMRAW,  "PBMRAW",  "Portable Bitmap (BINARY) (*.PBM)")
THEA_DEF_IMAGE_CODEC(CodecPCD,     "PCD",     "Kodak PhotoCD (*.PCD)")
THEA_DEF_IMAGE_CODEC(CodecPCX,     "PCX",     "Zsoft Paintbrush PCX bitmap format (*.PCX)")
THEA_DEF_IMAGE_CODEC(CodecPFM,     "PFM",     "Portable Floatmap (*.PFM)")
THEA_DEF_IMAGE_CODEC(CodecPGM,     "PGM",     "Portable Graymap (ASCII) (*.PGM)")
THEA_DEF_IMAGE_CODEC(CodecPGMRAW,  "PGMRAW",  "Portable Graymap (BINARY) (*.PGM)")
THEA_DEF_IMAGE_CODEC(CodecPNG,     "PNG",     "Portable Network Graphics (*.PNG)")
THEA_DEF_IMAGE_CODEC(CodecPPM,     "PPM",     "Portable Pixelmap (ASCII) (*.PPM)")
THEA_DEF_IMAGE_CODEC(CodecPPMRAW,  "PPMRAW",  "Portable Pixelmap (BINARY) (*.PPM)")
THEA_DEF_IMAGE_CODEC(CodecPSD,     "PSD",     "Adobe Photoshop (*.PSD)")
THEA_DEF_IMAGE_CODEC(CodecRAS,     "RAS",     "Sun Rasterfile (*.RAS)")
THEA_DEF_IMAGE_CODEC(CodecSGI,     "SGI",     "Silicon Graphics SGI image format (*.SGI)")
THEA_DEF_IMAGE_CODEC(CodecTARGA,   "TARGA",   "Truevision Targa files (*.TGA, *.TARGA)")
THEA_DEF_IMAGE_CODEC(CodecTIFF,    "TIFF",    "Tagged Image File Format (*.TIF, *.TIFF)")
THEA_DEF_IMAGE_CODEC(CodecWBMP,    "WBMP",    "Wireless Bitmap (*.WBMP)")
THEA_DEF_IMAGE_CODEC(CodecXBM,     "XBM",     "X11 Bitmap Format (*.XBM)")
THEA_DEF_IMAGE_CODEC(CodecXPM,     "XPM",     "X11 Pixmap Format (*.XPM)")

// 3D formats
THEA_DEF_IMAGE_CODEC(Codec3BM,     "3BM",     "3D Bitmap (*.3BM)")

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

    THEA_DEF_IMAGE_CODEC_BODY(CodecJPEG, "JPEG", "Independent JPEG Group (*.JPG, *.JIF, *.JPEG, *.JPE)")

  private:
    Options options;
};

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Image)

#endif
