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
class THEA_API Image : public virtual IImage, public Serializable
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
     * Get the number of bits assigned to a particular channel. If the image doesn't contain the specific channel (e.g.
     * luminance images don't have red, green or blue channels) a value of zero is returned.
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
     *   - For RGB images, passing \a channel = 3 returns 1, since the image is assumed to be non-transparent by common
     *     convention. In all other cases, passing an invalid channel index returns 0.
     *
     * This is a relatively slow way to iterate over pixel values and is provided only for convenience.
     *
     * This function <b>does not support channels smaller than a byte</b> (e.g. Type::LUMINANCE_1U), and returns a value of
     * zero for such images.
     */
    double getNormalizedValue(void const * pixel, int channel) const;

    /**
     * Reorder the channels of the image. \a order is a sequence of numChannels() indices that specifies the new channel
     * ordering. Duplicates are allowed. E.g.
     * \code
     * int order[] = { 2, 1, 0 }
     * \endcode
     * converts RGB to BGR, and
     * \code
     * int order[] = { 0, 0, 0 }
     * \endcode
     * converts an RGB image to greyscale by copying the red channel to green and blue.
     */
    bool reorderChannels(int const * order);

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
THEA_DEF_IMAGE_CODEC(CodecBmp,     "BMP",     "Windows or OS/2 Bitmap File (*.bmp)")
THEA_DEF_IMAGE_CODEC(CodecCut,     "CUT",     "Dr. Halo (*.cut)")
THEA_DEF_IMAGE_CODEC(CodecDds,     "DDS",     "DirectDraw Surface (*.dds)")
THEA_DEF_IMAGE_CODEC(CodecExr,     "EXR",     "ILM OpenEXR (*.exr)")
THEA_DEF_IMAGE_CODEC(CodecFaxg3,   "FAXG3",   "Raw Fax format CCITT G3 (*.g3)")
THEA_DEF_IMAGE_CODEC(CodecGif,     "GIF",     "Graphics Interchange Format (*.gif)")
THEA_DEF_IMAGE_CODEC(CodecHdr,     "HDR",     "High Dynamic Range (*.hdr)")
THEA_DEF_IMAGE_CODEC(CodecIco,     "ICO",     "Windows Icon (*.ico)")
THEA_DEF_IMAGE_CODEC(CodecIff,     "IFF",     "Amiga IFF (*.iff, *.lbm)")
THEA_DEF_IMAGE_CODEC(CodecJ2k,     "J2K",     "JPEG-2000 codestream (*.j2k, *.j2c)")
THEA_DEF_IMAGE_CODEC(CodecJng,     "JNG",     "JPEG Network Graphics (*.jng)")
THEA_DEF_IMAGE_CODEC(CodecJp2,     "JP2",     "JPEG-2000 File Format (*.jp2)")
THEA_DEF_IMAGE_CODEC(CodecKoala,   "KOALA",   "Commodore 64 Koala format (*.koa)")
THEA_DEF_IMAGE_CODEC(CodecMng,     "MNG",     "Multiple Network Graphics (*.mng)")
THEA_DEF_IMAGE_CODEC(CodecPbm,     "PBM",     "Portable Bitmap (ASCII) (*.pbm)")
THEA_DEF_IMAGE_CODEC(CodecPbmraw,  "PBMRAW",  "Portable Bitmap (BINARY) (*.pbm)")
THEA_DEF_IMAGE_CODEC(CodecPcd,     "PCD",     "Kodak PhotoCD (*.pcd)")
THEA_DEF_IMAGE_CODEC(CodecPcx,     "PCX",     "Zsoft Paintbrush PCX bitmap format (*.pcx)")
THEA_DEF_IMAGE_CODEC(CodecPfm,     "PFM",     "Portable Floatmap (*.pfm)")
THEA_DEF_IMAGE_CODEC(CodecPgm,     "PGM",     "Portable Graymap (ASCII) (*.pgm)")
THEA_DEF_IMAGE_CODEC(CodecPgmraw,  "PGMRAW",  "Portable Graymap (BINARY) (*.pgm)")
THEA_DEF_IMAGE_CODEC(CodecPng,     "PNG",     "Portable Network Graphics (*.png)")
THEA_DEF_IMAGE_CODEC(CodecPpm,     "PPM",     "Portable Pixelmap (ASCII) (*.ppm)")
THEA_DEF_IMAGE_CODEC(CodecPpmraw,  "PPMRAW",  "Portable Pixelmap (BINARY) (*.ppm)")
THEA_DEF_IMAGE_CODEC(CodecPsd,     "PSD",     "Adobe Photoshop (*.psd)")
THEA_DEF_IMAGE_CODEC(CodecRas,     "RAS",     "Sun Rasterfile (*.ras)")
THEA_DEF_IMAGE_CODEC(CodecSgi,     "SGI",     "Silicon Graphics SGI image format (*.sgi)")
THEA_DEF_IMAGE_CODEC(CodecTarga,   "TARGA",   "Truevision Targa files (*.tga, *.targa)")
THEA_DEF_IMAGE_CODEC(CodecTiff,    "TIFF",    "Tagged Image File Format (*.tif, *.tiff)")
THEA_DEF_IMAGE_CODEC(CodecWbmp,    "WBMP",    "Wireless Bitmap (*.wbmp)")
THEA_DEF_IMAGE_CODEC(CodecXbm,     "XBM",     "X11 Bitmap Format (*.xbm)")
THEA_DEF_IMAGE_CODEC(CodecXpm,     "XPM",     "X11 Pixmap Format (*.xpm)")

// 3D formats
THEA_DEF_IMAGE_CODEC(Codec3bm,     "3BM",     "3D Bitmap (*.3bm)")

/** JPEG image codec. */
class THEA_API CodecJpeg : public ImageCodec
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
    CodecJpeg() : options(Options::defaults()) {}

    /** Constructor to set encoding options. */
    CodecJpeg(Options const & options_) : options(options_) {}

    THEA_DEF_IMAGE_CODEC_BODY(CodecJpeg, "JPEG", "Independent JPEG Group (*.jpg, *.jif, *.jpeg, *.jpe)")

  private:
    Options options;
};

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Image)

#endif
