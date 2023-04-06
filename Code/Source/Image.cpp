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
// =====
//
// This product uses software from the G3D project (http://g3d-cpp.sf.net).
// Copyright (c) 2000-2008, Morgan McGuire
// For the full G3D license see LICENSE.txt in the documentation.
//
// =====
//
// This software is based in part on the work of the Independent JPEG Group.
//
//============================================================================

#include "Image.hpp"
#include "Array.hpp"
#include "FilePath.hpp"
#include "UnorderedMap.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>

#if THEA_ENABLE_FREEIMAGE
#  include <FreeImagePlus.h>
   // Work around a bug where FreeImage.h defines _WINDOWS_ and this messes with platform detection in other libs e.g. wxWindows
#  if !THEA_WINDOWS
#    ifdef _WINDOWS_
#      undef _WINDOWS_
#    endif
#  endif
#else
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
#    define STBI_NO_STDIO
#    define STB_IMAGE_IMPLEMENTATION
#    define STB_IMAGE_RESIZE_IMPLEMENTATION
#    define STB_IMAGE_WRITE_IMPLEMENTATION
#    include "ThirdParty/stb/stb_image.h"
#    include "ThirdParty/stb/stb_image_resize.h"
#    include "ThirdParty/stb/stb_image_write.h"
#  pragma clang diagnostic pop
#endif

THEA_INSTANTIATE_SMART_POINTERS(Thea::Image)

namespace Thea {

#if THEA_ENABLE_FREEIMAGE

#  define THEA_DEF_IMAGE_CODEC_FREEIMAGE(name, fif, read_flags, write_flags)                                                  \
      int name::getImplFormat() const { return (fif); }                                                                       \
      int name::getImplReadFlags() const { return (read_flags); }                                                             \
      int name::getImplWriteFlags() const { return (write_flags); }

THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecBmp,    FIF_BMP,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecCut,    FIF_CUT,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecDds,    FIF_DDS,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecExr,    FIF_EXR,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecFaxg3,  FIF_FAXG3,    0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecGif,    FIF_GIF,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecHdr,    FIF_HDR,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecIco,    FIF_ICO,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecIff,    FIF_IFF,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecJ2k,    FIF_J2K,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecJng,    FIF_JNG,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecJp2,    FIF_JP2,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecJpg,    FIF_JPEG,     0, (options.quality | (options.progressive ? JPEG_PROGRESSIVE : 0)))
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecKoa,    FIF_KOALA,    0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecMng,    FIF_MNG,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPbm,    (binary ? FIF_PBMRAW : FIF_PBM), 0, 0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPcd,    FIF_PCD,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPcx,    FIF_PCX,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPfm,    FIF_PFM,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPgm,    (binary ? FIF_PGMRAW : FIF_PGM), 0, 0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPng,    FIF_PNG,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPpm,    (binary ? FIF_PPMRAW : FIF_PPM), 0, 0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecPsd,    FIF_PSD,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecRas,    FIF_RAS,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecSgi,    FIF_SGI,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecTga,    FIF_TARGA,    0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecTif,    FIF_TIFF,     0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecWbmp,   FIF_WBMP,     0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecXbm,    FIF_XBM,      0,  0)
THEA_DEF_IMAGE_CODEC_FREEIMAGE(CodecXpm,    FIF_XPM,      0,  0)

THEA_DEF_IMAGE_CODEC_FREEIMAGE(Codec3bm,    FIF_UNKNOWN,  0,  0)

#endif

namespace ImageInternal {

// Bytes covered by one or more pixels (including partial bytes)
int64
pixelBytes(int bpp, int64 num_pixels = 1)
{
  auto bits = num_pixels * (int64)bpp;
  return bits / 8 + (bits % 8 == 0 ? 0 : 1);
}

// Bytes consumed by an aligned row
int64
strideBytes(int64 width, int bpp, int32 alignment)
{
  auto row_bytes = pixelBytes(width, bpp);
  return row_bytes + (alignment - (row_bytes % alignment)) % alignment;
}

// Bits (not bytes) per pixel
static int
typeToBpp(Image::Type type)
{
  switch (type)
  {
    case Image::Type::LUMINANCE_1U  : return 1;
    case Image::Type::LUMINANCE_2U  : return 2;
    case Image::Type::LUMINANCE_4U  : return 4;
    case Image::Type::LUMINANCE_8U  : return 8;
    case Image::Type::LUMINANCE_16  :
    case Image::Type::LUMINANCE_16U : return 16;
    case Image::Type::LUMINANCE_32  :
    case Image::Type::LUMINANCE_32U : return 32;
    case Image::Type::LUMINANCE_32F : return 32;
    case Image::Type::LUMINANCE_64F : return 64;
    case Image::Type::LA_8U         : return 16;
    case Image::Type::LA_16U        : return 32;
    case Image::Type::LA_32F        : return 64;
    case Image::Type::RGB_8U        : return 24;
    case Image::Type::RGBA_8U       : return 32;
    case Image::Type::RGB_16U       : return 48;
    case Image::Type::RGBA_16U      : return 64;
    case Image::Type::RGB_32F       : return 96;
    case Image::Type::RGBA_32F      : return 128;
    case Image::Type::COMPLEX_64F   : return 128;
    default                         : return 0;
  }
}

#if THEA_ENABLE_FREEIMAGE

static FREE_IMAGE_TYPE
typeToFreeImageType(Image::Type type)
{
  switch (type)
  {
    case Image::Type::LUMINANCE_1U  :
    case Image::Type::LUMINANCE_2U  :
    case Image::Type::LUMINANCE_4U  :
    case Image::Type::LUMINANCE_8U  : return FIT_BITMAP;
    case Image::Type::LUMINANCE_16  : return FIT_INT16;
    case Image::Type::LUMINANCE_16U : return FIT_UINT16;
    case Image::Type::LUMINANCE_32  : return FIT_INT32;
    case Image::Type::LUMINANCE_32U : return FIT_UINT32;
    case Image::Type::LUMINANCE_32F : return FIT_FLOAT;
    case Image::Type::LUMINANCE_64F : return FIT_DOUBLE;
    case Image::Type::RGB_8U        : return FIT_BITMAP;
    case Image::Type::RGBA_8U       : return FIT_BITMAP;
    case Image::Type::RGB_16U       : return FIT_RGB16;
    case Image::Type::RGBA_16U      : return FIT_RGBA16;
    case Image::Type::RGB_32F       : return FIT_RGBF;
    case Image::Type::RGBA_32F      : return FIT_RGBAF;
    case Image::Type::COMPLEX_64F   : return FIT_COMPLEX;
    default                         : return FIT_UNKNOWN;
  }
}

static Image::Type
typeFromFreeImageTypeAndBpp(FREE_IMAGE_TYPE fi_type, WORD fi_bpp)
{
  switch (fi_type)
  {
    case FIT_BITMAP:
    {
      switch (fi_bpp)
      {
        case 1:  return Image::Type::LUMINANCE_1U;
        case 2:  return Image::Type::LUMINANCE_2U;
        case 4:  return Image::Type::LUMINANCE_4U;
        case 8:  return Image::Type::LUMINANCE_8U;
        case 24: return Image::Type::RGB_8U;
        case 32: return Image::Type::RGBA_8U;
      }
      break;
    }
    case FIT_INT16:   return Image::Type::LUMINANCE_16;
    case FIT_UINT16:  return Image::Type::LUMINANCE_16U;
    case FIT_INT32:   return Image::Type::LUMINANCE_32;
    case FIT_UINT32:  return Image::Type::LUMINANCE_32U;
    case FIT_FLOAT:   return Image::Type::LUMINANCE_32F;
    case FIT_DOUBLE:  return Image::Type::LUMINANCE_64F;
    case FIT_RGB16:   return Image::Type::RGB_16U;
    case FIT_RGBA16:  return Image::Type::RGBA_16U;
    case FIT_RGBF:    return Image::Type::RGB_32F;
    case FIT_RGBAF:   return Image::Type::RGBA_32F;
    case FIT_COMPLEX: return Image::Type::COMPLEX_64F;
    default:          break;  // to eliminate a GCC warning about unhandled FIT_UNKNOWN
  }

  return Image::Type::UNKNOWN;
}

#else // THEA_ENABLE_FREEIMAGE

static Image::Type
typeFromStbChannelsAndBpc(int nc, int bpc)
{
  switch (bpc)
  {
    case 16:
      switch (nc)
      {
        case 1:  return Image::Type::LUMINANCE_16U;
        case 2:  return Image::Type::LA_16U;
        case 3:  return Image::Type::RGB_16U;
        case 4:  return Image::Type::RGBA_16U;
        default: return Image::Type::UNKNOWN;
      }

    case 32:
      switch (nc)
      {
        case 1:  return Image::Type::LUMINANCE_32F;
        case 2:  return Image::Type::LA_32F;
        case 3:  return Image::Type::RGB_32F;
        case 4:  return Image::Type::RGBA_32F;
        default: return Image::Type::UNKNOWN;
      }

    default:  // 8-bit
      switch (nc)
      {
        case 1:  return Image::Type::LUMINANCE_8U;
        case 2:  return Image::Type::LA_8U;
        case 3:  return Image::Type::RGB_8U;
        case 4:  return Image::Type::RGBA_8U;
        default: return Image::Type::UNKNOWN;
      }
  }
}

#endif // THEA_ENABLE_FREEIMAGE

// Only called for saving images to files when a codec should be auto-inferred from the filename
static ImageCodec const *
codecFromPath(std::string const & path)
{
  typedef std::shared_ptr<ImageCodec> ICPtr;
  typedef UnorderedMap<std::string, ICPtr > CodecMap;

  static CodecMap codec_map;
  if (codec_map.empty())  // only on the first call
  {
    // 2D formats

    codec_map["BMP"  ] = ICPtr(new CodecBmp());
    codec_map["GIF"  ] = ICPtr(new CodecGif());
    codec_map["HDR"  ] = ICPtr(new CodecHdr());
    codec_map["JPG"  ] =
    codec_map["JIF"  ] =
    codec_map["JPEG" ] =
    codec_map["JPE"  ] = ICPtr(new CodecJpg());
    codec_map["PGM"  ] = ICPtr(new CodecPgm());
    codec_map["PNG"  ] = ICPtr(new CodecPng());
    codec_map["PPM"  ] = ICPtr(new CodecPpm());
    codec_map["PSD"  ] = ICPtr(new CodecPsd());
    codec_map["TGA"  ] =
    codec_map["TARGA"] = ICPtr(new CodecTga());

#if THEA_ENABLE_FREEIMAGE

    codec_map["CUT"  ] = ICPtr(new CodecCut());
    codec_map["DDS"  ] = ICPtr(new CodecDds());
    codec_map["EXR"  ] = ICPtr(new CodecExr());
    codec_map["ICO"  ] = ICPtr(new CodecIco());
    codec_map["IFF"  ] =
    codec_map["LBM"  ] = ICPtr(new CodecIff());
    codec_map["J2K"  ] =
    codec_map["J2C"  ] = ICPtr(new CodecJ2k());
    codec_map["JNG"  ] = ICPtr(new CodecJng());
    codec_map["JP2"  ] = ICPtr(new CodecJp2());
    codec_map["KOA"  ] = ICPtr(new CodecKoa());
    codec_map["MNG"  ] = ICPtr(new CodecMng());
    codec_map["PBM"  ] = ICPtr(new CodecPbm());
    codec_map["PCD"  ] = ICPtr(new CodecPcd());
    codec_map["PCX"  ] = ICPtr(new CodecPcx());
    codec_map["PFM"  ] = ICPtr(new CodecPfm());
    codec_map["RAS"  ] = ICPtr(new CodecRas());
    codec_map["SGI"  ] = ICPtr(new CodecSgi());
    codec_map["TIF"  ] =
    codec_map["TIFF" ] = ICPtr(new CodecTif());
    codec_map["WBMP" ] = ICPtr(new CodecWbmp());
    codec_map["XBM"  ] = ICPtr(new CodecXbm());
    codec_map["XPM"  ] = ICPtr(new CodecXpm());

#else // THEA_ENABLE_FREEIMAGE

    codec_map["PIC"  ] = ICPtr(new CodecPic());

#endif // THEA_ENABLE_FREEIMAGE

    // 3D formats
    codec_map["3BM"  ] = ICPtr(new Codec3bm());
  }

  CodecMap::const_iterator existing = codec_map.find(toUpper(FilePath::extension(path)));
  if (existing == codec_map.end())
    return nullptr;
  else
    return existing->second.get();
}

// Currently only detects 3D formats, the rest are handled by FreeImage
ImageCodec const *
codecFromMagic(int64 num_bytes, uint8 const * buf)
{
  static Codec3bm const CODEC_3BM;

  if (num_bytes >= 4 && buf[0] == (uint8)'3' && buf[1] == (uint8)'B' && buf[2] == (uint8)'M' && buf[3] == (uint8)'\0')
    return &CODEC_3BM;

  return nullptr;
}

#if THEA_ENABLE_FREEIMAGE

// FreeImage can store 8-bit channels in arbitrary order, accessed via FI_RGBA_RED, FI_RGBA_GREEN etc indices. This function
// reorders each pixel in-place so that the order is RGBA.
bool
canonicalizeChannelOrder(Image & img)
{
  auto bpc = img.getBitsPerChannel();
  auto nc = img.numChannels();
  if (bpc == 8 && nc > 1)
  {
    static int const CANONICAL_ORDER[4] = { FI_RGBA_RED, FI_RGBA_GREEN, FI_RGBA_BLUE, FI_RGBA_ALPHA };
    return img.reorderChannels(CANONICAL_ORDER);
  }

  return true;
}

// Undo the operation of canonicalizeChannelOrder()
bool
decanonicalizeChannelOrder(Image & img)
{
  auto bpc = img.getBitsPerChannel();
  auto nc = img.numChannels();
  if (bpc == 8 && nc > 1)
  {
    int decanonical_order[4];
    decanonical_order[FI_RGBA_RED  ] = 0;
    decanonical_order[FI_RGBA_GREEN] = 1;
    decanonical_order[FI_RGBA_BLUE ] = 2;
    decanonical_order[FI_RGBA_ALPHA] = 3;
    return img.reorderChannels(decanonical_order);
  }

  return true;
}

#else // THEA_ENABLE_FREEIMAGE

// Fill 'data' with 'size' bytes. Return number of bytes actually read.
int
stbStreamRead(void * user, char * data, int size)
{
  auto bis = (BinaryInputStream *)user;
  if (!bis) { return 0; }

  int64 num_read = 0;
  bis->readBytes(size, data, num_read);

  return (int)num_read;
}

// Skip the next 'n' bytes, or 'unget' the last -n bytes if negative.
void
stbStreamSkip(void * user, int n)
{
  auto bis = (BinaryInputStream *)user;
  if (!bis) { return; }

  bis->skip(n);
}

// Returns nonzero if we are at end of file/data.
int
stbStreamEof(void * user)
{
  auto bis = (BinaryInputStream *)user;
  if (!bis) { return 1; }

  return !bis->hasMore();
}

// Write a chunk of data to a stream.
void
stbStreamWrite(void * user, void * data, int size)
{
  auto bos = (BinaryOutputStream *)user;
  if (!bos) { return; }

  bos->writeBytes(size, data);
}

// Global set of callbacks, all access is read-only so it's threadsafe.
static stbi_io_callbacks const STB_CALLBACKS = { stbStreamRead, stbStreamSkip, stbStreamEof };

#endif // THEA_ENABLE_FREEIMAGE

} // namespace ImageInternal

int
Image::Type::numChannels() const
{
  switch (value)
  {
    case UNKNOWN   :  return -1;

    case LA_8U     :
    case LA_16U    :
    case LA_32F    :  return 2;

    case RGB_8U    :
    case RGB_16U   :
    case RGB_32F   :  return 3;

    case RGBA_8U   :
    case RGBA_16U  :
    case RGBA_32F  :  return 4;

    default        :  return 1;
  }
}

bool
Image::Type::isComplex() const
{
  return value == COMPLEX_64F;
}

bool
Image::Type::isFloatingPoint() const
{
  return value == LUMINANCE_32F
      || value == LUMINANCE_64F
      || value == LA_32F
      || value == RGB_32F
      || value == RGBA_32F
      || value == COMPLEX_64F;
}

int
Image::Type::getBitsPerPixel() const
{
  return ImageInternal::typeToBpp(*this);
}

int
Image::Type::getBitsPerChannel() const
{
  // Currently all supported types have the same number of bits per channel
  int nc = numChannels();
  return nc > 0 ? (getBitsPerPixel() / nc) : 0;
}

int
Image::Type::getBitsInChannel(int channel) const
{
  // Currently all supported types have the same number of bits per channel
  switch (numChannels())
  {
    case 1  : return channel == 0 ? getBitsPerPixel()   : 0;
    case 3  : return channel < 3  ? getBitsPerChannel() : 0;
    case 4  : return channel < 4  ? getBitsPerChannel() : 0;
    default : return 0;
  }
}

bool
Image::Type::hasByteAlignedPixels() const
{
  return getBitsPerPixel() >= 8;
}

bool
Image::Type::hasByteAlignedChannels() const
{
  return getBitsPerPixel() >= 8;
}

Image::Image()
: type(Type::UNKNOWN), width(0), height(0), depth(0), impl(nullptr), data(nullptr), owns_data(true), data_stride(0)
{
  cacheTypeProperties();
}

Image::Image(Type type_, int64 width_, int64 height_, int64 depth_)
: type(Type::UNKNOWN), width(0), height(0), depth(0), impl(nullptr), data(nullptr), owns_data(true), data_stride(0)
{
  cacheTypeProperties();
  resize(type_, width_, height_, depth_);
}

Image::Image(void * buf, Type type_, int64 width_, int64 height_, int64 depth_, int64 stride_bytes_)
: type(type_), width(width_), height(height_), depth(depth_), impl(nullptr), data(buf), owns_data(false),
  data_stride(stride_bytes_)
{
  alwaysAssertM(buf, "Image: Cannot wrap a null buffer");
  alwaysAssertM(type_ != Type::UNKNOWN, "Image: Cannot wrap buffer of unknown pixel type");
  alwaysAssertM(width_ >= 0 && height_ >= 0 && depth_ >= 0, "Image: Cannot wrap buffer with negative dimensions");

  cacheTypeProperties();

  if (stride_bytes_ <= 0)
    data_stride = ImageInternal::strideBytes(width_, type_.getBitsPerPixel(), /* row alignment = */ 1);
}

Image::Image(BinaryInputStream & input, Codec const & codec, bool read_block_header)
: type(Type::UNKNOWN), width(0), height(0), depth(0), impl(nullptr), data(nullptr), owns_data(true), data_stride(0)
{
  cacheTypeProperties();
  read(input, codec, read_block_header);
}

Image::Image(std::string const & path, Codec const & codec)
: type(Type::UNKNOWN), width(0), height(0), depth(0), impl(nullptr), data(nullptr), owns_data(true), data_stride(0)
{
  cacheTypeProperties();
  load(path, codec);
}

Image::Image(Image const & src)
: type(Type::UNKNOWN), width(0), height(0), depth(0), impl(nullptr), data(nullptr), owns_data(true), data_stride(0)
{
  cacheTypeProperties();
  *this = src;
}

Image &
Image::operator=(Image const & src)
{
  alwaysAssertM(!src.impl || !src.data, "Image: Source image has both externally and internally managed implementations");

  if (src.impl)
  {
#if THEA_ENABLE_FREEIMAGE
    if (!impl) { impl = new fipImage(*src.impl); }
    else       { *impl = *src.impl;              }
#else
    resizeImpl(src.type, src.width, src.height, src.depth);

    auto src_bytes = src.minTotalBytes();
    alwaysAssertM(src_bytes == minTotalBytes(), "Image: Sizes of source and destination image buffers don't match");

    std::memcpy(impl, src.impl, (size_t)src_bytes);
#endif
  }
  else
    clearImpl();

  if (src.data)
  {
    resizeData(src.type, src.width, src.height, src.depth);

    // Copy row by row since strides may not match
    auto row_bytes = ImageInternal::pixelBytes(src.type.getBitsPerPixel(), src.width);
    for (int64 i = 0; i < src.height; ++i)
    {
      auto src_scanline = src.getScanLine(i);
      auto dst_scanline = getScanLine(i);
      std::memcpy(dst_scanline, src_scanline, (size_t)row_bytes);
    }
  }
  else
    clearData();

  return *this;
}

Image::~Image()
{
  clear();
}

int8
Image::isValid() const
{
  theaAssertM(!impl || !data, "Image: Image has both externally and internally managed implementations");

  return type != Type::UNKNOWN && width >= 0 && height >= 0 && depth >= 0
      && ((width <= 0 || height <= 0 || depth <= 0)
#if THEA_ENABLE_FREEIMAGE
       || (impl && impl->isValid() != 0)
#else
       || impl
#endif
       || data);
}

void
Image::clearImpl()
{
  if (!impl) { return; }

#if THEA_ENABLE_FREEIMAGE
  delete impl;
#else
  stbi_image_free(impl);
#endif

  impl = nullptr;
  setAttribs(Type::UNKNOWN, 0, 0, 0, /* data_stride = */ 0);
}

void
Image::clearData()
{
  if (!data) { return; }

  if (owns_data) { data_allocator.deallocate((unsigned char *)data); }
  data = nullptr;
  owns_data = true;

  setAttribs(Type::UNKNOWN, 0, 0, 0, /* data_stride = */ 0);
}

int8
Image::clear()
{
  clearImpl();
  clearData();

  return true;
}

int64
Image::minTotalBytes() const
{
  return minTotalBytes(type, width, height, depth, getStrideBytes());
}

int64
Image::minTotalBytes(Type t, int64 w, int64 h, int64 d, int64 stride_bytes_, bool last_row_occupies_stride)
{
  if (t == Type::UNKNOWN || w <= 0 || h <= 0 || d <= 0) return 0;

  alwaysAssertM(stride_bytes_ > 0, "Image: Stride for non-empty image must be positive");

  auto bytes_pp = ImageInternal::pixelBytes(t.getBitsPerPixel());
  auto bytes = last_row_occupies_stride ? (int64)(stride_bytes_ * h * d)
                                        : (int64)(stride_bytes_ * (h * d - 1) + w * bytes_pp);

  theaAssertM(bytes >= 0, "Image: Total number of bytes computed as negative");

  return bytes;
}

int8
Image::resize(int64 type_, int64 width_, int64 height_, int64 depth_)
{
  Type t(type_);
  alwaysAssertM(type_ >= 0 && t != Type::UNKNOWN, "Image: Cannot resize to unknown type");

  // If an externally managed image which supports the desired type exists, resize it. Else if an internal data buffer exists,
  // resize it. If neither exists, resize the external one for 2D images and the internal one for 3D images.
  if (impl
   || (!data && depth_ == 1
#if !THEA_ENABLE_FREEIMAGE  // stb doesn't support all the formats in Image::Type
    && (t == Type::LUMINANCE_8U || t == Type::LUMINANCE_16U || t == Type::LUMINANCE_32F || t == Type::RGB_8U
     || t == Type::RGBA_8U || t == Type::RGB_16U || t == Type::RGBA_16U || t == Type::RGB_32F || t == Type::RGBA_32F)
#endif
      ))
    return resizeImpl(type_, width_, height_, depth_);
  else
    return resizeData(type_, width_, height_, depth_);
}

bool
Image::resizeImpl(int64 type_, int64 width_, int64 height_, int64 depth_)
{
  Type t(type_);
  alwaysAssertM(t != Type::UNKNOWN && width_ >= 0 && height_ >= 0 && depth_ >= 0,
                "Image: Cannot resize to unknown type or negative dimensions");

  // Data and impl cannot co-exist
  clearData();

#if THEA_ENABLE_FREEIMAGE

  if (width_ > 0 && height_ > 0 && depth_ > 0)
  {
    if (!impl)
    {
      impl = new fipImage;
      if (!impl) { THEA_ERROR << "Image: Could not create new FreeImage object"; return false; }
    }

    if (!impl->setSize(ImageInternal::typeToFreeImageType(t), (int)width_, (int)height_, ImageInternal::typeToBpp(t)))
    { THEA_ERROR << "Image: Could not allocate memory for resizing image"; return false; }
  }
  else if (impl)  // else just leave it as a null pointer
    impl->clear();

#else

  auto old_bytes = minTotalBytes();
  auto new_bytes = minTotalBytes(t, width_, height_, depth_,
                                 ImageInternal::strideBytes(width_, t.getBitsPerPixel(), /* stb row alignment = */ 1));
  if (new_bytes > old_bytes)  // like std::vector, don't deallocate if new size is smaller
  {
    impl = STBI_REALLOC(impl, new_bytes);
    if (!impl) { THEA_ERROR << "Image: Could not allocate memory for resizing image"; return false; }
  }

#endif

  setAttribs(t, width_, height_, depth_);

  if (!isValid())
  {
    THEA_ERROR << "Image: Could not resize to the specified type and dimensions " << width_ << 'x' << height_ << 'x' << depth_;
    return false;
  }

  return true;
}

bool
Image::resizeData(int64 type_, int64 width_, int64 height_, int64 depth_, int64 data_stride_)
{
  Type t(type_);
  alwaysAssertM(t != Type::UNKNOWN && width_ >= 0 && height_ >= 0 && depth_ >= 0,
                "Image: Cannot resize to unknown type or negative dimensions");

  // Impl and data cannot co-exist
  clearImpl();

  // Compute a default stride if the supplied one is non-positive
  if (data_stride_ <= 0)
    data_stride_ = ImageInternal::strideBytes(width_, t.getBitsPerPixel(), ROW_ALIGNMENT);

  auto old_bytes = minTotalBytes();
  auto new_bytes = minTotalBytes(t, width_, height_, depth_, data_stride_);
  if (new_bytes > old_bytes)  // like std::vector, don't deallocate if new size is smaller
  {
    if (!owns_data) { THEA_ERROR << "Image: Wrapped external buffer cannot be resized"; return false; }

    data_allocator.deallocate((unsigned char *)data);
    data = data_allocator.allocate(new_bytes);
    if (!data) { THEA_ERROR << "Image: Could not allocate memory for resizing image"; return false; }
    owns_data = true;
  }

  setAttribs(t, width_, height_, depth_, data_stride_);

  if (!isValid())
  {
    THEA_ERROR << "Image: Could not resize to the specified type and dimensions " << width_ << 'x' << height_ << 'x' << depth_;
    return false;
  }

  return true;
}

void const *
Image::getData() const
{
  return const_cast<Image *>(this)->getData();
}

void *
Image::getData()
{
  theaAssertM(!impl || !data, "Image: Image has both externally and internally managed implementations");

  if (impl)
  {
#if THEA_ENABLE_FREEIMAGE
    return impl->isValid() ? impl->accessPixels() : nullptr;
#else
    return impl;
#endif
  }
  else
    return data;
}

void const *
Image::getScanLine(int64 row, int64 z) const
{
  return const_cast<Image *>(this)->getScanLine(row, z);
}

void *
Image::getScanLine(int64 row, int64 z)
{
  theaAssertM(row >= 0 && row < height, "Image: Scanline row index out of bounds");
  theaAssertM(z   >= 0 && z   < depth,  "Image: Scanline Z index out of bounds");
  theaAssertM(isValid(),                "Image: Can't get scanline of an invalid image");

  if (impl)
  {
#if THEA_ENABLE_FREEIMAGE
    return impl->isValid() ? impl->getScanLine(row) : nullptr;
#else
    return (uint8 *)impl + (z * height + row) * getStrideBytes();
#endif
  }
  else
    return (uint8 *)data + (z * height + row) * getStrideBytes();
}

int64
Image::getStrideBytes() const
{
  if (impl)
  {
#if THEA_ENABLE_FREEIMAGE
    return impl->getScanWidth();
#else
    return ImageInternal::strideBytes(width, getBitsPerPixel(), 1);  // stb always has zero padding
#endif
  }
  else
    return data_stride;
}

double
Image::getNormalizedValue(void const * pixel, int channel) const
{
  if (channel < 0 || channel >= 4) { return 0; }

  switch (type)
  {
    case Image::Type::LUMINANCE_8U  : return channel == 0 ? *((uint8 const *)pixel) / 255.0 : 0.0;

    case Image::Type::LUMINANCE_16  : return channel == 0 ? *((int16  const *)pixel) / 32768.0 : 0.0;
    case Image::Type::LUMINANCE_16U : return channel == 0 ? *((uint16 const *)pixel) / 65535.0 : 0.0;

    case Image::Type::LUMINANCE_32  : return channel == 0 ? *((int32  const *)pixel) / 2147483648.0 : 0.0;
    case Image::Type::LUMINANCE_32U : return channel == 0 ? *((uint32 const *)pixel) / 4294967295.0 : 0.0;

    case Image::Type::LUMINANCE_32F : return channel == 0 ? *((float32 const *)pixel) : 0.0;
    case Image::Type::LUMINANCE_64F : return channel == 0 ? *((float64 const *)pixel) : 0.0;

    case Image::Type::LA_8U         : return channel < 2 ? ((uint8   const *)pixel)[channel] / 255.0 : 0.0;
    case Image::Type::LA_16U        : return channel < 2 ? ((uint16  const *)pixel)[channel] / 65535.0 : 0.0;
    case Image::Type::LA_32F        : return channel < 2 ? ((float32 const *)pixel)[channel] : 0.0;

    // Alpha of RGB images is assumed to be 1
    case Image::Type::RGB_8U        : return channel < 3 ? ((uint8   const *)pixel)[channel] / 255.0 : 1.0;
    case Image::Type::RGB_16U       : return channel < 3 ? ((uint16  const *)pixel)[channel] / 65535.0 : 1.0;
    case Image::Type::RGB_32F       : return channel < 3 ? ((float32 const *)pixel)[channel] : 1.0;

    case Image::Type::RGBA_8U       : return ((uint8   const *)pixel)[channel] / 255.0;
    case Image::Type::RGBA_16U      : return ((uint16  const *)pixel)[channel] / 65535.0;
    case Image::Type::RGBA_32F      : return ((float32 const *)pixel)[channel];

    case Image::Type::COMPLEX_64F   :
    {
      if (channel == 0)
      {
        float64 re = ((float64 const *)pixel)[0];
        float64 im = ((float64 const *)pixel)[1];
        return std::sqrt(re * re + im * im);
      }
      else
        return 0.0;
    }

    default                         : return 0.0;
  }
}

void
Image::setAttribs(Type type_, int64 w, int64 h, int64 d, int64 data_stride_)
{
  type = type_;
  cacheTypeProperties();

  if (w >= 0) { width  = w; }
  if (h >= 0) { height = h; }
  if (d >= 0) { depth  = d; }
  if (data_stride_ >= 0) { data_stride = data_stride_; }
}

void
Image::cacheTypeProperties()
{
  num_channels               =  type.numChannels();
  is_complex                 =  type.isComplex();
  is_floating_point          =  type.isFloatingPoint();
  bits_per_pixel             =  type.getBitsPerPixel();
  bits_per_channel           =  type.getBitsPerChannel();
  has_byte_aligned_pixels    =  type.hasByteAlignedPixels();
  has_byte_aligned_channels  =  type.hasByteAlignedChannels();
}

bool
Image::invert()
{
  if (impl)
  {
#if THEA_ENABLE_FREEIMAGE
    return !impl->isValid() || impl->invert() == TRUE;
#else
    THEA_ERROR << "Image: Cannot invert stb image";
    return false;
#endif
  }
  else
  {
    THEA_ERROR << "Image: Cannot invert internally managed pixel buffer";
    return false;
  }
}

bool
Image::convert(Type dst_type)
{
  return convert(dst_type, *this);
}

bool
Image::convert(Type dst_type, Image & dst) const
{
  bool status = false;
  if (type == dst_type)
  {
    if (&dst != this) { dst = *this; }
    status = true;
  }
  else if (!impl)
  {
    // TODO
    THEA_ERROR << "Image: Format conversion of internally managed buffers currently not supported";
  }
  else
  {
#if THEA_ENABLE_FREEIMAGE

    FREE_IMAGE_TYPE src_fitype = ImageInternal::typeToFreeImageType(type);
    FREE_IMAGE_TYPE dst_fitype = ImageInternal::typeToFreeImageType(dst_type);

    if (dst_fitype == FIT_BITMAP)
    {
      auto dst_fibpp = ImageInternal::typeToBpp(dst_type);
      switch (dst_fibpp)
      {
        case 4:
        {
          if (&dst != this) { dst = *this; }
          if (!ImageInternal::decanonicalizeChannelOrder(dst)) { return false; }

          if (src_fitype != FIT_BITMAP) { status = (dst.impl->convertToType(FIT_BITMAP) == TRUE); }
          if (status)                   { status = (dst.impl->convertTo4Bits() == TRUE);          }

          break;
        }

        case 8:
        {
          if (&dst != this) { dst = *this; }
          if (!ImageInternal::decanonicalizeChannelOrder(dst)) { return false; }

          if (src_fitype != FIT_BITMAP) { status = (dst.impl->convertToType(FIT_BITMAP) == TRUE); }
          else                          { status = true;                                          }

          if (status)
            status = (dst.impl->convertToGrayscale() == TRUE);  // convertTo8Bits() can palettize

          break;
        }

        case 24:
        {
          if (&dst != this) { dst = *this; }
          if (!ImageInternal::decanonicalizeChannelOrder(dst)) { return false; }

          status = (dst.impl->convertTo24Bits() == TRUE);
          break;
        }

        case 32:
        {
          if (&dst != this) { dst = *this; }
          if (!ImageInternal::decanonicalizeChannelOrder(dst)) { return false; }

          status = (dst.impl->convertTo32Bits() == TRUE);
          break;
        }
      }
    }
    else
    {
      if (&dst != this) { dst = *this; }
      if (!ImageInternal::decanonicalizeChannelOrder(dst)) { return false; }

      if (dst.impl->convertToType(dst_fitype) == TRUE)
        status = true;
    }

    if (status)
    {
      dst.setAttribs(dst_type);

      if (!ImageInternal::canonicalizeChannelOrder(dst)) { return false; }
    }

#else // THEA_ENABLE_FREEIMAGE

    // TODO
    THEA_ERROR << "Image: Format conversion of STB images currently not supported";

#endif // THEA_ENABLE_FREEIMAGE
  }

  return status;
}

namespace ImageInternal {

template <int Bits> struct UInt {};
template <> struct UInt<8>  { typedef uint8  type; };
template <> struct UInt<16> { typedef uint16 type; };
template <> struct UInt<32> { typedef uint32 type; };
template <> struct UInt<64> { typedef uint64 type; };

template <int BitsPerChannel>
bool
reorderChannels(Image & img, int const * order)
{
  alwaysAssertM(img.getBitsPerChannel() == BitsPerChannel, "Image: Bits-per-channel does not match template parameter");

  typedef typename UInt<BitsPerChannel>::type ChannelT;

  size_t const num_channels = static_cast<size_t>(img.numChannels());
  for (size_t i = 0; i < num_channels; ++i)
    if (order[i] < 0 || order[i] >= num_channels)
    {
      THEA_ERROR << "Image: Invalid channel index: order[" << i << "] = " << order[i];
      return false;
    }

  auto depth = img.getDepth(), height = img.getHeight(), width = img.getWidth();
  ChannelT reordered[256];  // one less thing to change if we support images with more than 4 channels later
  for (int64 i = 0; i < depth; ++i)
    for (int64 j = 0; j < height; ++j)
    {
      auto pixel = static_cast<ChannelT *>(img.getScanLine(j, i));
      for (int64 k = 0; k < width; ++k, pixel += num_channels)
      {
        for (size_t c = 0; c < num_channels; ++c)
          reordered[c] = pixel[order[c]];

        std::memcpy(pixel, reordered, sizeof(ChannelT) * num_channels);
      }
    }

  return true;
}

} // namespace ImageInternal

bool
Image::reorderChannels(int const * order)
{
  if (numChannels() <= 1)
  {
    if (numChannels() == 1 && order[0] != 0)
    {
      THEA_ERROR << "Image: Invalid channel index: order[0] = " << order[0];
      return false;
    }

    return true;
  }

  switch (getBitsPerChannel())
  {
    case   8: return ImageInternal::reorderChannels<  8>(*this, order);
    case  16: return ImageInternal::reorderChannels< 16>(*this, order);
    case  32: return ImageInternal::reorderChannels< 32>(*this, order);
    case  64: return ImageInternal::reorderChannels< 64>(*this, order);
    // No need to reorder 128-bit channels (Type::COMPLEX_64F) since this is a single-channel format
    default:
    {
      THEA_ERROR << "Image: Channel reordering not supported for " << getBitsPerChannel() << " bits per channel";
      return false;
    }
  }
}

namespace ImageInternal {

#if THEA_ENABLE_FREEIMAGE

FREE_IMAGE_FILTER
filterToImplFilter(Image::Filter filter)
{
  switch (filter)
  {
    case Image::Filter::BOX          :  return FILTER_BOX; break;
    case Image::Filter::BILINEAR     :  return FILTER_BILINEAR; break;
    case Image::Filter::BSPLINE      :  return FILTER_BSPLINE; break;
    case Image::Filter::BICUBIC      :  return FILTER_BICUBIC; break;
    case Image::Filter::CATMULL_ROM  :  return FILTER_CATMULLROM; break;
    case Image::Filter::LANCZOS3     :  return FILTER_LANCZOS3; break;
    default                          :  return FILTER_BICUBIC;
  }
}

#else // THEA_ENABLE_FREEIMAGE

stbir_filter
filterToImplFilter(Image::Filter filter)
{
  switch (filter)
  {
    case Image::Filter::AUTO:         return STBIR_FILTER_DEFAULT;
    case Image::Filter::BOX:          return STBIR_FILTER_BOX;
    case Image::Filter::BILINEAR:     return STBIR_FILTER_TRIANGLE;
    case Image::Filter::BSPLINE:      return STBIR_FILTER_CUBICBSPLINE;
    case Image::Filter::BICUBIC:      return STBIR_FILTER_MITCHELL;
    case Image::Filter::CATMULL_ROM:  return STBIR_FILTER_CATMULLROM;
    default: throw Error("Image: Unsupported rescaling filter");
  }
}

#endif // THEA_ENABLE_FREEIMAGE

} // namespace ImageInternal

bool
Image::rescale(int64 new_width, int64 new_height, int64 new_depth, Filter filter)
{
  if (depth != 1 || new_depth != 1)
  {
    THEA_ERROR << "Image: Rescaling non-2D images currently not supported";
    return false;
  }

  if (!isValid())
  {
    THEA_ERROR << "Image: Attempting to rescale an invalid image";
    return false;
  }

  if (width <= 0 || height <= 0 || new_width <= 0 || new_height <= 0)
  {
    THEA_ERROR << "Image: Attempting to rescale between invalid dimensions: " << width << " x " << height << " to "
                                                                              << new_width << " x " << new_height;
    return false;
  }

  if (!impl)
  {
    THEA_ERROR << "Image: Rescaling wrapped buffers currently not supported";  // FIXME: Can use stb
    return false;
  }

  auto impl_filter = ImageInternal::filterToImplFilter(filter);
  bool status = false;

#if THEA_ENABLE_FREEIMAGE

  status = (impl->rescale((unsigned int)new_width, (unsigned int)new_height, impl_filter) == TRUE);

#else

  auto new_stride = ImageInternal::strideBytes(new_width, getBitsPerPixel(), /* stb row alignment = */ 1);
  size_t rescaled_bytes = (size_t)minTotalBytes(type, new_width, new_height, new_depth, new_stride);
  void * new_impl = STBI_MALLOC(rescaled_bytes);
  if (!new_impl)
  {
    THEA_ERROR << "Image: Could not allocate memory for rescaled image";
    return false;
  }

  auto bpc = getBitsPerChannel();
  if (bpc == 8)
  {
    status = (bool)stbir_resize_uint8_generic((uint8 const *)impl, (int)width, (int)height, (int)getStrideBytes(),
                                              (uint8 *)new_impl, (int)new_width, (int)new_height, (int)new_stride,
                                              numChannels(), /* alpha_channel = */ numChannels() - 1, /* flags = */ 0,
                                              STBIR_EDGE_CLAMP, impl_filter, STBIR_COLORSPACE_LINEAR,
                                              /* alloc_context = */ nullptr);
  }
  else if (bpc == 16)
  {
    status = (bool)stbir_resize_uint16_generic((uint16 const *)impl, (int)width, (int)height, (int)getStrideBytes(),
                                               (uint16 *)new_impl, (int)new_width, (int)new_height, (int)new_stride,
                                               numChannels(), /* alpha_channel = */ numChannels() - 1, /* flags = */ 0,
                                               STBIR_EDGE_CLAMP, impl_filter, STBIR_COLORSPACE_LINEAR,
                                               /* alloc_context = */ nullptr);
  }
  else if (bpc == 32 && isFloatingPoint())
  {
    status = (bool)stbir_resize_float_generic((float const *)impl, (int)width, (int)height, (int)getStrideBytes(),
                                              (float *)new_impl, (int)new_width, (int)new_height, (int)new_stride,
                                              numChannels(), /* alpha_channel = */ numChannels() - 1, /* flags = */ 0,
                                              STBIR_EDGE_CLAMP, impl_filter, STBIR_COLORSPACE_LINEAR,
                                              /* alloc_context = */ nullptr);
  }
  else
  {
    THEA_ERROR << "Image: Unsupported format for rescaling (must be 8/16-bit integer or 32-bit float per channel)";
    stbi_image_free(new_impl);
    return false;
  }

#endif

  if (!status)
  {
    THEA_ERROR << "Image: Could not rescale " << width << " x " << height << " image to " << new_width << " x " << new_height;

#if !THEA_ENABLE_FREEIMAGE
    stbi_image_free(new_impl);
#endif

    return false;
  }

#if !THEA_ENABLE_FREEIMAGE
  stbi_image_free(impl);
  impl = new_impl;
#endif

  width = new_width;
  height = new_height;
  depth = new_depth;

  return true;
}

void
Image::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  // Codecs with custom read functions
  if (codec == Codec3bm())
  {
    read3bm(codec, input, read_block_header);
    return;
  }

#if THEA_ENABLE_FREEIMAGE

  // Get the size of the image block in bytes
  auto len = (read_block_header ? (int64)Codec::BlockHeader(input).data_size : input.size());

  // Read the image block into a memory buffer (optimization possible when the data has already been buffered within the input
  // stream?)
  Array<uint8> img_block((size_t)len);
  input.readBytes(len, img_block.data());

  // Decode the image
  clearData();
  if (!impl)
  {
    impl = new fipImage;
    if (!impl) { throw Error("Image: Could not create FreeImage object"); }
  }

  fipMemoryIO mem((BYTE *)img_block.data(), (DWORD)len);
  if (codec == CodecAuto())
  {
    if (impl->loadFromMemory(mem) != TRUE)
      throw Error("Image: Could not decode image from memory stream");
  }
  else
  {
    ImageCodec const * img_codec = dynamic_cast<ImageCodec const *>(&codec);
    if (!img_codec) { Error("Image: Codec specified for image reading is not an image codec"); }

    if (impl->loadFromMemory((FREE_IMAGE_FORMAT)img_codec->getImplFormat(), mem, img_codec->getImplReadFlags()) != TRUE)
      throw Error("Image: Could not deserialize image from memory stream");
  }

  auto t = ImageInternal::typeFromFreeImageTypeAndBpp(impl->getImageType(), impl->getBitsPerPixel());
  if (t == Type::UNKNOWN)
  {
    clear();
    throw Error("Image: Image successfully deserialized but has unsupported format");
  }

  setAttribs(t, impl->getWidth(), impl->getHeight(), /* d = */ 1);

  if (!ImageInternal::canonicalizeChannelOrder(*this))
    throw Error("Image: Could not canonicalize channel order");

#else // THEA_ENABLE_FREEIMAGE

  // Decode via stb (codec is always auto-detected regardless of what was passed to this function)
  clear();

  // If a header is present, use it to create a new stream that wraps just the relevant part of the main input stream, to
  // prevent reading past the end of the buffer
  BinaryInputStream::Ptr tmp_in;
  BinaryInputStream * img_in = nullptr;
  if (read_block_header)
  {
    auto len = (int64)Codec::BlockHeader(input).data_size;
    tmp_in = std::make_shared<BinaryInputStream>(input, len);
    img_in = tmp_in.get();
  }
  else
    img_in = &input;

  stbi_set_flip_vertically_on_load(true);  // origin at bottom-left

  auto bpc = getBitsPerChannel();  // try to respect existing bit depth setting, fall back on 8 bits per channel
  int w, h, c;
  switch (bpc)
  {
    case 16: impl = stbi_load_16_from_callbacks(&ImageInternal::STB_CALLBACKS, img_in, &w, &h, &c, 0); break;
    case 32: impl = stbi_loadf_from_callbacks  (&ImageInternal::STB_CALLBACKS, img_in, &w, &h, &c, 0); break;
    default: impl = stbi_load_from_callbacks   (&ImageInternal::STB_CALLBACKS, img_in, &w, &h, &c, 0);
  }

  if (!impl)
    throw Error("Image: Could not read image");

  auto t = ImageInternal::typeFromStbChannelsAndBpc(c, bpc);
  if (t == Type::UNKNOWN)
  {
    clear();
    throw Error("Image: Deserialized image has unsupported type");
  }

  setAttribs(t, w, h, /* d = */ 1);

  if (!isValid())  // check dimensions
    throw Error("Image: Deserialized image has invalid dimensions");

#endif // THEA_ENABLE_FREEIMAGE
}

void
Image::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  // Codecs with custom write functions
  if (codec == Codec3bm())
  {
    write3bm(codec, output, write_block_header);
    return;
  }

  if (!impl || !isValid())
    throw Error("Image: Can't write an invalid or empty image");

  if (codec == CodecAuto())
    throw Error("Image: You must explicitly choose a codec for writing images");

#if THEA_ENABLE_FREEIMAGE

  ImageCodec const * img_codec = dynamic_cast<ImageCodec const *>(&codec);
  if (!img_codec)
    throw Error("Image: Codec specified for image writing is not an image codec");

  if (!ImageInternal::decanonicalizeChannelOrder(const_cast<Image &>(*this)))
    throw Error("Image: Could not decanonicalize channel order");

  fipMemoryIO mem;
  if (const_cast<Impl *>(impl)->saveToMemory((FREE_IMAGE_FORMAT)img_codec->getImplFormat(), mem,
                                             img_codec->getImplWriteFlags()) != TRUE)
    throw Error(toString(codec.getName()) + ": Could not save image to memory stream");

  if (!ImageInternal::canonicalizeChannelOrder(const_cast<Image &>(*this)))
    throw Error("Image: Could not canonicalize channel order");

  BYTE * buf;
  DWORD size_in_bytes;
  if (!mem.acquire(&buf, &size_in_bytes))
    throw Error(toString(codec.getName()) + ": Error accessing the FreeImage memory stream");

  if (write_block_header)
  {
    Codec::BlockHeader header(codec.getMagic(), (uint64)size_in_bytes);
    header.write(output);
  }

  output.writeBytes(size_in_bytes, buf);

#else // THEA_ENABLE_FREEIMAGE

  stbi_flip_vertically_on_write(true);  // had flipped when loading

  int nc = (int)numChannels();
  int status = 0;
  if (codec == CodecBmp())
  {
    if (getBitsPerChannel() != 8) { throw Error(toString(codec.getName()) + ": Writing BMP requires 8-bit channels"); }

    status = stbi_write_bmp_to_func(ImageInternal::stbStreamWrite, &output, (int)width, (int)height, nc, impl);
  }
  else if (codec == CodecHdr())
  {
    if (getBitsPerChannel() != 32 || !isFloatingPoint())
    { throw Error(toString(codec.getName()) + ": Writing HDR requires 32-bit float channels"); }

    status = stbi_write_hdr_to_func(ImageInternal::stbStreamWrite, &output, (int)width, (int)height, nc, (float const *)impl);
  }
  else if (codec == CodecJpg())
  {
    if (getBitsPerChannel() != 8) { throw Error(toString(codec.getName()) + ": Writing JPEG requires 8-bit channels"); }

    int q = (int)dynamic_cast<CodecJpg const &>(codec).getOptions().quality;
    if (q < 0) { q = 90; }
    else { alwaysAssertM(q >= 1 && q <= 100, toString(codec.getName()) + ": JPEG quality should be in the range [1, 100]"); }

    status = stbi_write_jpg_to_func(ImageInternal::stbStreamWrite, &output, (int)width, (int)height, nc, impl, q);
  }
  else if (codec == CodecPng())
  {
    if (getBitsPerChannel() != 8) { throw Error(toString(codec.getName()) + ": Writing PNG requires 8-bit channels"); }

    status = stbi_write_png_to_func(ImageInternal::stbStreamWrite, &output, (int)width, (int)height, nc, impl,
                                    (int)getStrideBytes());
  }
  else if (codec == CodecTga())
  {
    if (getBitsPerChannel() != 8) { throw Error(toString(codec.getName()) + ": Writing TGA requires 8-bit channels"); }

    status = stbi_write_tga_to_func(ImageInternal::stbStreamWrite, &output, (int)width, (int)height, nc, impl);
  }
  else
    throw Error("Image: Codec '" + toString(codec.getName()) + "' not supported for writing images");

  if (!status)
    throw Error(toString(codec.getName()) + ": Could not write image");

#endif // THEA_ENABLE_FREEIMAGE
}

void
Image::load(std::string const & path, Codec const & codec)
{
  BinaryInputStream in(path, Endianness::LITTLE);
  int64 file_size = in.size();
  if (file_size <= 0)
    throw Error("Image: File does not exist or is empty");

  read(in, codec, false);
}

void
Image::save(std::string const & path, Codec const & codec) const
{
  if (!isValid())
    throw Error("Image: Can't save an invalid image");

  Codec const * c = nullptr;
  if (codec == CodecAuto())
  {
    c = ImageInternal::codecFromPath(path);
    if (!c)
      throw Error("Image: Could not detect codec from path");
  }
  else
    c = &codec;

  BinaryOutputStream out(path, Endianness::LITTLE);
  if (!out.ok())
    throw Error("Image: Could not open image file for writing");

  write(out, *c, false);
}

void
Image::read3bm(Codec const & codec, BinaryInputStream & input, bool read_block_header)
{
  // Get the size of the image block in bytes
  uint64 prefixed_len = (read_block_header ? Codec::BlockHeader(input).data_size : 0);

  BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

  int64 start_position = input.getPosition();

  // Read file header (32 bytes):
  //   4 bytes: Magic string "3BM\0"
  //   4 bytes: Version
  //   8 bytes: Size of the whole file in bytes
  //   8 bytes: Reserved
  //   8 bytes: Starting offset (from the beginning of the file) of the bitmap image data

  uint8 magic[4];
  input.readBytes(4, &magic);
  if (magic[0] != (uint8)'3' || magic[1] != (uint8)'B' || magic[2] != (uint8)'M' || magic[3] != (uint8)'\0')
    throw Error(toString(codec.getName()) + ": Image is not in 3BM format");

  input.skip(4);
  uint64 len = input.readUInt64();
  if (prefixed_len > 0 && len != prefixed_len)
    throw Error(toString(codec.getName()) + ": Prefixed size does not match image size from header");

  input.skip(8);
  uint64 bitmap_offset = input.readUInt64();

  // Read info header (72 bytes):
  //   8 bytes: Size of this header (72 bytes)
  //   8 bytes: Width of the bitmap in voxels
  //   8 bytes: Height of the bitmap in voxels
  //   8 bytes: Depth of the bitmap in voxels
  //   4 bytes: Number of bits per pixel (1, 4, 8, 16, 24, 32), currently only 8, 16, 24 and 32 are supported
  //   4 bytes: Compression method (0 for no compression, currently this is the only one supported)
  //   8 bytes: Size of raw image data
  //   8 bytes: Width resolution in voxels/meter
  //   8 bytes: Height resolution in voxels/meter
  //   8 bytes: Depth resolution in voxels/meter

  input.skip(8);
  int64 w  =  (int64)input.readUInt64();
  int64 h  =  (int64)input.readUInt64();
  int64 d  =  (int64)input.readUInt64();

  int bpp = (int)input.readUInt32();
  if (bpp != 8 && bpp != 16 && bpp != 24 && bpp != 32)
    throw Error(toString(codec.getName()) + ": Only 8, 16, 24 and 32-bit bitmaps currently supported");

  uint32 compression = input.readUInt32();
  if (compression != 0)
    throw Error(toString(codec.getName()) + ": Unsupported compression method");

  int64 stream_stride = ImageInternal::strideBytes(w, bpp, 16);  // 16-byte row alignment for SSE compatibility
  uint64 data_len = input.readUInt64();
  if (data_len != (uint64)stream_stride * (uint64)h * (uint64)d)
    throw Error(toString(codec.getName()) + ": Pixel data block size does not match image dimensions");

  // Read image data as L, BGR or BGRA voxels, one 16-byte aligned scanline at a time
  input.setPosition(start_position + (int64)bitmap_offset);

  Type t = Type::UNKNOWN;
  switch (bpp)
  {
    case 8   :  t = Type::LUMINANCE_8U; break;
    case 16  :  t = Type::LA_8U; break;
    case 24  :  t = Type::RGB_8U; break;
    case 32  :  t = Type::RGBA_8U; break;
    default  :  throw Error(format("%s: Bit d %d not supported", codec.getName(), bpp));
  }

  resizeData(t, w, h, d);
  if (w > 0 && h > 0 && d > 0 && bpp > 0)
  {
    int bytes_pp = bpp / 8;
    Array<uint8> row_buf((size_t)stream_stride);
    for (int64 i = 0; i < d; ++i)
      for (int64 j = 0; j < h; ++j)
      {
        input.readBytes(stream_stride, row_buf.data());

        uint8 const * in_pixel = row_buf.data();
        uint8 * out_pixel = (uint8 *)getScanLine(j, i);

        switch (bytes_pp)
        {
          case 1:  // no channel..,
          case 2:  // ... reordering
          {
            std::copy(in_pixel, in_pixel + w * bytes_pp, out_pixel);
            break;
          }
          case 3:
          {
            for (int64 k = 0; k < w; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
            {
              out_pixel[0] = in_pixel[2];
              out_pixel[1] = in_pixel[1];
              out_pixel[2] = in_pixel[0];
            }
            break;
          }
          default: // 4
          {
            for (int64 k = 0; k < w; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
            {
              out_pixel[0] = in_pixel[2];
              out_pixel[1] = in_pixel[1];
              out_pixel[2] = in_pixel[0];
              out_pixel[3] = in_pixel[3];
            }
            break;
          }
        }
      }
  }
}

void
Image::write3bm(Codec const & codec, BinaryOutputStream & output, bool write_block_header) const
{
  if (!isValid())
    throw Error(toString(codec.getName()) + ": Cannot write an invalid image");

  if (type != Type::LUMINANCE_8U
   && type != Type::LA_8U
   && type != Type::RGB_8U
   && type != Type::RGBA_8U)
    throw Error(toString(codec.getName()) + ": Can only write 8-bit luminance, RGB or RGBA images");

  auto bpp = getBitsPerPixel();
  int64 stream_stride = ImageInternal::strideBytes(width, bpp, 16);  // 16-byte row alignment for SSE compatibility

  static uint64 const FILE_HEADER_SIZE = 32;  // Remember to change these if the format changes!
  static uint64 const INFO_HEADER_SIZE = 72;
  uint64 data_len = (uint64)stream_stride * (uint64)height * (uint64)depth;
  uint64 size_in_bytes = FILE_HEADER_SIZE + INFO_HEADER_SIZE + data_len;

  if (write_block_header)
    Codec::BlockHeader(codec.getMagic(), size_in_bytes).write(output);

  BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

  // Write file header (32 bytes):
  //   4 bytes: Magic string "3BM\0"
  //   4 bytes: Version
  //   8 bytes: Size of the whole file in bytes
  //   8 bytes: Reserved
  //   8 bytes: Starting offset (from the beginning of the file) of the bitmap image data

  static uint8 const MAGIC[4] = { (uint8)'3', (uint8)'B', (uint8)'M', (uint8)'\0' };
  output.writeBytes(4, MAGIC);
  output.writeUInt32(0x0001);
  output.writeUInt64(size_in_bytes);
  output.writeUInt64(0);
  output.writeUInt64(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

  // Write info header (72 bytes):
  //   8 bytes: Size of this header (72 bytes)
  //   8 bytes: Width of the bitmap in voxels
  //   8 bytes: Height of the bitmap in voxels
  //   8 bytes: Depth of the bitmap in voxels
  //   4 bytes: Number of bits per pixel (1, 4, 8, 16, 24, 32), currently only 8, 16, 24 and 32 are supported
  //   4 bytes: Compression method (0 for no compression, currently this is the only one supported)
  //   8 bytes: Size of raw image data
  //   8 bytes: Width resolution in voxels/meter
  //   8 bytes: Height resolution in voxels/meter
  //   8 bytes: Depth resolution in voxels/meter

  output.writeUInt64(INFO_HEADER_SIZE);
  output.writeUInt64((uint64)width);
  output.writeUInt64((uint64)height);
  output.writeUInt64((uint64)depth);
  output.writeUInt32((uint32)bpp);
  output.writeUInt32(0);
  output.writeUInt64(data_len);
  output.writeUInt64(3937);  // 100dpi
  output.writeUInt64(3937);  // 100dpi
  output.writeUInt64(3937);  // 100dpi

  // Write image data as L, BGR or BGRA voxels, one 16-byte aligned scanline at a time
  if (width <= 0 || height <= 0 || depth <= 0 || bpp <= 0)
    return;

  int bytes_pp = bpp / 8;
  Array<uint8> row_buf((bytes_pp > 2 ? (size_t)stream_stride : 0), 0);
  auto row_padding = stream_stride - (width * bytes_pp);
  for (int64 i = 0; i < depth; ++i)
    for (int64 j = 0; j < height; ++j)
    {
      uint8 const * in_pixel = (uint8 const *)getScanLine(j, i);

      switch (bytes_pp)
      {
        case 1:  // no channel...
        case 2:  // ... reordering
        {
          output.writeBytes(width * bytes_pp, in_pixel);
          output.skip(row_padding);  // writes undefined values
          break;
        }
        case 3:
        {
          uint8 * out_pixel = row_buf.data();
          for (int64 k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
          {
            out_pixel[0] = in_pixel[2];
            out_pixel[1] = in_pixel[1];
            out_pixel[2] = in_pixel[0];
          }
          output.writeBytes((int64)stream_stride, row_buf.data());
          break;
        }
        default: // 4
        {
          uint8 * out_pixel = row_buf.data();
          for (int64 k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
          {
            out_pixel[0] = in_pixel[2];
            out_pixel[1] = in_pixel[1];
            out_pixel[2] = in_pixel[0];
            out_pixel[3] = in_pixel[3];
          }
          output.writeBytes((int64)stream_stride, row_buf.data());
          break;
        }
      }
    }
}

} // namespace Thea
