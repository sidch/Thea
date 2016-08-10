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
#include <FreeImagePlus.h>
#include <cmath>

THEA_INSTANTIATE_SMART_POINTERS(Thea::Image)

namespace Thea {

int const AbstractImage::Channel::RED    =  FI_RGBA_RED;
int const AbstractImage::Channel::GREEN  =  FI_RGBA_GREEN;
int const AbstractImage::Channel::BLUE   =  FI_RGBA_BLUE;
int const AbstractImage::Channel::ALPHA  =  FI_RGBA_ALPHA;

namespace ImageInternal {

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

static WORD
typeToFreeImageBPP(Image::Type type)
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

static Image::Type
typeFromFreeImageTypeAndBPP(FREE_IMAGE_TYPE fi_type, WORD fi_bpp)
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

static ImageCodec const *
codecFromPath(std::string const & path)
{
  typedef shared_ptr<ImageCodec> ICPtr;
  typedef TheaUnorderedMap<std::string, ICPtr > CodecMap;

  static CodecMap codec_map;
  if (codec_map.empty())  // only on the first call
  {
    codec_map["BMP"  ] = ICPtr(new CodecBMP());
    codec_map["CUT"  ] = ICPtr(new CodecCUT());
    codec_map["DDS"  ] = ICPtr(new CodecDDS());
    codec_map["EXR"  ] = ICPtr(new CodecEXR());
    codec_map["GIF"  ] = ICPtr(new CodecGIF());
    codec_map["HDR"  ] = ICPtr(new CodecHDR());
    codec_map["ICO"  ] = ICPtr(new CodecICO());
    codec_map["IFF"  ] = ICPtr(new CodecIFF());
    codec_map["LBM"  ] = ICPtr(new CodecIFF());
    codec_map["J2K"  ] = ICPtr(new CodecJ2K());
    codec_map["J2C"  ] = ICPtr(new CodecJ2K());
    codec_map["JNG"  ] = ICPtr(new CodecJNG());
    codec_map["JP2"  ] = ICPtr(new CodecJP2());
    codec_map["JPG"  ] = ICPtr(new CodecJPEG());
    codec_map["JIF"  ] = ICPtr(new CodecJPEG());
    codec_map["JPEG" ] = ICPtr(new CodecJPEG());
    codec_map["JPE"  ] = ICPtr(new CodecJPEG());
    codec_map["KOA"  ] = ICPtr(new CodecKOALA());
    codec_map["MNG"  ] = ICPtr(new CodecMNG());
    codec_map["PBM"  ] = ICPtr(new CodecPBM());
    codec_map["PBM"  ] = ICPtr(new CodecPBMRAW());
    codec_map["PCD"  ] = ICPtr(new CodecPCD());
    codec_map["PCX"  ] = ICPtr(new CodecPCX());
    codec_map["PFM"  ] = ICPtr(new CodecPFM());
    codec_map["PGM"  ] = ICPtr(new CodecPGM());
    codec_map["PGM"  ] = ICPtr(new CodecPGMRAW());
    codec_map["PNG"  ] = ICPtr(new CodecPNG());
    codec_map["PPM"  ] = ICPtr(new CodecPPM());
    codec_map["PPM"  ] = ICPtr(new CodecPPMRAW());
    codec_map["PSD"  ] = ICPtr(new CodecPSD());
    codec_map["RAS"  ] = ICPtr(new CodecRAS());
    codec_map["SGI"  ] = ICPtr(new CodecSGI());
    codec_map["TGA"  ] = ICPtr(new CodecTARGA());
    codec_map["TARGA"] = ICPtr(new CodecTARGA());
    codec_map["TIF"  ] = ICPtr(new CodecTIFF());
    codec_map["TIFF" ] = ICPtr(new CodecTIFF());
    codec_map["WBMP" ] = ICPtr(new CodecWBMP());
    codec_map["XBM"  ] = ICPtr(new CodecXBM());
    codec_map["XPM"  ] = ICPtr(new CodecXPM());
  }

  CodecMap::const_iterator existing = codec_map.find(toUpper(FilePath::extension(path)));
  if (existing == codec_map.end())
    return NULL;
  else
    return existing->second.get();
}

} // namespace ImageInternal

#define THEA_DEF_SERIALIZE_IMAGE(codec, fip_format, flags)                                                                    \
long                                                                                                                          \
codec::serializeImage(Image const & image, BinaryOutputStream & output, bool prefix_info) const                               \
{                                                                                                                             \
  fipMemoryIO mem;                                                                                                            \
  image._getFreeImage()->saveToMemory(fip_format, mem, flags);                                                                \
                                                                                                                              \
  BYTE * data;                                                                                                                \
  DWORD size_in_bytes;                                                                                                        \
  if (!mem.acquire(&data, &size_in_bytes))                                                                                    \
    throw Error(std::string(getName()) + ": Error accessing the FreeImage memory stream");                                    \
                                                                                                                              \
  if (prefix_info)                                                                                                            \
  {                                                                                                                           \
    output.setEndianness(Endianness::LITTLE);                                                                                 \
    output.writeUInt64(static_cast<uint64>(size_in_bytes));                                                                   \
  }                                                                                                                           \
                                                                                                                              \
  output.writeBytes(size_in_bytes, data);                                                                                     \
                                                                                                                              \
  return (long)(prefix_info ? size_in_bytes + 4 : size_in_bytes);                                                             \
}

// TODO: Add options to all the ones that support them
THEA_DEF_SERIALIZE_IMAGE(CodecBMP,     FIF_BMP,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecCUT,     FIF_CUT,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecDDS,     FIF_DDS,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecEXR,     FIF_EXR,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecGIF,     FIF_GIF,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecHDR,     FIF_HDR,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecICO,     FIF_ICO,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecIFF,     FIF_IFF,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecJ2K,     FIF_J2K,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecJNG,     FIF_JNG,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecJP2,     FIF_JP2,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecJPEG,    FIF_JPEG,    options.quality | (options.progressive ? JPEG_PROGRESSIVE : 0))
THEA_DEF_SERIALIZE_IMAGE(CodecKOALA,   FIF_KOALA,   0)
THEA_DEF_SERIALIZE_IMAGE(CodecMNG,     FIF_MNG,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPBM,     FIF_PBM,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPBMRAW,  FIF_PBMRAW,  0)
THEA_DEF_SERIALIZE_IMAGE(CodecPCD,     FIF_PCD,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPCX,     FIF_PCX,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPFM,     FIF_PFM,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPGM,     FIF_PGM,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPGMRAW,  FIF_PGMRAW,  0)
THEA_DEF_SERIALIZE_IMAGE(CodecPNG,     FIF_PNG,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPPM,     FIF_PPM,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecPPMRAW,  FIF_PPMRAW,  0)
THEA_DEF_SERIALIZE_IMAGE(CodecPSD,     FIF_PSD,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecRAS,     FIF_RAS,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecSGI,     FIF_SGI,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecTARGA,   FIF_TARGA,   0)
THEA_DEF_SERIALIZE_IMAGE(CodecTIFF,    FIF_TIFF,    0)
THEA_DEF_SERIALIZE_IMAGE(CodecWBMP,    FIF_WBMP,    0)
THEA_DEF_SERIALIZE_IMAGE(CodecXBM,     FIF_XBM,     0)
THEA_DEF_SERIALIZE_IMAGE(CodecXPM,     FIF_XPM,     0)

#define THEA_DEF_DESERIALIZE_IMAGE(codec, fip_format, flags)                                                                  \
void                                                                                                                          \
codec::deserializeImage(Image & image, BinaryInputStream & input, bool read_prefixed_info) const                              \
{                                                                                                                             \
  /* Get the size of the image block in bytes */                                                                              \
  input.setEndianness(Endianness::LITTLE);                                                                                    \
  uint64 size = read_prefixed_info ? input.readUInt64() : input.size();                                                       \
                                                                                                                              \
  /* Read the image block into a memory buffer (optimization possible when the data has already been buffered within the */   \
  /* input stream?)  */                                                                                                       \
  TheaArray<uint8> img_block((array_size_t)size);                                                                             \
  input.readBytes((int64)size, &img_block[0]);                                                                                \
                                                                                                                              \
  /* Decode the image */                                                                                                      \
  fipMemoryIO mem((BYTE *)&img_block[0], (DWORD)size);                                                                        \
  FIBITMAP * bitmap = mem.load(fip_format, flags);                                                                            \
  if (!bitmap)                                                                                                                \
    throw Error(std::string(getName()) + ": Could not decode image from memory stream");                                      \
                                                                                                                              \
  fipImage * fip_img = image._getFreeImage();                                                                                 \
  debugAssertM(fip_img, std::string(getName()) + ": Image does not wrap a valid FreeImage bitmap");                           \
                                                                                                                              \
  *fip_img = bitmap;  /* the FIP object will now manage the destruction of the bitmap */                                      \
                                                                                                                              \
  Image::Type type = ImageInternal::typeFromFreeImageTypeAndBPP(fip_img->getImageType(), fip_img->getBitsPerPixel());         \
  if (type == Image::Type::UNKNOWN)                                                                                           \
  {                                                                                                                           \
    image.clear();                                                                                                            \
    throw Error(std::string(getName())                                                                                        \
            + ": Image was successfully decoded but it has a format for which this library does not provide an interface");   \
  }                                                                                                                           \
  image._setType(type);                                                                                                       \
}

// TODO: Add options to all the ones that support them
THEA_DEF_DESERIALIZE_IMAGE(CodecBMP,     FIF_BMP,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecCUT,     FIF_CUT,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecDDS,     FIF_DDS,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecEXR,     FIF_EXR,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecGIF,     FIF_GIF,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecHDR,     FIF_HDR,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecICO,     FIF_ICO,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecIFF,     FIF_IFF,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecJ2K,     FIF_J2K,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecJNG,     FIF_JNG,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecJP2,     FIF_JP2,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecJPEG,    FIF_JPEG,    0)
THEA_DEF_DESERIALIZE_IMAGE(CodecKOALA,   FIF_KOALA,   0)
THEA_DEF_DESERIALIZE_IMAGE(CodecMNG,     FIF_MNG,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPBM,     FIF_PBM,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPBMRAW,  FIF_PBMRAW,  0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPCD,     FIF_PCD,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPCX,     FIF_PCX,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPFM,     FIF_PFM,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPGM,     FIF_PGM,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPGMRAW,  FIF_PGMRAW,  0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPNG,     FIF_PNG,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPPM,     FIF_PPM,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPPMRAW,  FIF_PPMRAW,  0)
THEA_DEF_DESERIALIZE_IMAGE(CodecPSD,     FIF_PSD,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecRAS,     FIF_RAS,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecSGI,     FIF_SGI,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecTARGA,   FIF_TARGA,   0)
THEA_DEF_DESERIALIZE_IMAGE(CodecTIFF,    FIF_TIFF,    0)
THEA_DEF_DESERIALIZE_IMAGE(CodecWBMP,    FIF_WBMP,    0)
THEA_DEF_DESERIALIZE_IMAGE(CodecXBM,     FIF_XBM,     0)
THEA_DEF_DESERIALIZE_IMAGE(CodecXPM,     FIF_XPM,     0)

int
Image::Type::numChannels() const
{
  switch (value)
  {
    case UNKNOWN   :  return -1;

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
      || value == RGB_32F
      || value == RGBA_32F
      || value == COMPLEX_64F;
}

int
Image::Type::getBitsPerPixel() const
{
  return (int)ImageInternal::typeToFreeImageBPP(*this);
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
    case 1  : return channel == Channel::ALPHA ? getBitsPerPixel() : 0;
    case 3  : return channel != Channel::ALPHA ? getBitsPerChannel() : 0;
    case 4  : return getBitsPerChannel();
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
: fip_img(new fipImage), type(Type::UNKNOWN)
{
  cacheTypeProperties();
}

Image::Image(Type type_, int width, int height)
: fip_img(NULL), type(type_)
{
  if (type == Type::UNKNOWN || width <= 0 || height <= 0)
    throw Error("Cannot initialize image of unknown type or non-positive size");

  fip_img = new fipImage(ImageInternal::typeToFreeImageType(type), width, height, ImageInternal::typeToFreeImageBPP(type));
  if (!isValid())
    throw Error("Could not create an image of the specified type and dimensions");

  cacheTypeProperties();
}

Image::Image(BinaryInputStream & input, Codec const & codec)
: fip_img(new fipImage), type(Type::UNKNOWN)
{
  deserialize(input, codec);
}

Image::Image(std::string const & filename, Codec const & codec)
: fip_img(new fipImage), type(Type::UNKNOWN)
{
  load(filename, codec);
}

Image::Image(Image const & src)
: fip_img(new fipImage(*src.fip_img)), type(src.type)
{
  cacheTypeProperties();
}

Image &
Image::operator=(Image const & src)
{
  *fip_img = *src.fip_img;
  type = src.type;
  cacheTypeProperties();

  return *this;
}

Image::~Image()
{
  delete fip_img;
}

bool
Image::isValid() const
{
  return fip_img->isValid() != 0;
}

void
Image::clear()
{
  fip_img->clear();
  type = Type::UNKNOWN;
}

void
Image::resize(Type type_, int width_, int height_)
{
  if (type_ == Type::UNKNOWN || width_ <= 0 || height_ <= 0)
    throw Error("Cannot resize image to unknown type or non-positive size (use clear() function to destroy data)");

  if (type_ == type && width_ == getWidth() && height_ == getHeight())
    return;

  fip_img->setSize(ImageInternal::typeToFreeImageType(type_), width_, height_, ImageInternal::typeToFreeImageBPP(type_));
  if (!isValid())
    throw Error("Could not resize the image to the specified type and dimensions");

  type = type_;

  cacheTypeProperties();
}

int
Image::getWidth() const
{
  return isValid() ? fip_img->getWidth() : 0;
}

int
Image::getHeight() const
{
  return isValid() ? fip_img->getHeight() : 0;
}

void const *
Image::getData() const
{
  return isValid() ? fip_img->accessPixels() : NULL;
}

void *
Image::getData()
{
  return isValid() ? fip_img->accessPixels() : NULL;
}

double
Image::getNormalizedValue(void const * pixel, int channel) const
{
  switch (type)
  {
    case Image::Type::LUMINANCE_8U  : return channel == Channel::ALPHA ? *((uint8 const *)pixel) / 255.0 : 0.0;

    case Image::Type::LUMINANCE_16  : return channel == Channel::ALPHA ? *((int16  const *)pixel) / 32768.0 : 0.0;
    case Image::Type::LUMINANCE_16U : return channel == Channel::ALPHA ? *((uint16 const *)pixel) / 65535.0 : 0.0;

    case Image::Type::LUMINANCE_32  : return channel == Channel::ALPHA ? *((int32  const *)pixel) / 2147483648.0 : 0.0;
    case Image::Type::LUMINANCE_32U : return channel == Channel::ALPHA ? *((uint32 const *)pixel) / 4294967295.0 : 0.0;

    case Image::Type::LUMINANCE_32F : return channel == Channel::ALPHA ? *((float32 const *)pixel) : 0.0;
    case Image::Type::LUMINANCE_64F : return channel == Channel::ALPHA ? *((float64 const *)pixel) : 0.0;

    case Image::Type::RGB_8U        : return channel != Channel::ALPHA ? ((uint8   const *)pixel)[channel] / 255.0 : 0.0;
    case Image::Type::RGB_16U       : return channel != Channel::ALPHA ? ((uint16  const *)pixel)[channel] / 65535.0 : 0.0;
    case Image::Type::RGB_32F       : return channel != Channel::ALPHA ? ((float32 const *)pixel)[channel] : 0.0;

    case Image::Type::RGBA_8U       : return ((uint8   const *)pixel)[channel] / 255.0;
    case Image::Type::RGBA_16U      : return ((uint16  const *)pixel)[channel] / 65535.0;
    case Image::Type::RGBA_32F      : return ((float32 const *)pixel)[channel];

    case Image::Type::COMPLEX_64F   :
    {
      if (channel == Channel::ALPHA)
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

void const *
Image::getScanLine(int row) const
{
  return isValid() ? fip_img->getScanLine(row) : NULL;
}

void *
Image::getScanLine(int row)
{
  return isValid() ? fip_img->getScanLine(row) : NULL;
}

int
Image::getScanWidth() const
{
  return isValid() ? fip_img->getScanWidth() : 0;
}

int
Image::getRowAlignment() const
{
  return 4;  // the current FreeImage default
}

namespace ImageInternal {

FREE_IMAGE_FILTER
filterToFreeImageFilter(Image::Filter filter)
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

} // namespace ImageInternal

bool
Image::rescale(int new_width, int new_height, Filter filter)
{
  if (!isValid())
  {
    THEA_ERROR << "Image: Attempting to rescale an invalid image";
    return false;
  }

  if (new_width <= 0 || new_height <= 0)
  {
    THEA_ERROR << "Image: Attempting to rescale to invalid dimensions: " << new_width << " x " << new_height;
    return false;
  }

  FREE_IMAGE_FILTER fi_filter = ImageInternal::filterToFreeImageFilter(filter);
  return (fip_img->rescale((unsigned int)new_width, (unsigned int)new_height, fi_filter) == TRUE);
}

void
Image::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  if (!isValid())
    throw Error("Can't serialize an invalid image");

  if (codec == Codec_AUTO())
    throw Error("You must explicitly choose a codec for serializing images");

  ImageCodec const * img_codec = dynamic_cast<ImageCodec const *>(&codec);
  if (!img_codec)
    throw Error("Codec specified for image serialization is not an image codec.");

  img_codec->serializeImage(*this, output, true);
}

void
Image::deserialize(BinaryInputStream & input, Codec const & codec)
{
  if (codec == Codec_AUTO())
    deserialize_AUTO(input, true);
  else
  {
    ImageCodec const * img_codec = dynamic_cast<ImageCodec const *>(&codec);
    if (!img_codec)
      throw Error("Codec specified for image deserialization is not an image codec.");

    img_codec->deserializeImage(*this, input, true);
  }
}

void
Image::save(std::string const & filename, Codec const & codec) const
{
  if (!isValid())
    throw Error("Can't save an invalid image");

  ImageCodec const * c = NULL;
  if (codec == Codec_AUTO())
  {
    c = ImageInternal::codecFromPath(filename);
    if (!c)
      throw Error("Could not detect codec from filename");
  }

  if (!c)
  {
    c = dynamic_cast<ImageCodec const *>(&codec);
    if (!c)
      throw Error("Codec specified for saving image is not an image codec.");
  }

  BinaryOutputStream out(filename, Endianness::LITTLE);
  if (!out.ok())
    throw Error("Could not open image file for writing");

  c->serializeImage(*this, out, false);

  out.commit();
  if (!out.ok())
    throw Error("Could not save image file");
}

void
Image::load(std::string const & filename, Codec const & codec)
{
  BinaryInputStream in(filename, Endianness::LITTLE);
  int64 file_size = in.size();
  if (file_size <= 0)
    throw Error("Image file does not exist or is empty");

  if (codec == Codec_AUTO())
    deserialize_AUTO(in, false);
  else
  {
    try
    {
      ImageCodec const & img_codec = dynamic_cast<ImageCodec const &>(codec);
      img_codec.deserializeImage(*this, in, false);
    }
    catch (std::bad_cast &)
    {
      throw Error("Codec specified for loading image is not an image codec.");
    }
  }
}

void
Image::deserialize_AUTO(BinaryInputStream & input, bool read_prefixed_info)
{
  // Get the size of the image block in bytes
  input.setEndianness(Endianness::LITTLE);
  uint32 size = read_prefixed_info ? input.readUInt32() : input.size();

  // Read the image block into a memory buffer (optimization possible when the data has already been buffered within the input
  // stream?)
  TheaArray<uint8> img_block((array_size_t)size);
  input.readBytes((int64)size, &img_block[0]);

  // Decode the image
  fipMemoryIO mem((BYTE *)&img_block[0], (DWORD)size);
  if (!fip_img->loadFromMemory(mem))
    throw Error("Could not load image from memory stream");

  type = ImageInternal::typeFromFreeImageTypeAndBPP(fip_img->getImageType(), fip_img->getBitsPerPixel());
  if (type == Type::UNKNOWN)
  {
    clear();
    throw Error("Image was successfully decoded but it has a format for which this library does not provide an interface");
  }

  cacheTypeProperties();
}

void
Image::_setType(Type type_)
{
  type = type_;
  cacheTypeProperties();
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
  return (fip_img->invert() == TRUE);
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
    if (&dst != this) dst = *this;
    status = true;
  }
  else
  {
    FREE_IMAGE_TYPE src_fitype = ImageInternal::typeToFreeImageType(type);
    FREE_IMAGE_TYPE dst_fitype = ImageInternal::typeToFreeImageType(dst_type);

    if (dst_fitype == FIT_BITMAP)
    {
      WORD dst_fibpp = ImageInternal::typeToFreeImageBPP(dst_type);
      switch (dst_fibpp)
      {
        case 4:
        {
          if (&dst != this) dst = *this;

          if (src_fitype != FIT_BITMAP) status = (dst.fip_img->convertToType(FIT_BITMAP) == TRUE);
          if (status)                   status = (dst.fip_img->convertTo4Bits() == TRUE);

          break;
        }

        case 8:
        {
          if (&dst != this) dst = *this;

          if (src_fitype != FIT_BITMAP)
            status = (dst.fip_img->convertToType(FIT_BITMAP) == TRUE);
          else
            status = true;

          if (status)
            status = (dst.fip_img->convertToGrayscale() == TRUE);  // convertTo8Bits() can palletize

          break;
        }

        case 24:
        {
          if (&dst != this) dst = *this;
          status = (dst.fip_img->convertTo24Bits() == TRUE);
          break;
        }

        case 32:
        {
          if (&dst != this) dst = *this;
          status = (dst.fip_img->convertTo32Bits() == TRUE);
          break;
        }
      }
    }
    else
    {
      if (&dst != this) dst = *this;

      if (dst.fip_img->convertToType(dst_fitype) == TRUE)
        status = true;
    }

    if (status)
    {
      dst.type = dst_type;
      dst.cacheTypeProperties();
    }
  }

  return status;
}

} // namespace Thea
