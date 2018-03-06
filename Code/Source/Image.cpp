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
#include "Algorithms/FastCopy.hpp"
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

// Returns bytes consumed by an aligned row
int
scanWidth(int width, int bpp, int alignment)
{
  int width_bits = width * bpp;
  int row_bytes = width_bits / 8 + (width_bits % 8 == 0 ? 0 : 1);
  return row_bytes + (alignment - (row_bytes % alignment)) % alignment;
}

static ImageCodec const *
codecFromPath(std::string const & path)
{
  typedef shared_ptr<ImageCodec> ICPtr;
  typedef TheaUnorderedMap<std::string, ICPtr > CodecMap;

  static CodecMap codec_map;
  if (codec_map.empty())  // only on the first call
  {
    // 2D formats
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

    // 3D formats
    codec_map["3BM"  ] = ICPtr(new Codec3BM());
  }

  CodecMap::const_iterator existing = codec_map.find(toUpper(FilePath::extension(path)));
  if (existing == codec_map.end())
    return NULL;
  else
    return existing->second.get();
}

// Currently only detects 3D formats, the rest are handled by FreeImage
ImageCodec const *
codecFromMagic(int64 num_bytes, uint8 const * buf)
{
  static Codec3BM const CODEC_3BM;

  if (num_bytes >= 4 && buf[0] == (uint8)'3' && buf[1] == (uint8)'B' && buf[2] == (uint8)'M' && buf[3] == (uint8)'\0')
    return &CODEC_3BM;

  return NULL;
}

} // namespace ImageInternal

//=============================================================================================================================
// 2D formats
//=============================================================================================================================

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
  TheaArray<uint8> img_block((size_t)size);                                                                             \
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

//=============================================================================================================================
// 3D formats
//=============================================================================================================================

long
Codec3BM::serializeImage(Image const & image, BinaryOutputStream & output, bool prefix_info) const
{
  if (!image.isValid())
    throw Error(std::string(getName()) + ": Cannot serialize an invalid image");

  Image::Type type = image.getType();
  int width   =  image.getWidth();
  int height  =  image.getHeight();
  int depth   =  image.getDepth();

  if (type != Image::Type::LUMINANCE_8U
   && type != Image::Type::RGB_8U
   && type != Image::Type::RGBA_8U)
    throw Error(std::string(getName()) + ": Can only serialize 8-bit luminance, RGB or RGBA images");

  output.setEndianness(Endianness::LITTLE);

  int bpp = image.getBitsPerPixel();
  int stream_scan_width = ImageInternal::scanWidth(width, bpp, 16);  // 16-byte row alignment for SSE compatibility

  static uint64 const FILE_HEADER_SIZE = 32;  // Remember to change these if the format changes!
  static uint64 const INFO_HEADER_SIZE = 72;
  uint64 data_size = (uint64)stream_scan_width * (uint64)height * (uint64)depth;
  uint64 size_in_bytes = FILE_HEADER_SIZE + INFO_HEADER_SIZE + data_size;

  if (prefix_info)
    output.writeUInt64(size_in_bytes);

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
  //   4 bytes: Number of bits per pixel (1, 4, 8, 16, 24, 32), currently only 8, 24 and 32 are supported
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
  output.writeUInt64(data_size);
  output.writeUInt64(3937);  // 100dpi
  output.writeUInt64(3937);  // 100dpi
  output.writeUInt64(3937);  // 100dpi

  // Write image data as L, BGR or BGRA voxels, one 16-byte aligned scanline at a time
  if (width > 0 && height > 0 && depth > 0 && bpp > 0)
  {
    int bytes_pp = bpp / 8;
    TheaArray<uint8> row_buf((size_t)stream_scan_width, 0);
    for (int i = 0; i < depth; ++i)
      for (int j = 0; j < height; ++j)
      {
        uint8 const * in_pixel = (uint8 const *)image.getScanLine(j, i);
        uint8 * out_pixel = &row_buf[0];

        switch (bytes_pp)
        {
          case 1: Algorithms::fastCopy(in_pixel, in_pixel + width * bytes_pp, out_pixel); break;
          case 3:
          {
            for (int k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
            {
              out_pixel[0] = in_pixel[Image::Channel::BLUE ];
              out_pixel[1] = in_pixel[Image::Channel::GREEN];
              out_pixel[2] = in_pixel[Image::Channel::RED  ];
            }
            break;
          }
          default: // 4
          {
            for (int k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
            {
              out_pixel[0] = in_pixel[Image::Channel::BLUE ];
              out_pixel[1] = in_pixel[Image::Channel::GREEN];
              out_pixel[2] = in_pixel[Image::Channel::RED  ];
              out_pixel[3] = in_pixel[Image::Channel::ALPHA];
            }
            break;
          }
        }

        output.writeBytes((int64)stream_scan_width, &row_buf[0]);
      }
  }

  return (long)(prefix_info ? size_in_bytes + 4 : size_in_bytes);
}

void
Codec3BM::deserializeImage(Image & image, BinaryInputStream & input, bool read_prefixed_info) const
{
  // Get the size of the image block in bytes
  input.setEndianness(Endianness::LITTLE);
  uint64 prefixed_size = read_prefixed_info ? input.readUInt64() : 0 /* we'll fix this below */;

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
    throw Error(std::string(getName()) + ": Image is not in 3BM format");

  input.skip(4);
  uint64 size = input.readUInt64();
  if (prefixed_size > 0 && size != prefixed_size)
    throw Error(std::string(getName()) + ": Prefixed size does not match image size from header");

  input.skip(8);
  uint64 bitmap_offset = input.readUInt64();

  // Read info header (72 bytes):
  //   8 bytes: Size of this header (72 bytes)
  //   8 bytes: Width of the bitmap in voxels
  //   8 bytes: Height of the bitmap in voxels
  //   8 bytes: Depth of the bitmap in voxels
  //   4 bytes: Number of bits per pixel (1, 4, 8, 16, 24, 32), currently only 8, 24 and 32 are supported
  //   4 bytes: Compression method (0 for no compression, currently this is the only one supported)
  //   8 bytes: Size of raw image data
  //   8 bytes: Width resolution in voxels/meter
  //   8 bytes: Height resolution in voxels/meter
  //   8 bytes: Depth resolution in voxels/meter

  input.skip(8);
  int width   =  (int)input.readUInt64();
  int height  =  (int)input.readUInt64();
  int depth   =  (int)input.readUInt64();

  int bpp = (int)input.readUInt32();
  if (bpp != 8 && bpp != 24 && bpp != 32)
    throw Error(std::string(getName()) + ": Only 8, 24 and 32-bit bitmaps currently supported");

  uint32 compression = input.readUInt32();
  if (compression != 0)
    throw Error(std::string(getName()) + ": Unsupported compression method");

  int stream_scan_width = ImageInternal::scanWidth(width, bpp, 16);  // 16-byte row alignment for SSE compatibility
  uint64 data_size = input.readUInt64();
  if (data_size != (uint64)stream_scan_width * (uint64)height * (uint64)depth)
    throw Error(std::string(getName()) + ": Pixel data block size does not match image dimensions");

  // Read image data as L, BGR or BGRA voxels, one 16-byte aligned scanline at a time
  input.setPosition(start_position + (int64)bitmap_offset);

  Image::Type type = Image::Type::UNKNOWN;
  switch (bpp)
  {
    case 8   :  type = Image::Type::LUMINANCE_8U; break;
    case 24  :  type = Image::Type::RGB_8U; break;
    case 32  :  type = Image::Type::RGBA_8U; break;
    default  :  throw Error(format("%s: Bit depth %d not supported", getName(), bpp));
  }

  image.resize(type, width, height, depth);
  if (width <= 0 || height <= 0 || depth <= 0 || bpp <= 0)
    return;

  int bytes_pp = bpp / 8;
  TheaArray<uint8> row_buf((size_t)stream_scan_width);
  for (int i = 0; i < depth; ++i)
    for (int j = 0; j < height; ++j)
    {
      input.readBytes(stream_scan_width, &row_buf[0]);

      uint8 const * in_pixel = &row_buf[0];
      uint8 * out_pixel = (uint8 *)image.getScanLine(j, i);

      switch (bytes_pp)
      {
        case 1: Algorithms::fastCopy(in_pixel, in_pixel + width * bytes_pp, out_pixel); break;
        case 3:
        {
          for (int k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
          {
            out_pixel[Image::Channel::RED  ] = in_pixel[2];
            out_pixel[Image::Channel::GREEN] = in_pixel[1];
            out_pixel[Image::Channel::BLUE ] = in_pixel[0];
          }
          break;
        }
        default: // 4
        {
          for (int k = 0; k < width; ++k, in_pixel += bytes_pp, out_pixel += bytes_pp)
          {
            out_pixel[Image::Channel::RED  ] = in_pixel[2];
            out_pixel[Image::Channel::GREEN] = in_pixel[1];
            out_pixel[Image::Channel::BLUE ] = in_pixel[0];
            out_pixel[Image::Channel::ALPHA] = in_pixel[3];
          }
          break;
        }
      }
    }
}

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
: type(Type::UNKNOWN), width(0), height(0), depth(0), fip_img(NULL)
{
  cacheTypeProperties();
}

Image::Image(Type type_, int width_, int height_, int depth_)
: type(Type::UNKNOWN), width(0), height(0), depth(0), fip_img(NULL)
{
  resize(type_, width_, height_, depth_);
}

Image::Image(BinaryInputStream & input, Codec const & codec)
: type(Type::UNKNOWN), fip_img(NULL)
{
  deserialize(input, codec);
}

Image::Image(std::string const & path, Codec const & codec)
: type(Type::UNKNOWN), fip_img(NULL)
{
  load(path, codec);
}

Image::Image(Image const & src)
: type(Type::UNKNOWN), fip_img(NULL)
{
  *this = src;
  cacheTypeProperties();
}

Image &
Image::operator=(Image const & src)
{
  type = src.type;
  width = src.width;
  height = src.height;
  depth = src.depth;

  if (src.depth == 1)
  {
    if (src.fip_img)
    {
      if (!fip_img)
        fip_img = new fipImage(*src.fip_img);
      else
        *fip_img = *src.fip_img;
    }
    else
      fip_img = NULL;

    data.clear();
  }
  else
  {
    data = src.data;
    delete fip_img; fip_img = NULL;
  }

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
  return type != Type::UNKNOWN && width >= 0 && height >= 0 && depth >= 0
      && (depth != 1 || (fip_img && fip_img->isValid() != 0));
}

void
Image::clear()
{
  type = Type::UNKNOWN;
  width = height = depth = 0;

  if (fip_img)
  {
    fip_img->clear();
    fip_img = NULL;
  }

  data.clear();
}

void
Image::resize(Type type_, int width_, int height_, int depth_)
{
  if (type_ == Type::UNKNOWN || width_ <= 0 || height_ <= 0 || depth_ <= 0)
    throw Error("Cannot resize image to unknown type or non-positive size (use clear() function to destroy data)");

  if (type_ == type && width_ == getWidth() && height_ == getHeight() && depth_ == getDepth())
    return;

  if (depth_ == 1)
  {
    if (!fip_img)
      fip_img = new fipImage;

    fip_img->setSize(ImageInternal::typeToFreeImageType(type_), width_, height_, ImageInternal::typeToFreeImageBPP(type_));
  }
  else
  {
    if (type_.getBitsPerPixel() % 8 != 0)
      throw Error("Non-2D image must have byte-aligned pixels");

    int scan_width = ImageInternal::scanWidth(width_, type_.getBitsPerPixel(), (int)ROW_ALIGNMENT);
    int64 buf_size = scan_width * height_ * depth_;
    data.resize((size_t)buf_size);
  }

  type = type_;
  width = width_;
  height = height_;
  depth = depth_;

  cacheTypeProperties();

  if (!isValid())
    throw Error("Could not resize the image to the specified type and dimensions");
}

void const *
Image::getData() const
{
  return const_cast<Image *>(this)->getData();
}

void *
Image::getData()
{
  if (depth == 1)
    return isValid() ? fip_img->accessPixels() : NULL;
  else
    return data.empty() ? NULL : &data[0];
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
Image::getScanLine(int row, int z) const
{
  return const_cast<Image *>(this)->getScanLine(row, z);
}

void *
Image::getScanLine(int row, int z)
{
  alwaysAssertM(z >= 0 && z < depth, "Image: Z value out of bounds");

  if (!isValid())
    return NULL;

  if (depth == 1)
    return fip_img->getScanLine(row);
  else
  {
    int scan_width = ImageInternal::scanWidth(width, type.getBitsPerPixel(), (int)ROW_ALIGNMENT);
    return &data[(z * height + row) * scan_width];
  }
}

int
Image::getScanWidth() const
{
  if (depth == 1)
    return isValid() ? fip_img->getScanWidth() : 0;
  else
    return width + ((int)ROW_ALIGNMENT - (width % (int)ROW_ALIGNMENT)) % (int)ROW_ALIGNMENT;
}

int
Image::getRowAlignment() const
{
  return depth == 1 ? 4 : ROW_ALIGNMENT;  // the current FreeImage default is 4
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
Image::rescale(int new_width, int new_height, int new_depth, Filter filter)
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
Image::save(std::string const & path, Codec const & codec) const
{
  if (!isValid())
    throw Error("Can't save an invalid image");

  ImageCodec const * c = NULL;
  if (codec == Codec_AUTO())
  {
    c = ImageInternal::codecFromPath(path);
    if (!c)
      throw Error("Could not detect codec from path");
  }

  if (!c)
  {
    c = dynamic_cast<ImageCodec const *>(&codec);
    if (!c)
      throw Error("Codec specified for saving image is not an image codec.");
  }

  BinaryOutputStream out(path, Endianness::LITTLE);
  if (!out.ok())
    throw Error("Could not open image file for writing");

  c->serializeImage(*this, out, false);

  out.commit();
  if (!out.ok())
    throw Error("Could not save image file");
}

void
Image::load(std::string const & path, Codec const & codec)
{
  BinaryInputStream in(path, Endianness::LITTLE);
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

  if (size <= 0)
    throw Error("No image data found");

  // Read the image block into a memory buffer (optimization possible when the data has already been buffered within the input
  // stream?)
  TheaArray<uint8> img_block((size_t)size);
  input.readBytes((int64)size, &img_block[0]);

  ImageCodec const * detected_codec = ImageInternal::codecFromMagic((int64)size, &img_block[0]);
  if (detected_codec)
  {
    BinaryInputStream mem_stream(&img_block[0], (int64)size, Endianness::LITTLE, /* copy_memory = */ false);
    detected_codec->deserializeImage(*this, mem_stream, false);
  }
  else
  {
    // Decode via FreeImage
    if (!fip_img)
      fip_img = new fipImage;

    fipMemoryIO mem((BYTE *)&img_block[0], (DWORD)size);
    if (!fip_img->loadFromMemory(mem))
      throw Error("Could not load image from memory stream");

    type = ImageInternal::typeFromFreeImageTypeAndBPP(fip_img->getImageType(), fip_img->getBitsPerPixel());
    if (type == Type::UNKNOWN)
    {
      clear();
      throw Error("Image was successfully decoded but it has a format for which this library does not provide an interface");
    }

    width = fip_img->getWidth();
    height = fip_img->getHeight();
    depth = 1;
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
  else if (depth != 1)
  {
    // TODO
    THEA_ERROR << "Image: Format conversion of non-2D images currently not supported";
    return false;
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
