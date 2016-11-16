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

 @file BinaryInputStream.cpp

 @author Morgan McGuire, graphics3d.com
 Copyright 2001-2007, Morgan McGuire.  All rights reserved.

 @created 2001-08-09
 @edited  2010-03-05


  <PRE>
    {
    BinaryOutput b("c:/tmp/test.b", BinaryOutput::LITTLE_ENDIAN);

    float f = 3.1415926;
    int i = 1027221;
    std::string s = "Hello World!";

    b.writeFloat32(f);
    b.writeInt32(i);
    b.writeString(s);
    b.commit();


    BinaryInputStream in("c:/tmp/test.b", BinaryInputStream::LITTLE_ENDIAN);

    debugAssert(f == in.readFloat32());
    int ii = in.readInt32();
    debugAssert(i == ii);
    debugAssert(s == in.readString());
    }
  </PRE>
 */

#include "BinaryInputStream.hpp"
#include "FilePath.hpp"
#include "FileSystem.hpp"
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>

THEA_INSTANTIATE_SMART_POINTERS(Thea::BinaryInputStream)

namespace Thea {

bool const BinaryInputStream::NO_COPY = false;

// The initial buffer will be no larger than this (50 MB), but may grow if a large memory read occurs.
#define THEA_INITIAL_READ_BUFFER_LENGTH 50000000

void
BinaryInputStream::readBool8(int64 n, TheaArray<bool> & out)
{
  out.resize((array_size_t)n);

  // std::vector optimizes bool in a way that prevents fast reading
  for (int64 i = 0; i < n ; ++i)
    out[(array_size_t)i] = readBool8();
}

#define THEA_BINARY_INPUT_STREAM_DEFINE_READER(fname, tname) \
  void BinaryInputStream::read##fname(int64 n, TheaArray<tname> & out) \
  { \
    out.resize((array_size_t)n); \
    read##fname(n, &out[0]); \
  }

THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt8,               uint8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int8,                int8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt16,              uint16)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int16,               int16)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt32,              uint32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int32,               int32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt64,              uint64)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int64,               int64)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Float32,             float32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Float64,             float64)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector2,             Vector2)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector3,             Vector3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector4,             Vector4)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorL8,             ColorL8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorL,              ColorL)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGB8,           ColorRGB8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGB,            ColorRGB)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGBA8,          ColorRGBA8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGBA,           ColorRGBA)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix2,             Matrix2)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix3,             Matrix3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix4,             Matrix4)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(CoordinateFrame3,    CoordinateFrame3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Plane3,              Plane3)

#undef THEA_BINARY_INPUT_STREAM_DEFINE_READER

// Data structures that are one byte per element can be directly copied, regardles of endianness.
#define THEA_BINARY_INPUT_STREAM_DEFINE_READER(fname, tname) \
  void BinaryInputStream::read##fname(int64 n, tname * out) \
  { \
    if (sizeof(tname) == 1) \
      readBytes(n, out); \
    else \
    { \
      for (int64 i = 0; i < n ; ++i) \
        out[i] = read##fname(); \
    } \
  }

THEA_BINARY_INPUT_STREAM_DEFINE_READER(Bool8,               bool)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt8,               uint8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int8,                int8)

#undef THEA_BINARY_INPUT_STREAM_DEFINE_READER

#define THEA_BINARY_INPUT_STREAM_DEFINE_READER(fname, tname) \
  void BinaryInputStream::read##fname(int64 n, tname * out) \
  { \
    if (m_swapBytes) \
    { \
      for (int64 i = 0; i < n; ++i) \
        out[i] = read##fname(); \
    } \
    else \
      readBytes(sizeof(tname) * n, out); \
  }

THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt16,              uint16)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int16,               int16)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt32,              uint32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int32,               int32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(UInt64,              uint64)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Int64,               int64)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Float32,             float32)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Float64,             float64)

#undef THEA_BINARY_INPUT_STREAM_DEFINE_READER

#define THEA_BINARY_INPUT_STREAM_DEFINE_READER(fname, tname) \
  void BinaryInputStream::read##fname(int64 n, tname * out) \
  { \
    for (int64 i = 0; i < n; ++i) \
      out[i] = read##fname(); \
  }

THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector2,             Vector2)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector3,             Vector3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Vector4,             Vector4)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorL8,             ColorL8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorL,              ColorL)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGB8,           ColorRGB8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGB,            ColorRGB)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGBA8,          ColorRGBA8)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(ColorRGBA,           ColorRGBA)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix2,             Matrix2)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix3,             Matrix3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Matrix4,             Matrix4)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(CoordinateFrame3,    CoordinateFrame3)
THEA_BINARY_INPUT_STREAM_DEFINE_READER(Plane3,              Plane3)

#undef THEA_BINARY_INPUT_STREAM_DEFINE_READER

void
BinaryInputStream::loadIntoMemory(int64 start_position, int64 min_length)
{
  // Load the next section of the file
  debugAssertM(m_path != "<memory>", getNameStr() + ": Read past end of stream");
  int64 absPos = m_alreadyRead + m_pos;

  if (m_bufferLength < min_length)
  {
    // The current buffer isn't big enough to hold the chunk we want to read. This happens if there was little memory available
    // during the initial constructor read but more memory has since been freed.
    m_bufferLength = min_length;
    alwaysAssertM(m_freeBuffer, getNameStr() + ": Need to reallocate buffer, but it is read-only");

    m_buffer = (uint8 *)std::realloc(m_buffer, m_bufferLength);
    if (!m_buffer)
      throw Error(getNameStr() + ": Tried to read a larger memory chunk than could fit in memory");
  }

  m_alreadyRead = start_position;

#ifdef THEA_WINDOWS

  FILE * file = fopen(m_path.c_str(), "rb");
  if (!file)
    throw Error(getNameStr() + ": Could not open file");

  int ret_seek = fseek(file, (off_t)m_alreadyRead, SEEK_SET);
  if (ret_seek != 0)
    throw Error(getNameStr() + ": Seek failed");

  size_t toRead = (size_t)std::min(m_bufferLength, m_length - m_alreadyRead);
  size_t ret_read = fread(m_buffer, 1, toRead, file);
  if (ret_read != toRead)
    throw Error(getNameStr() + ": Read failed");

  fclose(file);
  file = NULL;

#else

  FILE * file = fopen(m_path.c_str(), "rb");
  if (!file)
    throw Error(getNameStr() + ": Could not open file");

  int ret = fseeko(file, (off_t)m_alreadyRead, SEEK_SET);
  if (ret != 0)
    throw Error(getNameStr() + ": Seek failed");

  size_t toRead = (size_t)std::min(m_bufferLength, m_length - m_alreadyRead);
  ret = fread(m_buffer, 1, toRead, file);
  if ((size_t)ret != (size_t)toRead)
    throw Error(getNameStr() + ": Read failed");

  fclose(file);
  file = NULL;

#endif

  m_pos = absPos - m_alreadyRead;
  if (m_pos < 0)
    throw Error(getNameStr() + ": Could not set initial read position");
}

void
BinaryInputStream::setEndianness(Endianness e)
{
  m_fileEndian = e;
  m_swapBytes = (e != Endianness::machine());
}

BinaryInputStream::BinaryInputStream(uint8 const * data, int64 data_len, Endianness data_endian, bool copy_memory)
: NamedObject("<memory>"),
  m_path("<memory>"),
  m_bitPos(0),
  m_bitString(0),
  m_beginEndBits(0),
  m_alreadyRead(0),
  m_bufferLength(0),
  m_pos(0)
{
  m_freeBuffer = copy_memory;
  setEndianness(data_endian);

  m_length = data_len;
  m_bufferLength = m_length;

  if (!copy_memory)
  {
    debugAssertM(!m_freeBuffer, "BinaryInputStream: Can't free buffer if it is not a copy");
    m_buffer = const_cast<uint8 *>(data);
  }
  else
  {
    debugAssertM(m_freeBuffer, "BinaryInputStream: Must free buffer if it is a copy");
    m_buffer = (uint8 *)std::malloc(m_length);
    std::memcpy(m_buffer, data, data_len);
  }
}

BinaryInputStream::BinaryInputStream(std::string const & path, Endianness file_endian)
: NamedObject(FilePath::objectName(path)),
  m_path(FileSystem::resolve(path)),
  m_bitPos(0),
  m_bitString(0),
  m_beginEndBits(0),
  m_alreadyRead(0),
  m_length(0),
  m_bufferLength(0),
  m_buffer(NULL),
  m_pos(0),
  m_freeBuffer(true)
{
  setEndianness(file_endian);

  // Figure out how big the file is and verify that it exists.
  m_length = FileSystem::fileSize(m_path);

  // Open the file
  FILE * file = fopen(m_path.c_str(), "rb");

  if (!file || (m_length == -1))
    throw Error("BinaryInputStream: File '" + m_path + "' not found");

  // Read part or all of the file into the memory buffer
  if (m_length > THEA_INITIAL_READ_BUFFER_LENGTH)
  {
    // Read only a subset of the file so we don't consume
    // all available memory.
    m_bufferLength = THEA_INITIAL_READ_BUFFER_LENGTH;
  }
  else
  {
    // Either the length is fine or the file is compressed
    // and requires us to read the whole thing for zlib.
    m_bufferLength = m_length;
  }

  alwaysAssertM(m_freeBuffer, "BinaryInputStream: Allocated buffer not set to be freed");
  m_buffer = (uint8 *)std::malloc(m_bufferLength);

  if (m_buffer == NULL)
  {
    // Try to allocate a small array; not much memory is available.
    // Give up if we can't allocate even 1k.
    while ((m_buffer == NULL) && (m_bufferLength > 1024))
    {
      m_bufferLength /= 2;
      m_buffer = (uint8 *)std::malloc(m_bufferLength);
    }
  }

  if (!m_buffer)
    throw Error(getNameStr() + ": Could not allocate buffer");

  size_t num_read = fread(m_buffer, sizeof(int8), (size_t)m_bufferLength, file);
  if (num_read != (size_t)m_bufferLength)
    throw Error(getNameStr() + ": Could not initialize file buffer");

  fclose(file);
  file = NULL;
}

BinaryInputStream::~BinaryInputStream()
{
  if (m_freeBuffer)
    std::free(m_buffer);
}

void
BinaryInputStream::readBytes(int64 n, void * bytes)
{
  prepareToRead(n);
  std::memcpy(bytes, m_buffer + m_pos, n);
  m_pos += n;
}

uint64
BinaryInputStream::readUInt64()
{
  prepareToRead(8);
  m_pos += 8;

  if (m_swapBytes)
  {
    uint8 out[8];
    out[0] = m_buffer[m_pos - 1];
    out[1] = m_buffer[m_pos - 2];
    out[2] = m_buffer[m_pos - 3];
    out[3] = m_buffer[m_pos - 4];
    out[4] = m_buffer[m_pos - 5];
    out[5] = m_buffer[m_pos - 6];
    out[6] = m_buffer[m_pos - 7];
    out[7] = m_buffer[m_pos - 8];

    return *(uint64 *)out;
  }
  else
  {
#ifdef THEA_ALLOW_UNALIGNED_READS
    return *(uint64 *)(m_buffer + m_pos);
#else
    uint8 out[8];
    out[0] = m_buffer[m_pos - 8];
    out[1] = m_buffer[m_pos - 7];
    out[2] = m_buffer[m_pos - 6];
    out[3] = m_buffer[m_pos - 5];
    out[4] = m_buffer[m_pos - 4];
    out[5] = m_buffer[m_pos - 3];
    out[6] = m_buffer[m_pos - 2];
    out[7] = m_buffer[m_pos - 1];

    return *(uint64 *)out;
#endif
  }
}

std::string
BinaryInputStream::readString(int64 n)
{
  prepareToRead(n);

  char * s = (char *)std::malloc(n + 1);
  if (!s)
    throw Error(getNameStr() + ": Could not allocate buffer for reading string");

  std::memcpy(s, m_buffer + m_pos, n);

  // There may not be a null, so make sure
  // we add one.
  s[n] = '\0';
  std::string out = s;

  std::free(s);
  s = NULL;

  m_pos += n;
  return out;
}

std::string
BinaryInputStream::readNullTerminatedString()
{
  int64 n = 0;

  if ((m_pos + m_alreadyRead + n) < (m_length - 1))
  {
    prepareToRead(1);
  }

  if ( ((m_pos + m_alreadyRead + n) < (m_length - 1)) &&
       (m_buffer[m_pos + n] != '\0') )
  {
    ++n;

    while ( ((m_pos + m_alreadyRead + n) < (m_length - 1)) &&
            (m_buffer[m_pos + n] != '\0') )
    {
      ++n;
      prepareToRead(n);
    }
  }

  // Consume NULL
  ++n;
  return readString(n);
}

std::string
BinaryInputStream::readLine()
{
  int64 n = 0;

  if ((m_pos + m_alreadyRead + n) < (m_length - 1))
  {
    prepareToRead(1);
  }

  if ( ((m_pos + m_alreadyRead + n) < (m_length - 1)) &&
       ! isNewline(m_buffer[m_pos + n]))
  {
    ++n;

    while ( ((m_pos + m_alreadyRead + n) < (m_length - 1)) &&
            ! isNewline(m_buffer[m_pos + n]))
    {
      ++n;
      prepareToRead(n);
    }
  }

  std::string const s = readString(n);

  // Consume the newline
  if (hasMore())  // if we haven't reached the end of the file, then we've encountered a newline
  {
    char firstNLChar = readInt8();

    // Consume the 2nd newline character, if it exists
    if (hasMore())
    {
      prepareToRead(1);
      if (isNewline(m_buffer[m_pos]) && (m_buffer[m_pos] != firstNLChar))
        skip(1);
    }
  }

  return s;
}

std::string
BinaryInputStream::readAlignedString(int alignment)
{
  int length = (int)readUInt32();  // string length
  std::string s = readString(length);  // string chars

  int padding = (alignment - (length % alignment)) % alignment;
  if (padding > 0)
    skip(padding);

  return s;
}

void
BinaryInputStream::beginBits()
{
  alwaysAssertM(m_beginEndBits == 0, getNameStr() + ": beginBits/endBits mismatch");

  m_beginEndBits = 1;
  m_bitPos = 0;

  if (!hasMore())
    throw Error(getNameStr() + ": Can't call beginBits() when at the end of a file");

  m_bitString = readUInt8();
}

uint32
BinaryInputStream::readBits(int num_bits)
{
  alwaysAssertM(m_beginEndBits == 1, getNameStr() + ": Can't call readBits outside a beginBits/endBits block");

  uint32 out = 0;
  int const total = num_bits;

  while (num_bits > 0)
  {
    if (m_bitPos > 7)
    {
      // Consume a new byte for reading.  We do this at the beginning
      // of the loop so that we don't try to read past the end of the file.
      m_bitPos = 0;
      m_bitString = readUInt8();
    }

    // Slide the lowest bit of the bitString into
    // the correct position.
    out |= (m_bitString & 1) << (total - num_bits);

    // Shift over to the next bit
    m_bitString = m_bitString >> 1;
    ++m_bitPos;
    --num_bits;
  }

  return out;
}

void
BinaryInputStream::endBits()
{
  alwaysAssertM(m_beginEndBits == 1, getNameStr() + ": beginBits/endBits mismatch");

  if (m_bitPos == 0)
  {
    // Put back the last byte we read
    --m_pos;
  }

  m_beginEndBits = 0;
  m_bitPos = 0;
}

} // namespace Thea
