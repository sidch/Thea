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

 @file BinaryOutputStream.h

 @maintainer Morgan McGuire, graphics3d.com

 @created 2001-08-09
 @edited  2008-01-24

 Copyright 2000-2006, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_BinaryOutputStream_hpp__
#define __Thea_BinaryOutputStream_hpp__

#ifdef _MSC_VER
// Disable conditional expression is constant, which occurs incorrectly on inlined functions
#  pragma warning(push)
#  pragma warning( disable : 4127 )
#endif

#include "Common.hpp"
#include "Array.hpp"
#include "Colors.hpp"
#include "CoordinateFrame3.hpp"
#include "MatrixMN.hpp"
#include "NamedObject.hpp"
#include "Noncopyable.hpp"
#include "Plane3.hpp"
#include "VectorN.hpp"
#include <algorithm>
#include <cstring>

namespace Thea {

/**
 * Sequential or random access output to binary files/memory. For every <code>writeX</code> method there are also versions that
 * operate on a whole TheaArray or C-array.
 *
 * Any method call can trigger an out-of-memory error when writing to a memory buffer if it can't be expanded any more.
 *
 * This class does <b>not</b> use overloading to distinguish between writing, say, a float and a double. Instead, the write
 * functions are explicitly named writeFloat32(), writeFloat64() etc. This is because it would be very hard to debug the error
 * sequence: <code>bo.write(1.0); ... float f = bi.readFloat32();</code> in which a double is serialized and then deserialized
 * as a float.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 *
 * @todo Reimplement using %<iostream%> for arbitrary seeking and safer performance?
 */
class THEA_API BinaryOutputStream : public virtual NamedObject, private Noncopyable
{
  private:
    /** Is the file big or little endian? */
    Endianness      m_fileEndian;

    /** Path to the opened file. */
    std::string     m_path;

    /** 0 outside of beginBits...endBits, 1 inside */
    int             m_beginEndBits;

    /**
     * The current string of bits being built up by beginBits() ... endBits(). This string is treated semantically, as if the
     * lowest bit was on the left and the highest was on the right.
     */
    uint8           m_bitString;

    /** Position (from the lowest bit) currently used in m_bitString. */
    int             m_bitPos;

    /** True if the file endianess does not match the machine endian. */
    bool            m_swapBytes;

    /** The buffer of bytes to be written to the output stream. */
    uint8     *     m_buffer;

    /** Number of elements in the buffer */
    int64           m_bufferLen;

    /** Underlying size of memory allocated. */
    int64           m_bufferCapacity;

    /** Next byte in file */
    int64           m_pos;

    /** Number of bytes already written to the file.*/
    int64           m_alreadyWritten;

    /** Error-check. */
    bool            m_ok;

    /** Reserve space by dumping buffer contents to disk if necessary. */
    void reserveBytesWhenOutOfMemory(size_t bytes);

    /** Resize the buffer. */
    void reallocBuffer(size_t bytes, size_t oldBufferLen);

    /** Make sure at least \a num_bytes bytes can be written, resizing if necessary. */
    void reserveBytes(int64 num_bytes)
    {
      debugAssertM(num_bytes > 0, getNameStr() + ": Can't reserve less than one byte");
      size_t oldBufferLen = (size_t)m_bufferLen;
      m_bufferLen = std::max(m_bufferLen, (m_pos + num_bytes));

      if (m_bufferLen > m_bufferCapacity)
        reallocBuffer((size_t)num_bytes, oldBufferLen);
    }

    /** Commit data to disk, optionally forcing a write even if the buffer is empty. */
    bool _commit(bool flush, bool force);

    // Not implemented on purpose, don't use.
    bool operator==(BinaryOutputStream const &);

  public:
    THEA_DEF_POINTER_TYPES(BinaryOutputStream, shared_ptr, weak_ptr)

    /** Construct a stream that writes to an (expanding, contiguous) memory buffer. */
    explicit BinaryOutputStream(Endianness endian = Endianness::LITTLE);

    /**
     * Construct a stream that writes to a file. Use "<memory>" as the path if you're going to commit to memory -- this has the
     * same effect as the default constructor. The file and its parent directories will be created if they do not exist, and the
     * file initialized to be blank. If the file cannot be constructed, ok() will return false.
     */
    BinaryOutputStream(std::string const & path, Endianness file_endian);

    /** Destructor. Automatically calls commit(). */
    ~BinaryOutputStream();

    /** True if no errors have been encountered.*/
    bool ok() const;

    /** Set the endianness of subsequent multi-byte write operations. */
    void setEndianness(Endianness endian);

    /** Get the endianness of current multi-byte write operations. */
    Endianness getEndianness() const
    {
      return m_fileEndian;
    }

    /** Get the path to the current file being written ("<memory>" for memory streams). */
    std::string getPath() const
    {
      return m_path;
    }

    /**
     * Write the buffered data to disk. Parent directories are created as needed if they do not exist. Calling commit() multiple
     * times without intermediate writes is safe -- the second and subsequent calls are no-ops.
     *
     * @note You cannot seek backwards in the file beyond the commit point, after calling commit().
     *
     * @param flush If true (default) the file is ready for reading when the method returns, otherwise the method returns
     *   immediately and writes the file in the background.
     *
     * @return True if the commit succeeded, else false.
     */
    bool commit(bool flush = true);

    /**
     * Copy data in a memory stream to a memory buffer (which must be preallocated to at least size() bytes). The stream buffer
     * remains unchanged.
     *
     * Triggers an error if this is a file stream.
     */
    void commit(uint8 * dst) const;

    /**
     * A memory stream may be reset so that it can be written to again without allocating new memory. The underlying array will
     * not be deallocated, but the reset structure will act like a newly intialized one.
     */
    void reset();

    /** Get the total number of bytes written. */
    int64 size() const
    {
      return m_bufferLen + m_alreadyWritten;
    }

    /**
     * Sets the length of the stream to \a n. If the new %size is greater than the old one, the content of the newly allocated
     * section is undefined. Does not change the position of the next byte to be written unless \a n < size(), in which case it
     * is set to the new stream end.
     *
     * Throws an error when resetting a huge file to be shorter than its current length.
     *
     * @note Currently, the increase in size cannot be more than the buffer capacity.
     */
    void setSize(int64 n)
    {
      n = n - m_alreadyWritten;

      if (n < 0)
        throw Error(getNameStr() + ": Cannot resize huge files to be shorter");
      else if (n < m_bufferLen)
      {
        m_bufferLen = n;
        m_pos = n;
      }
      else if (n > m_bufferLen)
        reserveBytes(n - m_pos);
    }

    /** Returns the current byte position in the file, where 0 is the beginning and size() - 1 is the end. */
    int64 getPosition() const
    {
      return m_pos + m_alreadyWritten;
    }

    /**
     * Sets the position. Can set past length, in which case the file is padded with undefined data up to one byte before the
     * next to be written.
     *
     * May throw an exception when seeking backwards in a huge file.
     *
     * @note Currently, the increase in size cannot be more than the buffer capacity.
     */
    void setPosition(int64 p)
    {
      int64 q = p - m_alreadyWritten;

      if (q < 0)
        throw Error(getNameStr() + ": Cannot seek too far backwards in a huge file");

      if (q > m_bufferLen)
        setSize(p);

      m_pos = p - m_alreadyWritten;
    }

    /**
     * Skips ahead \a n bytes (can skip past end-of-file, in which case the newly added bytes have undefined values). If \a n is
     * negative, seeks backwards from the current position
     */
    void skip(int64 n)
    {
      setPosition(m_pos + n);
    }

    /** Write a sequence of bytes. */
    void writeBytes(int64 n, void const * b)
    {
      reserveBytes(n);

      debugAssertM(m_pos >= 0, getNameStr() + ": Invalid write position");
      debugAssertM(m_bufferLen >= n, getNameStr() + format(": Could not reserve space to write %ld bytes", (long)n));

      std::memcpy(m_buffer + m_pos, b, n);
      m_pos += n;
    }

    /** Write an unsigned 8-bit integer. */
    void writeUInt8(uint8 i)
    {
      reserveBytes(1);
      m_buffer[m_pos] = i;
      m_pos++;
    }

    /** Write a signed 8-bit integer. */
    void writeInt8(int8 i)
    {
      reserveBytes(1);
      m_buffer[m_pos] = *(uint8 *)&i;
      m_pos++;
    }

    /** Write an 8-bit boolean value (0 for false, non-zero for true). */
    void writeBool8(bool b)
    {
      writeUInt8(b ? 1 : 0);
    }

    /** Write an unsigned 16-bit integer. */
    void writeUInt16(uint16 u);

    /** Write a signed 16-bit integer. */
    void writeInt16(int16 i)
    {
      writeUInt16(*(uint16 *)&i);
    }

    /** Write an unsigned 32-bit integer. */
    void writeUInt32(uint32 u);

    /** Write a signed 32-bit integer. */
    void writeInt32(int32 i)
    {
      writeUInt32(*(uint32 *)&i);
    }

    /** Write an unsigned 64-bit integer. */
    void writeUInt64(uint64 u);

    /** Write a signed 64-bit integer. */
    void writeInt64(int64 i)
    {
      writeUInt64(*(uint64 *)&i);
    }

    /** Write a 32-bit floating point number. */
    void writeFloat32(float32 f)
    {
      debugAssertM(m_beginEndBits == 0, getNameStr() + ": Byte-level writes not allowed in a beginBits/endBits block");

      union
      {
        float32 a;
        uint32 b;
      };
      a = f;
      writeUInt32(b);
    }

    /** Write a 64-bit floating point number. */
    void writeFloat64(float64 f)
    {
      debugAssertM(m_beginEndBits == 0, getNameStr() + ": Byte-level writes not allowed in a beginBits/endBits block");
      union
      {
        float64 a;
        uint64 b;
      };
      a = f;
      writeUInt64(b);
    }

    /**
     * Write a string. The format is:
     * - Length of string (32-bit integer)
     * - Characters of string ('length' bytes, no null termination)
     *
     * @see BinaryInputStream::readString
     *
     * @note This version explicitly writes the length of the string and does not rely on null termination, making this a safer
     * option than the original G3D version. To write a null-terminated string, use writeBytes.
     */
    void writeString(std::string const & s)
    {
      writeAlignedString(s, 1);
    }

    /**
     * Write an aligned string. The format is:
     * - Length of string (32-bit integer)
     * - Characters of string ('length' bytes, no null termination)
     * - Arbitrary padding bytes to align to \a alignment bytes boundary
     *
     * @see BinaryInputStream::readAlignedString()
     */
    void writeAlignedString(std::string const & s, int alignment = 4);

    /**
     * Begin writing individual bits. Call before a series of writeBits() calls. Only writeBits() can be called between
     * beginBits and endBits without corrupting the stream.
     */
    void beginBits();

    /**
     * Write \a num_bits bits from \a bit_string.  Bits are numbered from low to high.
     *
     * Can only be called between beginBits() and endBits().  Bits written are semantically little-endian, regardless of the
     * actual endianness of the system.  That is, <code>writeBits(0xABCD, 16)</code> writes 0xCD to the first byte and 0xAB to
     * the second byte.  However, if used with BinaryInput::readBits(), the ordering is transparent to the caller.
     */
    void writeBits(int num_bits, uint32 bit_string);

    /**
     * End writing individual bits. Call after a series of writeBits() calls. This will finish out with zeros the last byte into
     * which bits were written.
     */
    void endBits();

    /** Write a 2-vector. */
    void writeVector2(Vector2 const & v)
    {
      writeFloat32(v.x());
      writeFloat32(v.y());
    }

    /** Write a 3-vector. */
    void writeVector3(Vector3 const & v)
    {
      writeFloat32(v.x());
      writeFloat32(v.y());
      writeFloat32(v.z());
    }

    /** Write a 4-vector. */
    void writeVector4(Vector4 const & v)
    {
      writeFloat32(v.x());
      writeFloat32(v.y());
      writeFloat32(v.z());
      writeFloat32(v.w());
    }

    /** Write a color with 1 8-bit channel. */
    void writeColorL8(ColorL8 const & c)
    {
      writeUInt8(c.value());
    }

    /** Write a color with 1 floating-point channel. */
    void writeColorL(ColorL const & c)
    {
      writeFloat32(c.value());
    }

    /** Write a color with 3 8-bit channels. */
    void writeColorRGB8(ColorRGB8 const & c)
    {
      writeUInt8(c.r());
      writeUInt8(c.g());
      writeUInt8(c.b());
    }

    /** Write a color with 3 floating-point channels. */
    void writeColorRGB(ColorRGB const & c)
    {
      writeFloat32(c.r());
      writeFloat32(c.g());
      writeFloat32(c.b());
    }

    /** Write a color with 4 8-bit channels. */
    void writeColorRGBA8(ColorRGBA8 const & c)
    {
      writeUInt8(c.r());
      writeUInt8(c.g());
      writeUInt8(c.b());
      writeUInt8(c.a());
    }

    /** Write a color with 3 floating-point channels. */
    void writeColorRGBA(ColorRGBA const & c)
    {
      writeFloat32(c.r());
      writeFloat32(c.g());
      writeFloat32(c.b());
      writeFloat32(c.a());
    }

    /** Write a 2x2 matrix. */
    void writeMatrix2(Matrix2 const & m)
    {
      writeFloat32(m(0, 0));
      writeFloat32(m(0, 1));
      writeFloat32(m(1, 0));
      writeFloat32(m(1, 1));
    }

    /** Write a 3x3 matrix. */
    void writeMatrix3(Matrix3 const & m)
    {
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          writeFloat32(m(r, c));
    }

    /** Write a 4x4 matrix. */
    void writeMatrix4(Matrix4 const & m)
    {
      for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
          writeFloat32(m(r, c));
    }

    /** Write a coordinate frame. */
    void writeCoordinateFrame3(CoordinateFrame3 const & c)
    {
      writeMatrix3(c.getRotation());
      writeVector3(c.getTranslation());
    }

    /** Write a 3D plane. */
    void writePlane3(Plane3 const & plane)
    {
      Real a, b, c, d;
      plane.getEquation(a, b, c, d);
      writeFloat32(a);
      writeFloat32(b);
      writeFloat32(c);
      writeFloat32(d);
    }

#define THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(fname, tname)\
    void write##fname(int64 n, tname const * out); \
    void write##fname(int64 n, TheaArray<tname> const & out);

    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Bool8,               bool)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(UInt8,               uint8)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Int8,                int8)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(UInt16,              uint16)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Int16,               int16)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(UInt32,              uint32)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Int32,               int32)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(UInt64,              uint64)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Int64,               int64)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Float32,             float32)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Float64,             float64)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Vector2,             Vector2)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Vector3,             Vector3)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Vector4,             Vector4)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorL8,             ColorL8)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorL,              ColorL)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorRGB8,           ColorRGB8)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorRGB,            ColorRGB)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorRGBA8,          ColorRGBA8)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(ColorRGBA,           ColorRGBA)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Matrix2,             Matrix2)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Matrix3,             Matrix3)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Matrix4,             Matrix4)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(CoordinateFrame3,    CoordinateFrame3)
    THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER(Plane3,              Plane3)
#undef THEA_BINARY_OUTPUT_STREAM_DECLARE_WRITER

}; // class BinaryOutputStream

} // namespace Thea

#ifdef _MSC_VER
#   pragma warning (pop)
#endif

THEA_DECL_EXTERN_SMART_POINTERS(Thea::BinaryOutputStream)

#endif
