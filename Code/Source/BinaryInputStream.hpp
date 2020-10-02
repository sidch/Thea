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
// First version: 2013
//
//============================================================================

/*
 ORIGINAL HEADER

 @file BinaryInput.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2001-08-09
 @edited  2010-03-19

 Copyright 2000-2010, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_BinaryInputStream_hpp__
#define __Thea_BinaryInputStream_hpp__

#ifdef _MSC_VER
// Disable conditional expression is constant, which occurs incorrectly on inlined functions
#  pragma warning(push)
#  pragma warning( disable : 4127 )
#endif

#include "Common.hpp"
#include "Array.hpp"
#include "Codec.hpp"
#include "Colors.hpp"
#include "CoordinateFrame3.hpp"
#include "MatVec.hpp"
#include "NamedObject.hpp"
#include "Noncopyable.hpp"
#include "Plane3.hpp"

namespace Thea {

/**
 * Sequential or random access input from binary files/memory. For every <code>readX</code> method there are also versions that
 * operate on a whole Array or C-array. The first method resizes the Array to the appropriate size before reading. For a
 * C-array, they require the pointer to reference a memory block large enough to hold <I>n</I> elements.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 *
 * @todo Reimplement using %<iostream%> for arbitrary seeking and safer performance?
 */
class THEA_API BinaryInputStream : public virtual NamedObject, private Noncopyable
{
  private:
    /** Is the file big or little endian? */
    Endianness         m_fileEndian;

    /** Path to the opened file. */
    std::string        m_path;

    /** Swap bytes when writing? */
    bool               m_swapBytes;

    /** Next position to read from in bitString during readBits. */
    int                m_bitPos;

    /**
     * Bits currently being read by readBits. Contains at most 8 (low) bits.  Note that beginBits/readBits actually consumes one
     * extra byte, which will be restored by writeBits.*/
    uint32             m_bitString;

    /** 1 when between beginBits and endBits, 0 otherwise. */
    int                m_beginEndBits;

    /**
     * When operating on huge files, we cannot load the whole file into memory. This is the file position to which buffer[0]
     * corresponds.
     */
    int64              m_alreadyRead;

    /** Length of the entire file, in bytes. For the length of the buffer, see m_bufferLength */
    int64              m_length;

    /** Length of the array referenced by buffer. May go past the end of the file! */
    int64              m_bufferLength;

    /** The buffer of bytes from the input stream. */
    uint8            * m_buffer;

    /** Next byte in file, relative to buffer. */
    int64              m_pos;

    /** When true, the buffer is freed in the destructor. */
    bool               m_freeBuffer;

    /** Ensures that we are able to read at least min_length from start_position (relative to start of file). */
    void loadIntoMemory(int64 start_position, int64 min_length = 0);

    /** Verifies that at least this number of bytes can be read.*/
    void prepareToRead(int64 nbytes)
    {
      if (m_pos + nbytes + m_alreadyRead > m_length)
        throw Error(getNameStr() + ": Read past end of stream");

      if (m_pos + nbytes > m_bufferLength)
        loadIntoMemory(m_pos + m_alreadyRead, nbytes);
    }

    // Not implemented on purpose, don't use.
    bool operator==(BinaryInputStream const &) const;

  public:
    THEA_DECL_SMART_POINTERS(BinaryInputStream)

    /** Constant to use with the copy_memory option (evaluates to false). */
    static bool const NO_COPY;

    /**
     * An object that saves the current endianness state of a stream upon construction, and restores the previous state upon
     * destruction, using the RAII idiom.
     */
    class EndiannessScope
    {
      public:
        /** Constructor, saves the endianness state of a stream. */
        EndiannessScope(BinaryInputStream & stream_)
        : stream(stream_), saved_endian(stream_.getEndianness())
        {}

        /** Constructor, saves the endianness state of a stream and sets a new endianness. */
        EndiannessScope(BinaryInputStream & stream_, Endianness new_endian)
        : stream(stream_), saved_endian(stream_.getEndianness())
        {
          stream.setEndianness(new_endian);
        }

        /** Destructor, restores the saved endianness of the stream. */
        ~EndiannessScope() { stream.setEndianness(saved_endian); }

      private:
        BinaryInputStream & stream;  ///< The wrapped stream.
        Endianness saved_endian;     ///< The stream's saved endianness.

    }; // class EndiannessScope

    /** Open a file as a binary input stream. If the file cannot be opened, an error is thrown. */
    BinaryInputStream(std::string const & path, Endianness file_endian);

    /**
     * Wrap a block of in-memory data as an input stream. Unless you specify \a copy_memory = false, the data is copied from the
     * pointer, so you may deallocate it as soon as the object is constructed. It is an error to specify \a copy_memory = false.
     */
    BinaryInputStream(uint8 const * data, int64 data_len, Endianness data_endian, bool copy_memory = true);

    /**
     * Wrap the next \a block_len bytes of another input stream as a new, temporary input stream with the same endianness. The
     * source stream must exist as long as this object does, since the implementation does not guarantee that the block will be
     * fully copied to a new buffer. The block is marked as read in the source stream, and the next read position in that stream
     * is set to just after the block.
     */
    BinaryInputStream(BinaryInputStream & src, int64 block_len);

    /** Destructor. */
    ~BinaryInputStream();

    /**
     * Change the endianness for interpreting the file contents. This only changes the interpretation of the file for future
     * read calls; the underlying data is unmodified. To temporarily set the endianness within a scope using the RAII idiom,
     * and restore the original setting at the end of the scope, create an EndiannessScope object.
     */
    void setEndianness(Endianness endian);

    /** Get the current endianness setting for reading multi-byte data. */
    Endianness getEndianness() const
    {
      return m_fileEndian;
    }

    /** Get the path to the current file being read ("<memory>" for memory streams). */
    std::string getPath() const
    {
      return m_path;
    }

    /** Get the number of bytes in the stream. */
    int64 size() const
    {
      return m_length;
    }

    /** Check if there are more bytes to be read. Returns true if the current read position is not at the end of the file. */
    bool hasMore() const
    {
      return m_pos + m_alreadyRead < m_length;
    }

    /** Get the current byte position in the stream, where 0 is the beginning and size() - 1 is the end. */
    int64 getPosition() const
    {
      return m_pos + m_alreadyRead;
    }

    /**
     * Set the current byte position in the stream. May throw an assertion failure when seeking backwards more than 10MB on a
     * huge file.
     */
    void setPosition(int64 p)
    {
      if (p > m_length)
        throw Error("Read past end of stream");

      m_pos = p - m_alreadyRead;

      if ((m_pos < 0) || (m_pos > m_bufferLength))
        loadIntoMemory(m_pos + m_alreadyRead);
    }

    /** Set the read position at the beginning of the file. */
    void reset()
    {
      setPosition(0);
    }

    /** Skips ahead \a n bytes. */
    void skip(int64 n)
    {
      setPosition(m_pos + m_alreadyRead + n);
    }

    /** Read an 8-bit unsigned integer (an unsigned byte). */
    uint8 readUInt8()
    {
      prepareToRead(1);
      return ((uint8 *)m_buffer)[m_pos++];
    }

    /** Read an 8-bit signed integer (a signed byte). */
    int8 readInt8()
    {
      prepareToRead(1);
      return m_buffer[m_pos++];
    }

    /** Read an 8-bit boolean (0 is false, non-zero is true). */
    bool readBool8()
    {
      return (readInt8() != 0);
    }

    /** Read a 16-bit unsigned integer. */
    uint16 readUInt16()
    {
      prepareToRead(2);
      m_pos += 2;

      if (m_swapBytes)
      {
        uint8 out[2];
        out[0] = m_buffer[m_pos - 1];
        out[1] = m_buffer[m_pos - 2];
        return *(uint16 *)out;
      }
      else
      {
#ifdef THEA_ALLOW_UNALIGNED_READS
        return *(uint16 *)(&m_buffer[m_pos - 2]);
#else
        uint8 out[2];
        out[0] = m_buffer[m_pos - 2];
        out[1] = m_buffer[m_pos - 1];
        return *(uint16 *)out;
#endif
      }
    }

    /** Read a 16-bit signed integer. */
    int16 readInt16()
    {
      uint16 a = readUInt16();
      return *(int16 *)&a;
    }

    /** Read a 32-bit unsigned integer. */
    uint32 readUInt32()
    {
      prepareToRead(4);
      m_pos += 4;

      if (m_swapBytes)
      {
        uint8 out[4];
        out[0] = m_buffer[m_pos - 1];
        out[1] = m_buffer[m_pos - 2];
        out[2] = m_buffer[m_pos - 3];
        out[3] = m_buffer[m_pos - 4];
        return *(uint32 *)out;
      }
      else
      {
#ifdef THEA_ALLOW_UNALIGNED_READS
        return *(uint32 *)(&m_buffer[m_pos - 4]);
#else
        uint8 out[4];
        out[0] = m_buffer[m_pos - 4];
        out[1] = m_buffer[m_pos - 3];
        out[2] = m_buffer[m_pos - 2];
        out[3] = m_buffer[m_pos - 1];
        return *(uint32 *)out;
#endif
      }
    }

    /** Read a 32-bit signed integer. */
    int32 readInt32()
    {
      uint32 a = readUInt32();
      return *(int32 *)&a;
    }

    /** Read a 64-bit unsigned integer. */
    uint64 readUInt64();

    /** Read a 64-bit signed integer. */
    int64 readInt64()
    {
      uint64 a = readUInt64();
      return *(int64 *)&a;
    }

    /** Read a 32-bit floating point number. */
    float32 readFloat32()
    {
      union
      {
        uint32 a;
        float32 b;
      };

      a = readUInt32();
      return b;
    }

    /** Read a 64-bit floating point number. */
    float64 readFloat64()
    {
      union
      {
        uint64 a;
        float64 b;
      };

      a = readUInt64();
      return b;
    }

    /** Read a sequence of \a n bytes. */
    void readBytes(int64 n, void * bytes);

    /**
     * Reads until any newline character (\\r, \\r\\n, \\n\\r, \\n) or the end of the file is encountered. Consumes the newline.
     */
    std::string readLine();

    /**
     * Read a string. The format is:
     * - Length of string (32-bit integer)
     * - Characters of string ('length' bytes, no null termination)
     *
     * @see BinaryOutputStream::writeString()
     */
    std::string readString()
    {
      return readAlignedString(1);
    }

    /**
     * Read an aligned string. The format is:
     * - Length of string (32-bit integer)
     * - Characters of string ('length' bytes, no null termination)
     * - Arbitrary padding bytes to align to \a alignment bytes boundary
     *
     * @see BinaryOutputStream::writeAlignedString()
     */
    std::string readAlignedString(int alignment = 4);

    /**
     * Read an \a n character string.  The string is not required to end in nullptr in the file but will always be a proper
     * std::string when returned.
     */
    std::string readString(int64 n);

    /**
     * Read a null-terminated string. This is an <b>INHERENTLY UNSAFE</b> method and should be used with <b>GREAT CAUTION</b>!
     */
    std::string readNullTerminatedString();

    /**
     * Prepare for reading individual bits via readBits(). Only readBits() can be called between beginBits() and endBits()
     * without corrupting the data stream.
     */
    void beginBits();

    /** Read a sequence of up to 32 bits. Can only be called between beginBits() and endBits() */
    uint32 readBits(int num_bits);

    /** End bit-reading. */
    void endBits();

    /** Read a 2-vector. */
    Vector2 readVector2()
    {
      Vector2 v;
      v[0] = readFloat32();
      v[1] = readFloat32();
      return v;
    }

    /** Read a 3-vector. */
    Vector3 readVector3()
    {
      Vector3 v;
      v[0] = readFloat32();
      v[1] = readFloat32();
      v[2] = readFloat32();
      return v;
    }

    /** Read a 4-vector. */
    Vector4 readVector4()
    {
      Vector4 v;
      v[0] = readFloat32();
      v[1] = readFloat32();
      v[2] = readFloat32();
      v[3] = readFloat32();
      return v;
    }

    /** Read a color with 1 8-bit channel. */
    ColorL8 readColorL8()
    {
      return ColorL8(readUInt8());
    }

    /** Read a color with 1 floating-point channel. */
    ColorL readColorL()
    {
      return ColorL(readFloat32());
    }

    /** Read a color with 3 8-bit channels. */
    ColorRgb8 readColorRgb8()
    {
      ColorRgb8 c;
      c.r() = readUInt8();
      c.g() = readUInt8();
      c.b() = readUInt8();
      return c;
    }

    /** Read a color with 3 floating-point channels. */
    ColorRgb readColorRgb()
    {
      ColorRgb c;
      c.r() = readFloat32();
      c.g() = readFloat32();
      c.b() = readFloat32();
      return c;
    }

    /** Read a color with 4 8-bit channels. */
    ColorRgba8 readColorRgba8()
    {
      ColorRgba8 c;
      c.r() = readUInt8();
      c.g() = readUInt8();
      c.b() = readUInt8();
      c.a() = readUInt8();
      return c;
    }

    /** Read a color with 4 floating-point channels. */
    ColorRgba readColorRgba()
    {
      ColorRgba c;
      c.r() = readFloat32();
      c.g() = readFloat32();
      c.b() = readFloat32();
      c.a() = readFloat32();
      return c;
    }

    /** Read a 2x2 matrix. */
    Matrix2 readMatrix2()
    {
      Matrix2 m;
      m(0, 0) = readFloat32();
      m(0, 1) = readFloat32();
      m(1, 0) = readFloat32();
      m(1, 1) = readFloat32();
      return m;
    }

    /** Read a 3x3 matrix. */
    Matrix3 readMatrix3()
    {
      Matrix3 m;
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          m(r, c) = readFloat32();

      return m;
    }

    /** Read a 4x4 matrix. */
    Matrix4 readMatrix4()
    {
      Matrix4 m;
      for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
          m(r, c) = readFloat32();

      return m;
    }

    /** Read a coordinate frame. */
    CoordinateFrame3 readCoordinateFrame3()
    {
      Matrix3 rot = readMatrix3();
      Vector3 trn = readVector3();
      return CoordinateFrame3(RigidTransform3::_fromAffine(AffineTransform3(rot, trn)));
    }

    /** Read a 3D plane. */
    Plane3 readPlane3()
    {
      Real a = readFloat32();
      Real b = readFloat32();
      Real c = readFloat32();
      Real d = readFloat32();
      return Plane3::fromEquation(a, b, c, d);
    }

    /**
     * Read an arbitrary dense or sparse 2D matrix, using a codec such as CSV or HDF5 (some of which can be autodetected).
     *
     * @param read_block_header If true, first read a header section which stores the size and codec of the serialized matrix
     * @param m The matrix in which to store the deserialized input.
     * @param codec The codec to use (pass CodecAuto() to autodetect it from the input).
     *   data. Else, the matrix block is assumed to continue until the end of the input stream unless its end can be detected
     *   through some other means (e.g. end marker or embedded size field), and the codec will be autodetected if possible.
     *
     * @return Always returns 0 (except on error when it throws an exception), for consistency with the other version of this
     *   function which chooses between a variable number of matrix types depending on the input.
     *
     * @note Unlike most other deserialization functions that optionally read a block header, the \a read_block_header parameter
     *   is the first parameter instead of the last one for consistency with the other version of this function which chooses
     *   between a variable number of matrix types depending on the input.
     */
    template <typename MatrixT>
    intx readMatrix(bool read_block_header, MatrixT & m, Codec const & codec = CodecAuto());

    /**
     * Given a set of candidate matrices of different types, read that one whose format best matches the input data. The codec
     * will be auto-detected. The most common use is when either a dense matrix or a sparse matrix must be read, but it is not
     * known ahead of time which is present in the input. Sample use in such a situation could be something like:
     *
     * \code
     * MatrixXd dense;
     * SparseMatrix<double> sparse;
     * auto which = input.readMatrix(true, dense, sparse);
     * THEA_CONSOLE << "Read " << (which == 0 ? "dense" : "sparse") << " matrix";
     * \endcode
     *
     * @param read_block_header If true, first read a header section which stores the size and codec of the serialized matrix.
     * @param m0 A matrix of the first candidate type.
     * @param m1 A matrix of the second candidate type.
     * @param rest Optional matrices of more candidate types.
     *
     * @return The index of the matrix which was read. The first matrix (of type CandidateMatrixType0) has index 0, the second
     *   (of type CandidateMatrixType1) has index 1, and so on.
     */
    template < typename CandidateMatrixType0, typename CandidateMatrixType1, typename... MoreMatrixTypes,
               typename std::enable_if< !std::is_base_of<Codec, CandidateMatrixType1>::value, int >::type * = nullptr >
    intx readMatrix(bool read_block_header, CandidateMatrixType0 & m0, CandidateMatrixType1 & m1, MoreMatrixTypes & ... rest)
    {
      return readMatrixHelper(0, read_block_header, m0, m1, rest...);
    }

#   define THEA_BINARY_INPUT_STREAM_DECLARE_READER(fname, tname) \
    void read##fname(int64 n, tname * out); \
    void read##fname(int64 n, Array<tname> & out);

    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Bool8,               bool)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(UInt8,               uint8)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Int8,                int8)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(UInt16,              uint16)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Int16,               int16)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(UInt32,              uint32)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Int32,               int32)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(UInt64,              uint64)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Int64,               int64)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Float32,             float32)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Float64,             float64)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Vector2,             Vector2)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Vector3,             Vector3)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Vector4,             Vector4)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorL8,             ColorL8)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorL,              ColorL)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorRgb8,           ColorRgb8)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorRgb,            ColorRgb)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorRgba8,          ColorRgba8)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(ColorRgba,           ColorRgba)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Matrix2,             Matrix2)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Matrix3,             Matrix3)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Matrix4,             Matrix4)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(CoordinateFrame3,    CoordinateFrame3)
    THEA_BINARY_INPUT_STREAM_DECLARE_READER(Plane3,              Plane3)
#   undef THEA_BINARY_INPUT_STREAM_DECLARE_READER

  private:
    /**
     * Base case for reading a matrix of the type best matching the input. This base case takes no matrix arguments and always
     * throws an error.
     */
    intx readMatrixHelper(intx index, bool read_block_header);

    /** Helper function for reading a matrix of the type best matching the input. */
    template <typename CandidateMatrixT, typename... MoreMatrixTypes>
    intx readMatrixHelper(intx index, bool read_block_header, CandidateMatrixT & m, MoreMatrixTypes & ... rest);

}; // class BinaryInputStream

} // namespace Thea

#ifdef _MSC_VER
#   pragma  warning(pop)
#endif

THEA_DECL_EXTERN_SMART_POINTERS(Thea::BinaryInputStream)

#include "MatrixIO.hpp"

#endif
