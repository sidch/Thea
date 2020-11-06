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

#ifndef __Thea_Codec_hpp__
#define __Thea_Codec_hpp__

#include "Common.hpp"
#include "NamedObject.hpp"
#include <algorithm>
#include <array>

namespace Thea {

// Forward declarations
class BinaryInputStream;
class BinaryOutputStream;

/**
 * A serialization codec. Identified by an ID that is unique for a given run of the program (it is <b>not</b> guaranteed to
 * retain its value over different runs).
 */
class THEA_API Codec : public virtual INamedObject
{
  public:
    THEA_DECL_SMART_POINTERS(Codec)

    /** The standard length (in bytes) of the codec's magic string as used in BlockHeader. */
    static intx const MAGIC_LENGTH = 8;

    /** The type of the codec's magic string (an array of MAGIC_LENGTH bytes). */
    typedef std::array<int8, (size_t)MAGIC_LENGTH> MagicString;

    /**
     * A header preceding a data block serialized using a codec. The header contains the size of the data block, a magic string
     * identifying the serialization codec, and an optional 64-bit field that can be used for additional information.
     */
    struct BlockHeader
    {
      public:
        /** The length in bytes of a serialized block header. */
        static intx const SERIALIZED_LENGTH = MAGIC_LENGTH
                                            + /* sizeof(data_size) = */ 8
                                            + /* sizeof(custom) = */ 8
                                            + /* reserved */ 8;

        MagicString  magic;      ///< A magic string identifying the codec used to serialize the data.
        uint64       data_size;  ///< The size of the data block in bytes.
        uint64       custom;     ///< Additional custom data that can be optionally stored in the header.

        /** Construct with a magic string and a data block size. */
        BlockHeader(MagicString const & magic_ = zeroMagic(), uint64 data_size_ = 0)
        : magic(magic_), data_size(data_size_), custom(0), header_pos(0) {}

        /** Construct with a magic string specified as a std::string, and a data block size. */
        BlockHeader(std::string const & magic_, uint64 data_size_ = 0)
        : magic(toMagic(magic_)), data_size(data_size_), custom(0), header_pos(0) {}

        /** Construct by calling read() on an input stream. */
        BlockHeader(BinaryInputStream & in);

        /** Read the header from an input stream. This is guaranteed to read exactly SERIALIZED_LENGTH bytes. */
        void read(BinaryInputStream & in);

        /** Write the header to an output stream. This is guaranteed to write exactly SERIALIZED_LENGTH bytes. */
        void write(BinaryOutputStream & out) const;

        /**
         * Internally save the current location in an output stream, while reserving room to write a block header at that
         * position. This function should be called before writing a data block, and the return value passed to calcAndWrite()
         * after the data block is written.
         *
         * \code
         *   header.markAndSkip(out);
         *   // write data block
         *   header.calcAndWrite(out);
         * \endcode
         *
         * @return The saved location.
         *
         * @see calcAndWrite()
         */
        int64 markAndSkip(BinaryOutputStream & out);

        /**
         * Write the header to an output stream, after calculating the data block size based on where the header was supposed to
         * have been written (calculated using markAndSkip()) vs the current stream position, assumed to be at the end of the
         * data block. After writing the header, the next write position is moved back to the end of the data block. This is
         * guaranteed to write exactly Codec::BlockHeader::SERIALIZED_LENGTH bytes starting at \a header_pos.
         *
         * \code
         *   header.markAndSkip(out);
         *   // write data block
         *   header.calcAndWrite(out);
         * \endcode
         *
         * @note This function updates data_size.
         *
         * @warning This function seeks backwards in the stream, which may fail if a very large file is being written. This is
         *   because of the current implementation of BinaryOutputStream::setPosition(), which should be fixed at some point.
         *
         * @see markAndSkip()
         */
        void calcAndWrite(BinaryOutputStream & out);

      private:
        int64 header_pos;  ///< Used to store the marked stream position at which the header should be written.

    }; // struct BlockHeader

    /** Destructor. */
    virtual ~Codec() = 0;

    int8 THEA_ICALL setName(char const * s) { return false; /* codec name is read-only by default */ }

    /** Get the magic string for the codec, if it has one (else a string of all zeros). */
    virtual MagicString const & getMagic() const { static MagicString const ZERO = zeroMagic(); return ZERO; }

    /** Check if two codecs are equal. All instances of a codec class <b>must</b> be considered equal. */
    bool operator==(Codec const & other) const { return typeid(*this) == typeid(other); }

    /**
     * Implicitly convert to an integer value for use in switch statements etc. This value will be common to all instances of
     * the codec class
     */
    operator intx() const { return reinterpret_cast<intx>(&typeid(*this)); }

    /** Convenience function to convert a string literal to a magic string. */
    static MagicString toMagic(std::string const & s)
    {
      MagicString m;
      m.fill(0);  // padding, if needed
      std::copy_n(s.data(), std::min((size_t)MAGIC_LENGTH, s.length()), m.data());
      return m;
    }

  private:
    /** Construct a magic string with all zeros. */
    static MagicString zeroMagic() { MagicString z; z.fill(0); return z; }

}; // class Codec

// Pure virtual destructor should have a body
inline Codec::~Codec() {}

/** Write the name of the object to an output stream. */
inline std::ostream &
operator<<(std::ostream & os, Codec const & codec)
{
  return os << codec.getName() << " codec";
}

/** Indicates that the appropriate codec should be autodetected. */
class THEA_API CodecAuto : public Codec
{
  public:
    char const * THEA_ICALL getName() const { static char const * my_name = "Auto"; return my_name; }
};

/** Indicates that the codec is unknown. */
class THEA_API Codec_UNKNOWN : public Codec
{
  public:
    char const * THEA_ICALL getName() const { static char const * my_name = "Unknown"; return my_name; }
};

} // namespace Thea

#endif
