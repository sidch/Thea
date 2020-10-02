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

#ifndef __Thea_Serializable_hpp__
#define __Thea_Serializable_hpp__

#include "Common.hpp"
#include "Codec.hpp"
#include "Iostream.hpp"

namespace Thea {

/** Abstract base class for a serializable object. */
class THEA_API Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Serializable)

    /** Destructor. */
    virtual ~Serializable() {};

    /**
     * Read the object from a binary input stream.
     *
     * @param input The stream from which to read data.
     * @param codec The codec to use. If set to CodecAuto(), the codec will be autodetected (if possible) from the input.
     * @param read_block_header If true, a Codec::BlockHeader object containing information about the codec and size of the
     *   serialized block will be first read from the input, and used to aid codec detection etc. The implementation is <b>free
     *   to ignore this directive</b>, e.g. if the codec and input size can be detected through other means. This behavior must
     *   be synchronized with write(BinaryOutputStream &, Codec const &, bool): either both must omit block headers, or both
     *   must read/write them if directed to do so.
     */
    virtual void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false) = 0;

    /**
     * Write the object to a binary output stream.
     *
     * @param output The stream to which data will be written.
     * @param codec The codec to use. If set to CodecAuto(), an appropriate codec will be automatically selected.
     * @param write_block_header If true, a Codec::BlockHeader object containing information about the codec and size of the
     *   serialized block will be first written to the input. The implementation is <b>free to ignore this directive</b>, e.g.
     *   if the codec and input size are encoded through other means. This behavior must be synchronized with
     *   read(BinaryInputStream &, Codec const &, bool): either both must omit block headers, or both must read/write them if
     *   directed to do so.
     */
    virtual void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false)
                 const = 0;

    /** Read the object from a text input stream. */
    virtual void read(TextInputStream & input, Codec const & codec = CodecAuto())
    { throw Error("Deserialization from text stream not implemented"); }

    /** Write the object to a text output stream. */
    virtual void write(TextOutputStream & output, Codec const & codec = CodecAuto()) const
    { throw Error("Serialization to text stream not implemented"); }

    /** Get the default settings for parsing configuration text files. */
    static TextInputStream::Settings const & configReadSettings()
    {
      static const TextInputStream::Settings def = initConfigReadSettings();
      return def;
    }

    /** Get the default settings for writing configuration text files. */
    static TextOutputStream::Settings const & configWriteSettings()
    {
      static const TextOutputStream::Settings def = initConfigWriteSettings();
      return def;
    }

  private:
    /** Create default settings for parsing configuration text files. */
    static TextInputStream::Settings initConfigReadSettings();

    /** Create default settings for writing configuration text files. */
    static TextOutputStream::Settings initConfigWriteSettings();

}; // class Serializable

/**
 * The interface for a factory for creating serializable objects. Classes which provide factories should define a
 * <tt>getFactory</tt> function. Factories are useful for loading generic objects from a database, where the only parameters
 * that need to be passed to the database are the name of the object and a factory (accessed through a pointer of this base
 * class type) for the particular type of object to be loaded. This way, the database can load virtually any serializable object
 * without needing specialization for each class of object.
 */
class THEA_API SerializableFactory
{
  public:
    THEA_DECL_SMART_POINTERS(SerializableFactory)

    /** Destructor. */
    virtual ~SerializableFactory() {}

    /** Create an instance of the serializable class with a given name. */
    virtual Serializable * createSerializable(char const * name) const = 0;

    /** Destroy a serializable object created with createSerializable(). */
    virtual void destroySerializable(Serializable * serializable) = 0;
};

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Serializable)
THEA_DECL_EXTERN_SMART_POINTERS(Thea::SerializableFactory)

#endif
