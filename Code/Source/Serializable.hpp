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
//============================================================================

#ifndef __Thea_Serializable_hpp__
#define __Thea_Serializable_hpp__

#include "Common.hpp"
#include "IOStream.hpp"
#include "NamedObject.hpp"

namespace Thea {

/**
 * A serialization codec. Identified by an ID that is unique for a given run of the program (it is <b>not</b> guaranteed to
 * retain its value over different runs).
 */
class THEA_API Codec : public AbstractNamedObject
{
  public:
    /** Destructor. */
    virtual ~Codec() {}

    /** Check if two codecs are equal. All instances of a codec class <b>must</b> be considered equal. */
    bool operator==(Codec const & other) const { return typeid(*this) == typeid(other); }

    /**
     * Implicitly convert to an integer value for use in switch statements etc. This value will be common to all instances of
     * the codec class
     */
    operator long() const { return reinterpret_cast<long>(&typeid(*this)); }
};

/** Write the name of the object to an output stream. */
inline std::ostream &
operator<<(std::ostream & os, Codec const & codec)
{
  return os << codec.getName() << " codec";
}

/** Indicates that the appropriate codec should be auto-detected. */
class THEA_API Codec_AUTO : public Codec
{
  public:
    char const * getName() const { static char const * my_name = "Auto"; return my_name; }
};

/** Indicates that the codec is unknown. */
class THEA_API Codec_UNKNOWN : public Codec
{
  public:
    char const * getName() const { static char const * my_name = "Unknown"; return my_name; }
};

/** The interface for a serializable object. */
class THEA_API Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(Serializable, shared_ptr, weak_ptr)

    /** Destructor. */
    virtual ~Serializable() {};

    /** Serialize the object to a binary output stream. */
    virtual void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const = 0;

    /** Deserialize the object from a binary input stream. */
    virtual void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO()) = 0;

    /** Serialize the object to a text output stream. */
    virtual void serialize(TextOutputStream & output, Codec const & codec = Codec_AUTO()) const
    { throw Error("Serialization to text stream not implemented"); }

    /** Deserialize the object from a text input stream. */
    virtual void deserialize(TextInputStream & input, Codec const & codec = Codec_AUTO())
    { throw Error("Deserialization from text stream not implemented"); }

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
    THEA_DEF_POINTER_TYPES(SerializableFactory, shared_ptr, weak_ptr)

    /** Destructor. */
    virtual ~SerializableFactory() {}

    /** Create an instance of the serializable class with a given name. */
    virtual Serializable::Ptr createSerializable(std::string const & name) const = 0;
};

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Serializable)
THEA_DECL_EXTERN_SMART_POINTERS(Thea::SerializableFactory)

#endif
