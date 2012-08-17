//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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
#include "CoordinateFrameN.hpp"
#include "HyperplaneN.hpp"
#include "IOStream.hpp"
#include "MatrixMN.hpp"
#include "VectorN.hpp"

namespace Thea {

/**
 * A serialization codec. Identified by an ID that is unique for a given run of the program (it is <b>not</b> guaranteed to
 * retain its value over different runs).
 */
class THEA_API Codec
{
  public:
    /** Destructor. */
    virtual ~Codec() {}

    /** Get the name of the codec. */
    virtual std::string getName() const = 0;

    /** Check if two codecs are equal. All instances of a codec class <b>must</b> be considered equal. */
    bool operator==(Codec const & other) const { return typeid(*this) == typeid(other); }

    /**
     * Implicitly convert to an integer value for use in switch statements etc. This value will be common to all instances of
     * the codec class
     */
    operator long() const { return reinterpret_cast<long>(&typeid(*this)); }
};

/** Write the name of the object to an output stream. */
inline THEA_API std::ostream &
operator<<(std::ostream & os, Codec const & codec)
{
  return os << codec.getName() << " codec";
}

/** Indicates that the appropriate codec should be auto-detected. */
class THEA_API Codec_AUTO : public Codec
{
  public:
    std::string getName() const { static std::string const my_name = "Auto"; return my_name; }
};

/** Indicates that the codec is unknown. */
class THEA_API Codec_UNKNOWN : public Codec
{
  public:
    std::string getName() const { static std::string const my_name = "Unknown"; return my_name; }
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
};

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

// Serialization/deserialization functions for G3D-derived classes. All the functions have suffixes indicating the type of
// argument, to prevent user errors caused by implicit conversions.

/** Serialize a 2-vector to a binary output stream. */
inline THEA_API void
serializeVector2(Vector2 const & v, BinaryOutputStream & output)
{
  output.writeFloat32(v.x());
  output.writeFloat32(v.y());
}

/** Deserialize a 2-vector from a binary input stream. */
inline THEA_API void
deserializeVector2(Vector2 & v, BinaryInputStream & input)
{
  v.x() = input.readFloat32();
  v.y() = input.readFloat32();
}

/** Serialize a 3-vector to a binary output stream. */
inline THEA_API void
serializeVector3(Vector3 const & v, BinaryOutputStream & output)
{
  output.writeFloat32(v.x());
  output.writeFloat32(v.y());
  output.writeFloat32(v.z());
}

/** Deserialize a 3-vector from a binary input stream. */
inline THEA_API void
deserializeVector3(Vector3 & v, BinaryInputStream & input)
{
  v.x() = input.readFloat32();
  v.y() = input.readFloat32();
  v.z() = input.readFloat32();
}

/** Serialize a 4-vector to a binary output stream. */
inline THEA_API void
serializeVector4(Vector4 const & v, BinaryOutputStream & output)
{
  output.writeFloat32(v.x());
  output.writeFloat32(v.y());
  output.writeFloat32(v.z());
  output.writeFloat32(v.w());
}

/** Deserialize a 4-vector from a binary input stream. */
inline THEA_API void
deserializeVector4(Vector4 & v, BinaryInputStream & input)
{
  v.x() = input.readFloat32();
  v.y() = input.readFloat32();
  v.z() = input.readFloat32();
  v.w() = input.readFloat32();
}

/** Serialize a color with 1 8-bit channel to a binary output stream. */
inline THEA_API void
serializeColor1uint8(Color1uint8 const & c, BinaryOutputStream & output)
{
  output.writeUInt8(c.value);
}

/** Deserialize a color with 1 8-bit channel from a binary input stream. */
inline THEA_API void
deserializeColor1uint8(Color1uint8 & c, BinaryInputStream & input)
{
  c.value = input.readUInt8();
}

/** Serialize a color with 1 floating-point channel to a binary output stream. */
inline THEA_API void
serializeColor1(Color1 const & c, BinaryOutputStream & output)
{
  output.writeFloat32(c.value);
}

/** Deserialize a color with 1 floating-point channel from a binary input stream. */
inline THEA_API void
deserializeColor1(Color1 & c, BinaryInputStream & input)
{
  c.value = input.readFloat32();
}

/** Serialize a color with 3 8-bit channels to a binary output stream. */
inline THEA_API void
serializeColor3uint8(Color3uint8 const & c, BinaryOutputStream & output)
{
  output.writeUInt8(c.r);
  output.writeUInt8(c.g);
  output.writeUInt8(c.b);
}

/** Deserialize a color with 3 8-bit channels from a binary input stream. */
inline THEA_API void
deserializeColor3uint8(Color3uint8 & c, BinaryInputStream & input)
{
  c.r = input.readUInt8();
  c.g = input.readUInt8();
  c.b = input.readUInt8();
}

/** Serialize a color with 3 floating-point channels to a binary output stream. */
inline THEA_API void
serializeColor3(Color3 const & c, BinaryOutputStream & output)
{
  output.writeFloat32(c.r);
  output.writeFloat32(c.g);
  output.writeFloat32(c.b);
}

/** Deserialize a color with 3 floating-point channels from a binary input stream. */
inline THEA_API void
deserializeColor3(Color3 & c, BinaryInputStream & input)
{
  c.r = input.readFloat32();
  c.g = input.readFloat32();
  c.b = input.readFloat32();
}

/** Serialize a color with 4 8-bit channels to a binary output stream. */
inline THEA_API void
serializeColor4uint8(Color4uint8 const & c, BinaryOutputStream & output)
{
  output.writeUInt8(c.r);
  output.writeUInt8(c.g);
  output.writeUInt8(c.b);
  output.writeUInt8(c.a);
}

/** Deserialize a color with 4 8-bit channels from a binary input stream. */
inline THEA_API void
deserializeColor4uint8(Color4uint8 & c, BinaryInputStream & input)
{
  c.r = input.readUInt8();
  c.g = input.readUInt8();
  c.b = input.readUInt8();
  c.a = input.readUInt8();
}

/** Serialize a color with 3 floating-point channels to a binary output stream. */
inline THEA_API void
serializeColor4(Color4 const & c, BinaryOutputStream & output)
{
  output.writeFloat32(c.r);
  output.writeFloat32(c.g);
  output.writeFloat32(c.b);
  output.writeFloat32(c.a);
}

/** Deserialize a color with 4 floating-point channels from a binary input stream. */
inline THEA_API void
deserializeColor4(Color4 & c, BinaryInputStream & input)
{
  c.r = input.readFloat32();
  c.g = input.readFloat32();
  c.b = input.readFloat32();
  c.a = input.readFloat32();
}

/** Serialize a 2x2 matrix to a binary output stream. */
inline THEA_API void
serializeMatrix2(Matrix2 const & m, BinaryOutputStream & output)
{
  output.writeFloat32(m(0, 0));
  output.writeFloat32(m(0, 1));
  output.writeFloat32(m(1, 0));
  output.writeFloat32(m(1, 1));
}

/** Deserialize a 2x2 matrix from a binary input stream. */
inline THEA_API void
deserializeMatrix2(Matrix2 & m, BinaryInputStream & input)
{
  m(0, 0) = input.readFloat32();
  m(0, 1) = input.readFloat32();
  m(1, 0) = input.readFloat32();
  m(1, 1) = input.readFloat32();
}

/** Serialize a 3x3 matrix to a binary output stream. */
inline THEA_API void
serializeMatrix3(Matrix3 const & m, BinaryOutputStream & output)
{
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      output.writeFloat32(m(r, c));
}

/** Deserialize a 3x3 matrix from a binary input stream. */
inline THEA_API void
deserializeMatrix3(Matrix3 & m, BinaryInputStream & input)
{
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      m(r, c) = input.readFloat32();
}

/** Serialize a 4x4 matrix to a binary output stream. */
inline THEA_API void
serializeMatrix4(Matrix4 const & m, BinaryOutputStream & output)
{
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c)
      output.writeFloat32(m(r, c));
}

/** Deserialize a 4x4 matrix from a binary input stream. */
inline THEA_API void
deserializeMatrix4(Matrix4 & m, BinaryInputStream & input)
{
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c)
      m(r, c) = input.readFloat32();
}

/** Serialize a coordinate frame to a binary output stream. */
inline THEA_API void
serializeCoordinateFrame3(CoordinateFrame3 const & c, BinaryOutputStream & output)
{
  serializeMatrix3(c.getRotation(), output);
  serializeVector3(c.getTranslation(), output);
}

/** Deserialize a coordinate frame from a binary input stream. */
inline THEA_API void
deserializeCoordinateFrame3(CoordinateFrame3 & c, BinaryInputStream & input)
{
  Matrix3 rot;
  Vector3 trans;
  deserializeMatrix3(rot, input);
  deserializeVector3(trans, input);
  c = CoordinateFrame3(RigidTransform3::_fromAffine(AffineTransform3(rot, trans)));
}

/** Serialize a 3D plane to a binary output stream. */
inline THEA_API void
serializePlane3(Plane3 const & plane, BinaryOutputStream & output)
{
  Real a, b, c, d;
  plane.getEquation(a, b, c, d);
  output.writeFloat32(a);
  output.writeFloat32(b);
  output.writeFloat32(c);
  output.writeFloat32(d);
}

/** Deserialize a 3D plane from a binary input stream. */
inline THEA_API void
deserializePlane3(Plane3 & plane, BinaryInputStream & input)
{
  float a = input.readFloat32();
  float b = input.readFloat32();
  float c = input.readFloat32();
  float d = input.readFloat32();
  plane = Plane3::fromEquation(a, b, c, d);
}

/**
 * Serialize a string to a binary output stream, aligning the output to a 4-byte boundary. The format is:
 * - Length of string (32-bit integer)
 * - Characters of string ('length' bytes, no null termination)
 * - Arbitrary padding bytes to align to 4-bytes boundary
 *
 * @see deserializeAlignedString
 */
inline THEA_API void
serializeAlignedString(std::string const & s, BinaryOutputStream & output)
{
  static const char zero[4] = { 0, 0, 0, 0 };

  int length = (int)s.length();
  output.writeInt32(length);
  output.writeBytes(s.c_str(), length);
  output.writeBytes(zero, (4 - (length & 0x03)) & 0x03);  // padding, zero not necessary but let's use it anyway
}

/**
 * Read an aligned string from a binary input stream. The format is:
 * - Length of string (32-bit integer)
 * - Characters of string ('length' bytes, no null termination)
 * - Arbitrary padding bytes to align to 4-bytes boundary
 *
 * @see serializeAlignedString
 */
inline THEA_API std::string
deserializeAlignedString(BinaryInputStream & input)
{
  int length = (int)input.readUInt32();  // string length
  std::string s = input.readString(length);  // string chars
  input.skip((4 - (length & 0x03)) & 0x03);  // padding
  return s;
}

// [Internal] Reads a line of input from the stream, bounded by any of the standard newline characters.
inline THEA_API std::string
readLine(BinaryInputStream & input)
{
  std::ostringstream oss;
  char c = 0;
  while (input.hasMore())
  {
    c = (char)input.readInt8();
    if (c == '\r' || c == '\n' || c == '\f')
      break;

    oss << c;
  }

  if (input.hasMore() && c == '\r')  // extract the CR-LF pair if it is present
  {
    c = (char)input.readInt8();
    if (c != '\n')
      input.skip(-1);
  }

  return oss.str();
}

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Serializable)
THEA_DECL_EXTERN_SMART_POINTERS(Thea::SerializableFactory)

#endif
