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

#ifndef __Thea_EnumClass_hpp__
#define __Thea_EnumClass_hpp__

#include "Platform.hpp"
#include "BasicStringAlg.hpp"
#include <string>
#include <utility>

namespace Thea {

// Define the body of an enum class.
#define THEA_ENUM_CLASS_BODY(name)                                                                                            \
    public:                                                                                                                   \
      name() {}                                                                                                               \
      template <typename EnumClassSrcT> explicit name(EnumClassSrcT value_) : value(static_cast<Value>(value_)) {}            \
      name(Value value_) : value(value_) {}                                                                                   \
      operator Value() const { return value; }                                                                                \
      bool operator==(Value other) const { return value == other; }                                                           \
      bool operator!=(Value other) const { return value != other; }                                                           \
      bool operator==(name const & other) const { return value == other.value; }                                              \
      bool operator!=(name const & other) const { return value != other.value; }                                              \
                                                                                                                              \
    private:                                                                                                                  \
      Value value;                                                                                                            \
                                                                                                                              \
    public:

// Read and write enums from binary I/O streams (requires explicitly including Serialization.hpp)
#define THEA_ENUM_CLASS_SERIALIZATION(name)                                                                                   \
    public:                                                                                                                   \
      void read(BinaryInputStream & input, Codec const & codec = CodecAuto())                                                \
      {                                                                                                                       \
        BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);                                                  \
        int v = (int)input.readInt32();                                                                                       \
        value = static_cast<Value>(v);                                                                                        \
      }                                                                                                                       \
                                                                                                                              \
      void write(BinaryOutputStream & output, Codec const & codec = CodecAuto()) const                                       \
      {                                                                                                                       \
        BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);                                                \
        output.writeInt32(static_cast<int>(value));                                                                           \
      }

// Put quotes around the result of a macro expansion.
#define THEA_ENUM_CLASS_STRINGIFY_(x) #x
#define THEA_ENUM_CLASS_STRINGIFY(x) THEA_ENUM_CLASS_STRINGIFY_(x)

// Begin sequence of mappings of enum values to strings. There can be more than one string registered for the same value. When
// serializing, the first registered string is used.
#define THEA_ENUM_CLASS_STRINGS_BEGIN(name)                                                                                   \
    private:                                                                                                                  \
      static std::pair<Value, std::string> const * THEA_ENUM_CLASS_stringMapping(std::size_t i)                               \
      {                                                                                                                       \
        static std::pair<Value, std::string> const STRS[] = {

// Register an enum value with the string representation of its name.
#define THEA_ENUM_CLASS_STRING_SELF(value) std::pair<Value, std::string>(value, THEA_ENUM_CLASS_STRINGIFY(value)),

// Register an enum value to an arbitrary string.
#define THEA_ENUM_CLASS_STRING(value, str) std::pair<Value, std::string>(value, str),

// End sequence of mappings of enum values to strings, and define toString(), fromString() and construct-from-string functions.
#define THEA_ENUM_CLASS_STRINGS_END(name)                                                                                     \
        };                                                                                                                    \
                                                                                                                              \
        return (i >= (sizeof(STRS) / sizeof(std::pair<Value, std::string>))) ? nullptr : &STRS[i];                            \
      }                                                                                                                       \
                                                                                                                              \
    public:                                                                                                                   \
      std::string const & toString() const                                                                                    \
      {                                                                                                                       \
        std::pair<Value, std::string> const * mapping = THEA_ENUM_CLASS_stringMapping(0);                                     \
        std::size_t i = 0;                                                                                                    \
        while (mapping)                                                                                                       \
        {                                                                                                                     \
          if (mapping->first == value) return mapping->second;                                                                \
          mapping = THEA_ENUM_CLASS_stringMapping(++i);                                                                       \
        }                                                                                                                     \
                                                                                                                              \
        throw THEA_ENUM_CLASS_STRINGIFY(name) ": Enum value has no string representation";                                    \
      }                                                                                                                       \
                                                                                                                              \
      explicit name(std::string const & str, bool ignore_case = false)                                                        \
      {                                                                                                                       \
        if (!fromString(str, ignore_case))                                                                                    \
          throw THEA_ENUM_CLASS_STRINGIFY(name) ": Could not convert string to enum value";                                   \
      }                                                                                                                       \
                                                                                                                              \
      bool fromString(std::string str, bool ignore_case = false)                                                              \
      {                                                                                                                       \
        if (ignore_case)                                                                                                      \
          str = toLower(str);                                                                                                 \
                                                                                                                              \
        std::pair<Value, std::string> const * mapping = THEA_ENUM_CLASS_stringMapping(0);                                     \
        std::size_t i = 0;                                                                                                    \
        while (mapping)                                                                                                       \
        {                                                                                                                     \
          if ((!ignore_case && mapping->second == str) || (ignore_case && toLower(mapping->second) == str))                   \
          {                                                                                                                   \
            value = mapping->first;                                                                                           \
            return true;                                                                                                      \
          }                                                                                                                   \
                                                                                                                              \
          mapping = THEA_ENUM_CLASS_stringMapping(++i);                                                                       \
        }                                                                                                                     \
                                                                                                                              \
        return false;                                                                                                         \
      }

} // namespace Thea

#endif
