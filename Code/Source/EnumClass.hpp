//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Stanford University
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
      template <typename T> explicit name(T value_) : value(static_cast<Value>(value_)) {}                                    \
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

// Serialize and deserialize enums from binary I/O streams (requires explicitly including Serialization.hpp)
#define THEA_ENUM_CLASS_SERIALIZATION(name)                                                                                   \
    public:                                                                                                                   \
      void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const                                   \
      {                                                                                                                       \
        output.setEndianness(Endianness::LITTLE);                                                                             \
        output.writeInt32(static_cast<int>(value));                                                                           \
      }                                                                                                                       \
                                                                                                                              \
      void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO())                                         \
      {                                                                                                                       \
        input.setEndianness(Endianness::LITTLE);                                                                              \
        int v = (int)input.readInt32();                                                                                       \
        value = static_cast<Value>(v);                                                                                        \
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
        return (i >= (sizeof(STRS) / sizeof(std::pair<Value, std::string>))) ? NULL : &STRS[i];                               \
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
