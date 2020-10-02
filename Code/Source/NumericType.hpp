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

#ifndef __Thea_NumericType_hpp__
#define __Thea_NumericType_hpp__

#include "Platform.hpp"
#include "EnumClass.hpp"
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace Thea {

typedef  std::int8_t         int8;
typedef  std::int16_t        int16;
typedef  std::int32_t        int32;
typedef  std::int64_t        int64;

typedef  std::uint8_t        uint8;
typedef  std::uint16_t       uint16;
typedef  std::uint32_t       uint32;
typedef  std::uint64_t       uint64;

using    std::               size_t;

typedef  float               float32;  // assume IEEE 754
typedef  double              float64;  // assume IEEE 754
typedef  float32             Real;

/**
 * Type of a number (enum class). This values are guaranteed to have the following encoding:
 *   - Right-two hex digits: number of bits in type.
 *   - Third digit from right: signed (0) or unsigned (1).
 *   - Fourth digit from right: int (0) or float (1).
 */
struct NumericType
{
  enum Value
  {
    // See explanation of encoding in class description.
    INVALID     =   0x0000,
    INT8        =   0x0008,
    INT16       =   0x0010,
    INT32       =   0x0020,
    INT64       =   0x0040,
    UINT8       =   0x0108,
    UINT16      =   0x0110,
    UINT32      =   0x0120,
    UINT64      =   0x0140,
    FLOAT32     =   0x1020,
    FLOAT64     =   0x1040,
    REAL        =   (sizeof(Real) == 4 ? 0x1020 : 0x1040)  // duplicate values are ok in an enum
  };

  THEA_ENUM_CLASS_BODY(NumericType)

  THEA_ENUM_CLASS_STRINGS_BEGIN(NumericType)
    THEA_ENUM_CLASS_STRING(INVALID,   "invalid")
    THEA_ENUM_CLASS_STRING(INT8,      "int8")
    THEA_ENUM_CLASS_STRING(INT16,     "int16")
    THEA_ENUM_CLASS_STRING(INT32,     "int32")
    THEA_ENUM_CLASS_STRING(INT32,     "int64")
    THEA_ENUM_CLASS_STRING(UINT8,     "uint8")
    THEA_ENUM_CLASS_STRING(UINT16,    "uint16")
    THEA_ENUM_CLASS_STRING(UINT32,    "uint32")
    THEA_ENUM_CLASS_STRING(UINT64,    "uint64")
    THEA_ENUM_CLASS_STRING(FLOAT32,   "float32")
    THEA_ENUM_CLASS_STRING(FLOAT64,   "float64")
    // no string for REAL since it's a duplicate of either FLOAT32 or FLOAT64
  THEA_ENUM_CLASS_STRINGS_END(NumericType)

  /** Convert a numeric C++ type to the corresponding enum. Use as <code>NumericType::From<type>::value</code>. */
  template <typename T> struct From;

  /** Convert an enum to the corresponding numeric C++ type. Use as <code>NumericType::To<value>::type</code>. */
  template <Value> struct To;

  /** Get the number of bits in the type. */
  int numBits() const { return (int)(*this) & 0x00FF; }

  /** Is the type signed? */
  bool isSigned() const { return (int)(*this) & 0x0F00; }

  /** Is the type floating-point? */
  bool isFloat() const { return (int)(*this) & 0xF000; }

  /** Is the type integral? */
  bool isInteger() const { return !isFloat(); }

}; // struct NumericType

// Can't specialize these inside the outer class
template <> struct NumericType::From<int8   > { static Value const value = NumericType::INT8;    };
template <> struct NumericType::From<int16  > { static Value const value = NumericType::INT16;   };
template <> struct NumericType::From<int32  > { static Value const value = NumericType::INT32;   };
template <> struct NumericType::From<int64  > { static Value const value = NumericType::INT64;   };
template <> struct NumericType::From<uint8  > { static Value const value = NumericType::UINT8;   };
template <> struct NumericType::From<uint16 > { static Value const value = NumericType::UINT16;  };
template <> struct NumericType::From<uint32 > { static Value const value = NumericType::UINT32;  };
template <> struct NumericType::From<uint64 > { static Value const value = NumericType::UINT64;  };
template <> struct NumericType::From<float32> { static Value const value = NumericType::FLOAT32; };
template <> struct NumericType::From<float64> { static Value const value = NumericType::FLOAT64; };

template <> struct NumericType::To<NumericType::INT8   > { typedef int8    type; };
template <> struct NumericType::To<NumericType::INT16  > { typedef int16   type; };
template <> struct NumericType::To<NumericType::INT32  > { typedef int32   type; };
template <> struct NumericType::To<NumericType::INT64  > { typedef int64   type; };
template <> struct NumericType::To<NumericType::UINT8  > { typedef uint8   type; };
template <> struct NumericType::To<NumericType::UINT16 > { typedef uint16  type; };
template <> struct NumericType::To<NumericType::UINT32 > { typedef uint32  type; };
template <> struct NumericType::To<NumericType::UINT64 > { typedef uint64  type; };
template <> struct NumericType::To<NumericType::FLOAT32> { typedef float32 type; };
template <> struct NumericType::To<NumericType::FLOAT64> { typedef float64 type; };

/**
 * Check that <code>S</code> is a scalar type that can participate in arithmetic operations with field type <code>T</code>
 * without explicit casting.
 */
template <typename S, typename T>
struct IsCompatibleScalar
{
  static bool const value = std::is_arithmetic<S>::value || std::is_same<S, T>::value;

}; // struct IsCompatibleScalar

} // namespace Thea

#endif
