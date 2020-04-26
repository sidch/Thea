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
// First version: 2015
//
//============================================================================

#ifndef __Thea_CommonEnums_hpp__
#define __Thea_CommonEnums_hpp__

#include "Platform.hpp"
#include "EnumClass.hpp"
#include "NumericType.hpp"

namespace Thea {

/** Coordinate axis-aligned directions upto 4D (enum class). */
struct THEA_API AxisAlignedDirection
{
  /** Supported values. */
  enum Value
  {
    POS_X = 0,  ///< The positive X direction.
    NEG_X,      ///< The negative X direction.
    POS_Y,      ///< The positive Y direction.
    NEG_Y,      ///< The negative Y direction.
    POS_Z,      ///< The positive Z direction.
    NEG_Z,      ///< The negative Z direction.
    POS_W,      ///< The positive W direction.
    NEG_W       ///< The negative W direction.
  };

  THEA_ENUM_CLASS_BODY(AxisAlignedDirection)

  THEA_ENUM_CLASS_STRINGS_BEGIN(AxisAlignedDirection)
    THEA_ENUM_CLASS_STRING(POS_X,  "+X")
    THEA_ENUM_CLASS_STRING(NEG_X,  "-X")
    THEA_ENUM_CLASS_STRING(POS_Y,  "+Y")
    THEA_ENUM_CLASS_STRING(NEG_Y,  "-Y")
    THEA_ENUM_CLASS_STRING(POS_Z,  "+Z")
    THEA_ENUM_CLASS_STRING(NEG_Z,  "-Z")
    THEA_ENUM_CLASS_STRING(POS_W,  "+W")
    THEA_ENUM_CLASS_STRING(NEG_W,  "-W")
  THEA_ENUM_CLASS_STRINGS_END(AxisAlignedDirection)
};

/** Comparison operators (enum class). */
struct THEA_API CompareOp
{
  /** Supported values. */
  enum Value
  {
    EQUAL,      ///< Equality.
    NOT_EQUAL,  ///< Inequality.
    LESS,       ///< Less-than.
    LEQUAL,     ///< Less-than-or-equal.
    GREATER,    ///< Greater-than.
    GEQUAL      ///< Greater-than-or-equal.
  };

  THEA_ENUM_CLASS_BODY(CompareOp)

  THEA_ENUM_CLASS_STRINGS_BEGIN(CompareOp)
    THEA_ENUM_CLASS_STRING(EQUAL,      "==")
    THEA_ENUM_CLASS_STRING(NOT_EQUAL,  "!=")
    THEA_ENUM_CLASS_STRING(LESS,       "<")
    THEA_ENUM_CLASS_STRING(LEQUAL,     "<=")
    THEA_ENUM_CLASS_STRING(GREATER,    ">")
    THEA_ENUM_CLASS_STRING(GEQUAL,     ">=")
  THEA_ENUM_CLASS_STRINGS_END(CompareOp)
};

/** Common distance metrics. */
struct DistanceType
{
  /** Supported values. */
  enum Value
  {
    EUCLIDEAN,  ///< Euclidean distance.
    GEODESIC    ///< Geodesic distance on a surface.
  };

  THEA_ENUM_CLASS_BODY(DistanceType)

  THEA_ENUM_CLASS_STRINGS_BEGIN(DistanceType)
    THEA_ENUM_CLASS_STRING(EUCLIDEAN,  "euclidean")
    THEA_ENUM_CLASS_STRING(GEODESIC,   "geodesic")
  THEA_ENUM_CLASS_STRINGS_END(DistanceType)
};

/** %Endianness values (little-endian and big-endian) (enum class). Also has a function to check the machine endianness. */
struct THEA_API Endianness
{
  /** Supported values. */
  enum Value
  {
    BIG,
    LITTLE
  };

  THEA_ENUM_CLASS_BODY(Endianness)

  THEA_ENUM_CLASS_STRINGS_BEGIN(Endianness)
    THEA_ENUM_CLASS_STRING(BIG,    "big-endian")
    THEA_ENUM_CLASS_STRING(LITTLE, "little-endian")
  THEA_ENUM_CLASS_STRINGS_END(Endianness)

  /** Get the machine endian-ness. */
  static Endianness machine()
  {
    union
    {
      uint32 i;
      char c[4];
    } b = { 0x01020304 };

    return (b.c[0] == 1 ? BIG : LITTLE);
  }
};

} // namespace Thea

#endif
