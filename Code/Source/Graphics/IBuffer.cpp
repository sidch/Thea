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

#include "../MatVec.hpp"
#include "../Colors.hpp"

namespace Thea {

// We want to be able to interpret a buffer of (say) Vector3's as a tightly packed list of scalars. This isn't really necessary
// for the interface of IBuffer, but it makes things simpler and more efficient for the end user, so we'll make this a global
// requirement.
static_assert(sizeof(Vector2)    == 2 * sizeof(Real),   "IBuffer: Vector2 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(Vector3)    == 3 * sizeof(Real),   "IBuffer: Vector3 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(Vector4)    == 4 * sizeof(Real),   "IBuffer: Vector4 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorL)     == 1 * sizeof(Real),   "IBuffer: ColorL has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorL8)    == 1 * sizeof(uint8),  "IBuffer: ColorL8 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorL16)   == 1 * sizeof(uint16), "IBuffer: ColorL16 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorRgb)   == 3 * sizeof(Real),   "IBuffer: ColorRgb has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorRgb8)  == 3 * sizeof(uint8),  "IBuffer: ColorRgb8 has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorRgba)  == 4 * sizeof(Real),   "IBuffer: ColorRgba has padding, can't be tightly packed in a buffer");
static_assert(sizeof(ColorRgba8) == 4 * sizeof(uint8),  "IBuffer: ColorRgba8 has padding, can't be tightly packed in a buffer");

} // namespace Thea
