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

#include "GLVAR.hpp"
#include "GLCaps.hpp"
#include <cstring>

namespace Thea {
namespace Graphics {
namespace GL {

GLVAR::GLVAR()
: area(nullptr), capacity(0), pointer(nullptr), generation(-1), gl_type(-1), num_components(0), elem_size(0), gl_target(-1),
  num_elems(0)
{
}

GLVAR::GLVAR(GLVARArea * area_, int64 num_bytes)
: area(area_), capacity(num_bytes), pointer(nullptr), generation(-1), gl_type(-1), num_components(0), elem_size(0),
  gl_target(-1), num_elems(0)
{
  alwaysAssertM(area_, "GLVAR: Valid VAR area required");
  alwaysAssertM(num_bytes > 0, "GLVAR: Capacity must be greater than zero");

  area = area_;

  if (num_bytes > area_->getAvailableSize())
    throw Error("GLVAR: Not enough free space in VAR area to create VAR");

  generation = area_->getCurrentGeneration();
  pointer = (uint8 *)area_->getBasePointer() + area_->getAllocatedSize();

  area_->incrementAllocated(num_bytes);
}

std::string
GLVAR::toString() const
{
  std::ostringstream oss;
  oss << "OpenGL VAR (" << capacity << " bytes)";
  return oss.str();
}

// We need to be able to interpret an array of (say) Vector3's as a tightly packed list of scalars
static_assert(sizeof(Vector2)    == 2 * sizeof(Real),   "GLVAR: Vector2 has padding, can't be tightly packed in an array");
static_assert(sizeof(Vector3)    == 3 * sizeof(Real),   "GLVAR: Vector3 has padding, can't be tightly packed in an array");
static_assert(sizeof(Vector4)    == 4 * sizeof(Real),   "GLVAR: Vector4 has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorL)     == 1 * sizeof(Real),   "GLVAR: ColorL has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorL8)    == 1 * sizeof(uint8),  "GLVAR: ColorL8 has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorL16)   == 1 * sizeof(uint16), "GLVAR: ColorL16 has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorRGB)   == 3 * sizeof(Real),   "GLVAR: ColorRGB has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorRGB8)  == 3 * sizeof(uint8),  "GLVAR: ColorRGB8 has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorRGBA)  == 4 * sizeof(Real),   "GLVAR: ColorRGBA has padding, can't be tightly packed in an array");
static_assert(sizeof(ColorRGBA8) == 4 * sizeof(uint8),  "GLVAR: ColorRGBA8 has padding, can't be tightly packed in an array");

static_assert(sizeof(Real) == sizeof(float32) || sizeof(Real) == sizeof(float64),
              "GLVAR: Real number type must be either 32-bit float64 or 64-bit float64");
GLenum GL_REAL_TYPE = (sizeof(Real) == sizeof(float32) ? GL_FLOAT : GL_DOUBLE);

#define GLVAR_UPDATE_ARRAY(func_name, type_, gl_type_, num_components_, gl_target_)                                           \
  void                                                                                                                        \
  GLVAR::func_name(int64 start_elem, int64 num_elems_to_update, type_ const * array)                                          \
  {                                                                                                                           \
    if (num_elems_to_update <= 0) return;                                                                                     \
                                                                                                                              \
    alwaysAssertM(isValid(), "GLVAR: Can't update invalid VAR");                                                              \
                                                                                                                              \
    if (num_elems > 0)                                                                                                        \
    {                                                                                                                         \
      alwaysAssertM(gl_type == gl_type_ && num_components == num_components_ && gl_target == gl_target_,                      \
                    "GLVAR: Can't update non-empty VAR with elements of a different type");                                   \
    }                                                                                                                         \
    else                                                                                                                      \
    {                                                                                                                         \
      gl_type = gl_type_;                                                                                                     \
      num_components = num_components_;                                                                                       \
      elem_size = sizeof(type_);                                                                                              \
      gl_target = gl_target_;                                                                                                 \
    }                                                                                                                         \
                                                                                                                              \
    int64 offset_bytes = start_elem * elem_size;                                                                              \
    int64 num_bytes = num_elems_to_update * elem_size;                                                                        \
    int64 total_size = offset_bytes + num_bytes;                                                                              \
    if (total_size > capacity)                                                                                                \
      throw Error("GLVAR: Can't update beyond end of VAR");                                                                   \
                                                                                                                              \
    if (start_elem + num_elems_to_update > num_elems)                                                                         \
      num_elems = start_elem + num_elems_to_update;                                                                           \
                                                                                                                              \
    uploadToGraphicsSystem(offset_bytes, num_bytes, array);                                                                   \
  }

#define GLVAR_UPDATE_VECTOR_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateVectors, type_, gl_type_, num_components_, GL_ARRAY_BUFFER_ARB)

#define GLVAR_UPDATE_COLOR_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateColors, type_, gl_type_, num_components_, GL_ARRAY_BUFFER_ARB)

#define GLVAR_UPDATE_INDEX_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateIndices, type_, gl_type_, num_components_, GL_ELEMENT_ARRAY_BUFFER_ARB)

GLVAR_UPDATE_VECTOR_ARRAY(float32,    GL_FLOAT,           1)
GLVAR_UPDATE_VECTOR_ARRAY(float64,    GL_DOUBLE,          1)

GLVAR_UPDATE_VECTOR_ARRAY(Vector2,    GL_REAL_TYPE,       2)
GLVAR_UPDATE_VECTOR_ARRAY(Vector3,    GL_REAL_TYPE,       3)
GLVAR_UPDATE_VECTOR_ARRAY(Vector4,    GL_REAL_TYPE,       4)

GLVAR_UPDATE_COLOR_ARRAY(ColorL,      GL_REAL_TYPE,       1)
GLVAR_UPDATE_COLOR_ARRAY(ColorL8,     GL_UNSIGNED_BYTE,   1)
GLVAR_UPDATE_COLOR_ARRAY(ColorL16,    GL_UNSIGNED_SHORT,  1)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGB,    GL_REAL_TYPE,       3)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGB8,   GL_UNSIGNED_BYTE,   3)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGBA,   GL_REAL_TYPE,       4)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGBA8,  GL_UNSIGNED_BYTE,   4)

GLVAR_UPDATE_INDEX_ARRAY(uint8,       GL_UNSIGNED_BYTE,   1)
GLVAR_UPDATE_INDEX_ARRAY(uint16,      GL_UNSIGNED_SHORT,  1)
GLVAR_UPDATE_INDEX_ARRAY(uint32,      GL_UNSIGNED_INT,    1)

void
GLVAR::clear()
{
  gl_type = -1;
  num_components = -1;
  elem_size = -1;
  gl_target = -1;
  num_elems = 0;
}

void
GLVAR::uploadToGraphicsSystem(int64 offset_bytes, int64 num_bytes, void const * data)
{
  assert(isValid());

  void * ptr = (void *)((uint8 *)pointer + offset_bytes);

  if (area->inGPUMemory())
  {
    { GLClientScope scope(GL_CLIENT_VERTEX_ARRAY_BIT);

      glBindBufferARB(gl_target, area->getGLBuffer());
      glBufferSubDataARB(gl_target, (GLintptrARB)ptr, (GLsizeiptr)num_bytes, data);
      glBindBufferARB(gl_target, 0);
    }
  }
  else
    std::memcpy(ptr, data, num_bytes);
}

} // namespace GL
} // namespace Graphics
} // namespace Thea
