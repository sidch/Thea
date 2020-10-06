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

#include "GlBuffer.hpp"
#include "GlCaps.hpp"
#include <cstring>

namespace Thea {
namespace Graphics {
namespace Gl {

GlBuffer::GlBuffer()
: pool(nullptr), capacity(0), pointer(nullptr), generation(-1), gl_type(-1), value_size(0), gl_target(-1), num_values(0)
{
}

GlBuffer::GlBuffer(GlBufferPool * pool_, int64 num_bytes)
: pool(pool_), capacity(num_bytes), pointer(nullptr), generation(-1), gl_type(-1), value_size(0), gl_target(-1), num_values(0)
{
  alwaysAssertM(pool_, "GlBuffer: Valid buffer pool required");
  alwaysAssertM(num_bytes > 0, "GlBuffer: Capacity must be greater than zero");

  pool = pool_;

  if (num_bytes > pool_->getAvailableSize())
    throw Error("GlBuffer: Not enough free space in buffer pool to create buffer");

  generation = pool_->getCurrentGeneration();
  pointer = (uint8 *)pool_->getBasePointer() + pool_->getAllocatedSize();

  pool_->incrementAllocated(num_bytes);
}

std::string
GlBuffer::toString() const
{
  std::ostringstream oss;
  oss << "OpenGL buffer (" << capacity << " bytes)";
  return oss.str();
}

namespace GlInternal {

static_assert(sizeof(Real) == sizeof(float32) || sizeof(Real) == sizeof(float64),
              "GlBuffer: Real number type must be either float32 or float64");
GLenum GL_REAL_TYPE = (sizeof(Real) == sizeof(float32) ? GL_FLOAT : GL_DOUBLE);

GLenum
getGlType(NumericType type)
{
  switch (type)
  {
    case NumericType::INT8    : return GL_BYTE;
    case NumericType::INT16   : return GL_SHORT;
    case NumericType::INT32   : return GL_INT;
    case NumericType::UINT8   : return GL_UNSIGNED_BYTE;
    case NumericType::UINT16  : return GL_UNSIGNED_SHORT;
    case NumericType::UINT32  : return GL_UNSIGNED_INT;
    case NumericType::FLOAT32 : return GL_FLOAT;  // one of these will catch NumericType::REAL
    case NumericType::FLOAT64 : return GL_DOUBLE;
    default: return GL_INVALID_ENUM;
  }
}

} // namespace GlInternal

int8
GlBuffer::clear()
{
  gl_type = -1;
  value_size = -1;
  gl_target = -1;
  num_values = 0;

  return true;
}

int8
GlBuffer::updateAttributes(int64 start_value, int64 num_values_to_update, int32 ncomp, int32 type, void const * src)
{
  return update(start_value, num_values_to_update, ncomp, type, GL_ARRAY_BUFFER_ARB, src);
}

int8
GlBuffer::updateIndices(int64 start_value, int64 num_values_to_update, int32 type, void const * src)
{
  return update(start_value, num_values_to_update, 1, type, GL_ELEMENT_ARRAY_BUFFER_ARB, src);
}

int8
GlBuffer::update(int64 start_value, int64 num_values_to_update, int ncomp, int32 type, GLenum gl_target_, void const * src)
{
  if (num_values_to_update <= 0) return true;
  if (!isValid()) { THEA_ERROR << "GlBuffer: Can't update invalid buffer"; return GlCaps::setError(); }
  if (ncomp < 1)
  { THEA_ERROR << "GlBuffer: Can't update buffer with values of dimensionality <= 0"; return GlCaps::setError(); }

  GLenum gl_type_ = GlInternal::getGlType(NumericType(type));
  if (gl_type_ == GL_INVALID_ENUM) { THEA_ERROR << "GlBuffer: Invalid data type"; return GlCaps::setError(); }

  if (num_values > 0)
  {
    if (gl_type != gl_type_ || gl_target != gl_target_)
    {
      THEA_ERROR << "GlBuffer: Can't update non-empty buffer with values of a different type";
      return GlCaps::setError();
    }
  }
  else
  {
    num_components = ncomp;
    component_type = type;
    gl_type = gl_type_;
    value_size = ncomp * NumericType(type).numBits() / 8;
    gl_target = gl_target_;
  }

  int64 offset_bytes = start_value * value_size;
  int64 num_bytes = num_values_to_update * value_size;
  int64 total_size = offset_bytes + num_bytes;
  if (total_size > capacity)
  {
    THEA_ERROR << "GlBuffer: Can't update beyond end of buffer";
    return GlCaps::setError();
  }

  if (start_value + num_values_to_update > num_values)
    num_values = start_value + num_values_to_update;

  return uploadToGraphicsSystem(offset_bytes, num_bytes, src);
}

int8
GlBuffer::uploadToGraphicsSystem(int64 offset_bytes, int64 num_bytes, void const * data)
{
  if (!isValid()) { THEA_ERROR << "GlBuffer: Cannot upload invalid buffer to graphics system"; return GlCaps::setError(); }

  void * ptr = (void *)((uint8 *)pointer + offset_bytes);

  if (pool->inGpuMemory())
  {
    { GlClientScope scope(GL_CLIENT_VERTEX_ARRAY_BIT);

      glBindBufferARB(gl_target, pool->getGlBuffer());
      glBufferSubDataARB(gl_target, (GLintptrARB)ptr, (GLsizeiptr)num_bytes, data);
      glBindBufferARB(gl_target, 0);
    }
  }
  else
    std::memcpy(ptr, data, num_bytes);

  return true;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
