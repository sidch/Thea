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

#include "GLVAR.hpp"
#include "GLCaps.hpp"
#include <cstring>

namespace Thea {
namespace Graphics {
namespace GL {

GLVAR::GLVAR()
: area(NULL), capacity(0), pointer(NULL), generation(-1), gl_type(-1), num_components(0), elem_size(0), gl_target(-1),
  num_elems(0)
{
}

GLVAR::GLVAR(GLVARArea * area_, long num_bytes)
: area(area_), capacity(num_bytes), pointer(NULL), generation(-1), gl_type(-1), num_components(0), elem_size(0), gl_target(-1),
  num_elems(0)
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

#define GLVAR_UPDATE_ARRAY(func_name, type_, gl_type_, num_components_, gl_target_) \
  void \
  GLVAR::func_name(long start_elem, long num_elems_to_update, type_ const * array) \
  { \
    if (num_elems_to_update <= 0) return; \
    \
    alwaysAssertM(isValid(), "GLVAR: Can't update invalid VAR"); \
    \
    if (num_elems > 0) \
    { \
      alwaysAssertM(gl_type == gl_type_ && num_components == num_components_ && gl_target == gl_target_, \
                    "GLVAR: Can't update non-empty VAR with elements of a different type"); \
    } \
    else \
    { \
      gl_type = gl_type_; \
      num_components = num_components_; \
      elem_size = sizeof(type_); \
      gl_target = gl_target_; \
    } \
    \
    long offset_bytes = start_elem * elem_size; \
    long num_bytes = num_elems_to_update * elem_size; \
    long total_size = offset_bytes + num_bytes; \
    if (total_size > capacity) \
      throw Error("GLVAR: Can't update beyond end of VAR"); \
    \
    if (start_elem + num_elems_to_update > num_elems) \
      num_elems = start_elem + num_elems_to_update; \
    \
    uploadToGraphicsSystem(offset_bytes, num_bytes, array); \
  }

#define GLVAR_UPDATE_VECTOR_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateVectors, type_, gl_type_, num_components_, GL_ARRAY_BUFFER_ARB)

#define GLVAR_UPDATE_COLOR_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateColors, type_, gl_type_, num_components_, GL_ARRAY_BUFFER_ARB)

#define GLVAR_UPDATE_INDEX_ARRAY(type_, gl_type_, num_components_) \
  GLVAR_UPDATE_ARRAY(updateIndices, type_, gl_type_, num_components_, GL_ELEMENT_ARRAY_BUFFER_ARB)

GLVAR_UPDATE_VECTOR_ARRAY(float,    GL_FLOAT,  1)
GLVAR_UPDATE_VECTOR_ARRAY(Vector2,  GL_FLOAT,  2)
GLVAR_UPDATE_VECTOR_ARRAY(Vector3,  GL_FLOAT,  3)
GLVAR_UPDATE_VECTOR_ARRAY(Vector4,  GL_FLOAT,  4)

GLVAR_UPDATE_COLOR_ARRAY(ColorL,      GL_FLOAT,           1)
GLVAR_UPDATE_COLOR_ARRAY(ColorL8,     GL_UNSIGNED_BYTE,   1)
GLVAR_UPDATE_COLOR_ARRAY(ColorL16,    GL_UNSIGNED_SHORT,  1)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGB,    GL_FLOAT,           3)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGB8,   GL_UNSIGNED_BYTE,   3)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGBA,   GL_FLOAT,           4)
GLVAR_UPDATE_COLOR_ARRAY(ColorRGBA8,  GL_UNSIGNED_BYTE,   4)

GLVAR_UPDATE_INDEX_ARRAY(uint8,   GL_UNSIGNED_BYTE,   1)
GLVAR_UPDATE_INDEX_ARRAY(uint16,  GL_UNSIGNED_SHORT,  1)
GLVAR_UPDATE_INDEX_ARRAY(uint32,  GL_UNSIGNED_INT,    1)

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
GLVAR::uploadToGraphicsSystem(long offset_bytes, long num_bytes, void const * data)
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
