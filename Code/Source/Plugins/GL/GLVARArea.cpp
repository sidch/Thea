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

#include "GLVARArea.hpp"
#include "GLCaps.hpp"
#include "GLVAR.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

GLVARArea::GLVARArea(GLRenderSystem * render_system_, char const * name_, long capacity_, Usage usage, bool gpu_memory_)
: render_system(render_system_), name(name_), capacity(capacity_), gpu_memory(gpu_memory_), generation(0), allocated_size(0)
{
  if (gpu_memory && !THEA_GL_SUPPORTS(ARB_vertex_buffer_object))
    throw Error(std::string(getName()) + ": OpenGL vertex/index buffers in GPU memory are not supported");

  assert(capacity > 0);

  if (gpu_memory)
  {
    { GLClientScope scope(GL_CLIENT_VERTEX_ARRAY_BIT);

      glGenBuffersARB(1, &gl_buffer);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_buffer);
      THEA_CHECK_GL_OK

      GLenum gl_usage;
      switch (usage)
      {
        case Usage::WRITE_EVERY_FRAME:  gl_usage = GL_STREAM_DRAW_ARB;  break;
        case Usage::WRITE_OCCASIONALLY: gl_usage = GL_DYNAMIC_DRAW_ARB; break;
        case Usage::WRITE_ONCE:         gl_usage = GL_STATIC_DRAW_ARB;  break;
        default:                        gl_usage = GL_STREAM_DRAW_ARB;
      }

      // Load some (undefined) data to initialize the buffer
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, (GLsizei)capacity, NULL, gl_usage);
      THEA_CHECK_GL_OK

      base_pointer = NULL;

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
      THEA_CHECK_GL_OK
    }
  }
  else
  {
    gl_buffer = 0;
    base_pointer = new uint8[capacity];
  }
}

GLVARArea::~GLVARArea()
{
  if (capacity <= 0)  // how did this happen?
    return;

  destroyAllArrays();

  if (gpu_memory)
    glDeleteBuffersARB(1, &gl_buffer);
  else
    delete [] base_pointer;
}

void
GLVARArea::destroyAllArrays()
{
  for (VARSet::iterator vi = vars.begin(); vi != vars.end(); ++vi)
    delete *vi;

  vars.clear();
}

std::string
GLVARArea::toString() const
{
  std::ostringstream oss;
  oss << "OpenGL VARArea of capacity " << getCapacity() << " bytes in " << (gpu_memory ? "GPU" : "CPU") << " memory, of which "
      << getAllocatedSize() << " bytes are currently allocated";
  return oss.str();
}

void
GLVARArea::reset()
{
  destroyAllArrays();

  generation++;
  allocated_size = 0;
}

VAR *
GLVARArea::createArray(long num_bytes)
{
  VAR * v = new GLVAR(this, num_bytes);
  if (v) vars.insert(v);
  return v;
}

void
GLVARArea::destroyArray(VAR * array)
{
  vars.erase(array);
  delete array;
}

} // namespace GL
} // namespace Graphics
} // namespace Thea
