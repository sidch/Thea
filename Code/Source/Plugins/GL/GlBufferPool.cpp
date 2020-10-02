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

#include "GlBufferPool.hpp"
#include "GlCaps.hpp"
#include "GlBuffer.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

GlBufferPool::GlBufferPool(GlRenderSystem * render_system_, char const * name_, int64 capacity_, int32 usage, int8 gpu_memory_)
: render_system(render_system_), name(name_), capacity(capacity_), gpu_memory(gpu_memory_), generation(0), allocated_size(0)
{
  if (gpu_memory && !THEA_GL_SUPPORTS(ARB_vertex_buffer_object))
    throw Error(std::string(getName()) + ": OpenGL vertex/index buffers in GPU memory are not supported");

  assert(capacity > 0);

  if (gpu_memory)
  {
    { GlClientScope scope(GL_CLIENT_VERTEX_ARRAY_BIT);

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
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, (GLsizei)capacity, nullptr, gl_usage);
      THEA_CHECK_GL_OK

      base_pointer = nullptr;

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

GlBufferPool::~GlBufferPool()
{
  if (capacity <= 0)  // how did this happen?
    return;

  destroyAllBuffers();

  if (gpu_memory)
    glDeleteBuffersARB(1, &gl_buffer);
  else
    delete [] base_pointer;
}

void
GlBufferPool::destroyAllBuffers()
{
  for (BufferSet::iterator bi = buffers.begin(); bi != buffers.end(); ++bi)
    delete *bi;

  buffers.clear();
}

std::string
GlBufferPool::toString() const
{
  std::ostringstream oss;
  oss << "OpenGL BufferPool of capacity " << getCapacity() << " bytes in " << (gpu_memory ? "GPU" : "CPU")
      << " memory, of which " << getAllocatedSize() << " bytes are currently allocated";
  return oss.str();
}

int8
GlBufferPool::reset()
{
  destroyAllBuffers();

  generation++;
  allocated_size = 0;

  return true;
}

IBuffer *
GlBufferPool::createBuffer(int64 num_bytes)
{
  IBuffer * b = new GlBuffer(this, num_bytes);
  if (b) buffers.insert(b);
  return b;
}

int8
GlBufferPool::destroyBuffer(IBuffer * buf)
{
  buffers.erase(buf);
  delete buf;

  return true;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
