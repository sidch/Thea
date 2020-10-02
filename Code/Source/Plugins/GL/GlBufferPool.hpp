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

#ifndef __Thea_Graphics_Gl_GlBufferPool_hpp__
#define __Thea_Graphics_Gl_GlBufferPool_hpp__

#include "../../Graphics/IBufferPool.hpp"
#include "../../UnorderedSet.hpp"
#include "GlCommon.hpp"
#include "GlHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

// Forward declarations
class GlRenderSystem;

/** An OpenGL buffer pool in either main or GPU memory. Not threadsafe even when using main memory. */
class THEA_GL_DLL_LOCAL GlBufferPool : public IBufferPool
{
  public:
    /** Constructor. */
    GlBufferPool(GlRenderSystem * render_system, char const * name_, int64 capacity_, int32 usage, int8 gpu_memory_ = true);

    /** Destructor. Automatically destroys all associated buffers. */
    ~GlBufferPool();

    /** Get the parent rendersystem. */
    GlRenderSystem * getRenderSystem() const { return render_system; }

    /** Get a string describing the storage pool. */
    std::string toString() const;

    char const * THEA_ICALL getName() const { return name.c_str(); }
    int8 THEA_ICALL reset();
    int64 THEA_ICALL getCapacity() const { return capacity; }
    int64 THEA_ICALL getAllocatedSize() const { return allocated_size; }
    int64 THEA_ICALL getAvailableSize() const { return capacity - allocated_size; }
    int8 THEA_ICALL inGpuMemory() const { return gpu_memory; }
    IBuffer * THEA_ICALL createBuffer(int64 num_bytes);
    int8 THEA_ICALL destroyBuffer(IBuffer * buf);

    /** Get the OpenGL index for this buffer.*/
    GLuint getGlBuffer() const { return gl_buffer; }

    /** Get a pointer to the first byte of the storage pool.*/
    void * getBasePointer() const { return base_pointer; }

    /** Get the current generation (i.e. the number of times reset() has been called). */
    int32 getCurrentGeneration() const { return generation; }

    /** Mark an extra block as allocated. */
    void incrementAllocated(int64 inc) { allocated_size += inc; }

  private:
    /** Destroy all buffers that are allocated from this pool. */
    void destroyAllBuffers();

    typedef UnorderedSet<IBuffer *> BufferSet;

    GlRenderSystem * render_system;
    std::string name;
    int64 capacity;
    int8 gpu_memory;
    GLuint gl_buffer;
    uint8 * base_pointer;
    int32 generation;
    int64 allocated_size;
    BufferSet buffers;

}; // class GlBufferPool

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
