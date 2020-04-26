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

#ifndef __Thea_Graphics_GL_GLVARArea_hpp__
#define __Thea_Graphics_GL_GLVARArea_hpp__

#include "../../Graphics/VARArea.hpp"
#include "../../UnorderedSet.hpp"
#include "GLCommon.hpp"
#include "GLHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Forward declarations
class GLRenderSystem;

/** An OpenGL Vertex Area Range storage area in either main or GPU memory. Not threadsafe even when using main memory. */
class THEA_GL_DLL_LOCAL GLVARArea : public VARArea
{
  public:
    /** Constructor. */
    GLVARArea(GLRenderSystem * render_system, char const * name_, int64 capacity_, Usage usage, int8 gpu_memory_ = true);

    /** Destructor. Automatically destroys all associated vertex arrays. */
    ~GLVARArea();

    /** Get the parent rendersystem. */
    GLRenderSystem * getRenderSystem() const { return render_system; }

    /** Get a string describing the storage area. */
    std::string toString() const;

    char const * getName() const { return name.c_str(); }
    void reset();
    int64 getCapacity() const { return capacity; }
    int64 getAllocatedSize() const { return allocated_size; }
    int8 inGPUMemory() const { return gpu_memory; }

    VAR * createArray(int64 num_bytes);
    void destroyArray(VAR * array);

    /** Get the OpenGL index for this buffer.*/
    GLuint getGLBuffer() const { return gl_buffer; }

    /** Get a pointer to the first byte of the storage area.*/
    void * getBasePointer() const { return base_pointer; }

    /** Get the current generation (i.e. the number of times reset() has been called). */
    int32 getCurrentGeneration() const { return generation; }

    /** Mark an extra block as allocated. */
    void incrementAllocated(int64 inc) { allocated_size += inc; }

  private:
    /** Destroy all vertex arrays that are allocated from this area. */
    void destroyAllArrays();

    typedef UnorderedSet<VAR *> VARSet;

    GLRenderSystem * render_system;
    std::string name;
    int64 capacity;
    int8 gpu_memory;
    GLuint gl_buffer;
    uint8 * base_pointer;
    int32 generation;
    int64 allocated_size;
    VARSet vars;

}; // class GLVARArea

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
