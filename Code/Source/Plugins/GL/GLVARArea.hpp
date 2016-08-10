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
    GLVARArea(GLRenderSystem * render_system, char const * name_, long capacity_, Usage usage, bool gpu_memory_ = true);

    /** Destructor. Automatically destroys all associated vertex arrays. */
    ~GLVARArea();

    /** Get the parent rendersystem. */
    GLRenderSystem * getRenderSystem() const { return render_system; }

    /** Get a string describing the storage area. */
    std::string toString() const;

    char const * getName() const { return name.c_str(); }
    void reset();
    long getCapacity() const { return capacity; }
    long getAllocatedSize() const { return allocated_size; }
    bool inGPUMemory() const { return gpu_memory; }

    VAR * createArray(long num_bytes);
    void destroyArray(VAR * array);

    /** Get the OpenGL index for this buffer.*/
    GLuint getGLBuffer() const { return gl_buffer; }

    /** Get a pointer to the first byte of the storage area.*/
    void * getBasePointer() const { return base_pointer; }

    /** Get the current generation (i.e. the number of times reset() has been called). */
    int getCurrentGeneration() const { return generation; }

    /** Mark an extra block as allocated. */
    void incrementAllocated(long inc) { allocated_size += inc; }

  private:
    /** Destroy all vertex arrays that are allocated from this area. */
    void destroyAllArrays();

    typedef TheaUnorderedSet<VAR *> VARSet;

    GLRenderSystem * render_system;
    std::string name;
    long capacity;
    bool gpu_memory;
    GLuint gl_buffer;
    uint8 * base_pointer;
    int generation;
    long allocated_size;
    VARSet vars;

}; // class GLVARArea

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
