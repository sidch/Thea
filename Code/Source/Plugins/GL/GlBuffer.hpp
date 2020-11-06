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

#ifndef __Thea_Graphics_Gl_GlBuffer_hpp__
#define __Thea_Graphics_Gl_GlBuffer_hpp__

#include "../../Graphics/IBuffer.hpp"
#include "GlCommon.hpp"
#include "GlHeaders.hpp"
#include "GlBufferPool.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

/**
 * An OpenGL buffer, which may be in main or GPU memory. A GlBuffer object is safe to copy, all copies reference the same memory
 * area.
 */
class THEA_GL_DLL_LOCAL GlBuffer : public virtual IBuffer
{
  public:
    /** Default constructor. Creates an empty, invalid buffer. */
    GlBuffer();

    /**
     * Constructor. Creates an empty buffer of the specified size. The buffer is not valid until it has been initialized with
     * one of the <code>update...()</code> functions. \a pool_ must be non-null and \a num_bytes must be greater than zero.
     */
    GlBuffer(GlBufferPool * pool_, int64 num_bytes);

    /** Get a string describing the buffer. */
    std::string toString() const;

    int32 THEA_ICALL numComponents() const { return num_components; }
    int32 THEA_ICALL getComponentType() const { return component_type; }
    int64 THEA_ICALL numValues() const { return num_values; }
    int64 THEA_ICALL getCapacityInBytes() const { return capacity; }
    int8 THEA_ICALL isValid() const { return pool && capacity > 0 && generation == pool->getCurrentGeneration(); }
    int8 THEA_ICALL clear();
    int8 THEA_ICALL updateAttributes(int64 start_value, int64 num_values_to_update, int32 ncomp, int32 type, void const * src);
    int8 THEA_ICALL updateIndices(int64 start_value, int64 num_values_to_update, int32 type, void const * src);

    /** The OpenGL data type of a single component (eg GL_FLOAT). */
    GLenum getGlType() const { return gl_type; }

    /** The size of a value in bytes. */
    int32 getValueSize() const { return value_size; }

    /** The id of the OpenGL target. */
    int32 getGlTarget() const { return gl_target; }

    /** Get the BufferPool where this buffer is stored. */
    GlBufferPool * getPool() const { return pool; }

    /** A pointer to the first value of the buffer. */
    void * getBasePointer() const { return pointer; }

    /** The generation of the parent BufferPool when this buffer was created. */
    int32 getGeneration() const { return generation; }

  private:
    /** Helper function for updateAttributes() and updateIndices(). */
    int8 update(int64 start_value, int64 num_values_to_update, int ncomp, int32 type, GLenum gl_target_, void const * src);

    /** Upload source data to the graphics system. */
    int8 uploadToGraphicsSystem(int64 offset_bytes, int64 num_bytes, void const * data);

    GlBufferPool * pool;
    int64 capacity;
    void * pointer;
    int32 generation;

    int32 num_components;
    int32 component_type;
    GLenum gl_type;
    int32 value_size;
    GLenum gl_target;
    int64 num_values;

}; // class GlBuffer

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
