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

#ifndef __Thea_Graphics_GL_GLVAR_hpp__
#define __Thea_Graphics_GL_GLVAR_hpp__

#include "../../Graphics/VAR.hpp"
#include "GLCommon.hpp"
#include "GLHeaders.hpp"
#include "GLVARArea.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

/** An OpenGL Vertex Area Range object, which may be in main or GPU memory. */
class THEA_GL_DLL_LOCAL GLVAR : public VAR
{
  public:
    /** Default constructor. Creates an empty, invalid VAR. */
    GLVAR();

    /**
     * Constructor. Creates an empty VAR of the specified size. The VAR is not valid until it has been initialized with one of
     * the <code>update...()</code> functions. \a area_ must be non-null and \a num_bytes must be greater than zero.
     */
    GLVAR(GLVARArea * area_, int64 num_bytes);

    /** Get a string describing the VAR. */
    std::string toString() const;

    void updateVectors(int64 start_elem, int64 num_elems_to_update, float32 const * array);
    void updateVectors(int64 start_elem, int64 num_elems_to_update, float64 const * array);
    void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector2 const * array);
    void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector3 const * array);
    void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector4 const * array);

    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL8 const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL16 const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGB const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGB8 const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGBA const * array);
    void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGBA8 const * array);

    void updateIndices(int64 start_elem, int64 num_elems_to_update, uint8 const * array);
    void updateIndices(int64 start_elem, int64 num_elems_to_update, uint16 const * array);
    void updateIndices(int64 start_elem, int64 num_elems_to_update, uint32 const * array);

    void clear();

    int64 numElements() const { return num_elems; }
    int64 getCapacityInBytes() const { return capacity; }
    int8 isValid() const { return area && capacity > 0 && generation == area->getCurrentGeneration(); }

    /** The OpenGL data type of a single component (eg GL_FLOAT). */
    GLenum getGLType() const { return gl_type; }

    /** The number of components per element. */
    int32 getNumComponents() const { return num_components; }

    /** The size of an element in bytes. */
    int32 getElementSize() const { return elem_size; }

    /** The id of the OpenGL target. */
    int32 getGLTarget() const { return gl_target; }

    /** Get the VARArea where this VAR is stored. */
    GLVARArea * getArea() const { return area; }

    /** A pointer to the first element of the VAR. */
    void * getBasePointer() const { return pointer; }

    /** The generation of the parent VARArea when this VAR was created. */
    int32 getGeneration() const { return generation; }

  private:
    /** Upload source data to the graphics system. */
    void uploadToGraphicsSystem(int64 offset_bytes, int64 num_bytes, void const * data);

    GLVARArea * area;
    int64 capacity;
    void * pointer;
    int32 generation;

    GLenum gl_type;
    int32 num_components;
    int32 elem_size;
    GLenum gl_target;
    int64 num_elems;

}; // class GLVAR

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
