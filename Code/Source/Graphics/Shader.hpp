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

#ifndef __Thea_Graphics_Shader_hpp__
#define __Thea_Graphics_Shader_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Colors.hpp"
#include "../MatrixMN.hpp"
#include "../VectorN.hpp"
#include "../NamedObject.hpp"
#include "Texture.hpp"

namespace Thea {
namespace Graphics {

/** An interface for a shader. */
class THEA_API Shader : public AbstractNamedObject
{
  public:
    /** %Shader module types (enum class). */
    struct THEA_API ModuleType
    {
      /** Supported values. */
      enum Value
      {
        VERTEX    =  0,  ///< Vertex shader.
        FRAGMENT  =  1,  ///< Fragment (pixel) shader.
        GEOMETRY  =  2   ///< Geometry shader.
      };

      THEA_ENUM_CLASS_BODY(ModuleType)
    };

    /** Destructor. */
    virtual ~Shader() {}

    /**
     * Check if the shader is ready to be used for rendering or not. Typically this requires both a vertex and a fragment
     * program to be attached.
     */
    virtual bool isComplete() const = 0;

    /** Attach a program module to the shader from a file containing its source code. */
    virtual void attachModuleFromFile(ModuleType type, char const * path) = 0;

    /** Attach a program module to the shader from a string containing its source code. */
    virtual void attachModuleFromString(ModuleType type, char const * source) = 0;

    /** Check if the shader has an active uniform of the given name. */
    virtual bool hasUniform(char const * uniform_name) const = 0;

    /** Set a floating-point uniform. */
    virtual void setUniform(char const * uniform_name, float value) = 0;

    /** Set an integer uniform. */
    virtual void setUniform(char const * uniform_name, int value) = 0;

    /** Set a 2-vector uniform. */
    virtual void setUniform(char const * uniform_name, Vector2 const & value) = 0;

    /** Set a 3-vector uniform. */
    virtual void setUniform(char const * uniform_name, Vector3 const & value) = 0;

    /** Set a 4-vector uniform. */
    virtual void setUniform(char const * uniform_name, Vector4 const & value) = 0;

    /** Set a single-channel byte color uniform. */
    virtual void setUniform(char const * uniform_name, ColorL8 const & value) = 0;

    /** Set a single-channel floating-point color uniform. */
    virtual void setUniform(char const * uniform_name, ColorL const & value) = 0;

    /** Set a 3-channel byte color uniform. */
    virtual void setUniform(char const * uniform_name, ColorRGB8 const & value) = 0;

    /** Set a 3-channel floating-point color uniform. */
    virtual void setUniform(char const * uniform_name, ColorRGB const & value) = 0;

    /** Set a 4-channel byte color uniform. */
    virtual void setUniform(char const * uniform_name, ColorRGBA8 const & value) = 0;

    /** Set a 4-channel floating-point color uniform. */
    virtual void setUniform(char const * uniform_name, ColorRGBA const & value) = 0;

    /** Set a 2x2 matrix uniform. */
    virtual void setUniform(char const * uniform_name, Matrix2 const & value) = 0;

    /** Set a 3x3 matrix uniform. */
    virtual void setUniform(char const * uniform_name, Matrix3 const & value) = 0;

    /** Set a 4x4 matrix uniform. */
    virtual void setUniform(char const * uniform_name, Matrix4 const & value) = 0;

    /** Set a texture uniform. */
    virtual void setUniform(char const * uniform_name, Texture * value) = 0;

    /** Set a floating-point array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, float const * values) = 0;

    /** Set an integer array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, int const * values) = 0;

    /** Set a 2-vector array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Vector2 const * values) = 0;

    /** Set a 3-vector array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Vector3 const * values) = 0;

    /** Set a 4-vector array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Vector4 const * values) = 0;

    /** Set a single-channel byte color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorL8 const * values) = 0;

    /** Set a single-channel floating-point color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorL const * values) = 0;

    /** Set a 3-channel byte color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorRGB8 const * values) = 0;

    /** Set a 3-channel floating-point color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorRGB const * values) = 0;

    /** Set a 4-channel byte color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorRGBA8 const * values) = 0;

    /** Set a 4-channel floating-point color array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, ColorRGBA const * values) = 0;

    /** Set a 2x2 matrix array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Matrix2 const * values) = 0;

    /** Set a 3x3 matrix array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Matrix3 const * values) = 0;

    /** Set a 4x4 matrix array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Matrix4 const * values) = 0;

    /** Set a texture array uniform. */
    virtual void setUniform(char const * uniform_name, long num_values, Texture * const * values) = 0;

}; // class Shader

} // namespace Graphics
} // namespace Thea

#endif
