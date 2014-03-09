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

#ifndef __Thea_Graphics_GLShader_hpp__
#define __Thea_Graphics_GLShader_hpp__

#include "../../Graphics/Shader.hpp"
#include "../../Array.hpp"
#include "../../Map.hpp"
#include "GLCommon.hpp"
#include "GLTexture.hpp"
#include "GLHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Forward declarations
class GLRenderSystem;

/** An OpenGL shader. */
class THEA_GL_DLL_LOCAL GLShader : public Shader
{
  public:
    /** Constructor. */
    GLShader(GLRenderSystem * render_system_, char const * name_);

    /** Destructor. */
    ~GLShader();

    /** Get the parent rendersystem. */
    GLRenderSystem * getRenderSystem() const { return render_system; }

    char const * getName() const { return name.c_str(); }

    bool isComplete() const { return complete; }

    void attachModuleFromFile(ModuleType type, char const * path);
    void attachModuleFromString(ModuleType type, char const * source);

    bool hasUniform(char const * uniform_name) const { return uniforms.find(uniform_name) != uniforms.end(); }

    void setUniform(char const * uniform_name, float value);
    void setUniform(char const * uniform_name, int value);
    void setUniform(char const * uniform_name, Vector2 const & value);
    void setUniform(char const * uniform_name, Vector3 const & value);
    void setUniform(char const * uniform_name, Vector4 const & value);
    void setUniform(char const * uniform_name, ColorL8 const & value);
    void setUniform(char const * uniform_name, ColorL const & value);
    void setUniform(char const * uniform_name, ColorRGB8 const & value);
    void setUniform(char const * uniform_name, ColorRGB const & value);
    void setUniform(char const * uniform_name, ColorRGBA8 const & value);
    void setUniform(char const * uniform_name, ColorRGBA const & value);
    void setUniform(char const * uniform_name, Matrix2 const & value);
    void setUniform(char const * uniform_name, Matrix3 const & value);
    void setUniform(char const * uniform_name, Matrix4 const & value);
    void setUniform(char const * uniform_name, Texture * value);

    void setUniform(char const * uniform_name, long num_values, float const * values);
    void setUniform(char const * uniform_name, long num_values, int const * values);
    void setUniform(char const * uniform_name, long num_values, Vector2 const * values);
    void setUniform(char const * uniform_name, long num_values, Vector3 const * values);
    void setUniform(char const * uniform_name, long num_values, Vector4 const * values);
    void setUniform(char const * uniform_name, long num_values, ColorL8 const * values);
    void setUniform(char const * uniform_name, long num_values, ColorL const * values);
    void setUniform(char const * uniform_name, long num_values, ColorRGB8 const * values);
    void setUniform(char const * uniform_name, long num_values, ColorRGB const * values);
    void setUniform(char const * uniform_name, long num_values, ColorRGBA8 const * values);
    void setUniform(char const * uniform_name, long num_values, ColorRGBA const * values);
    void setUniform(char const * uniform_name, long num_values, Matrix2 const * values);
    void setUniform(char const * uniform_name, long num_values, Matrix3 const * values);
    void setUniform(char const * uniform_name, long num_values, Matrix4 const * values);
    void setUniform(char const * uniform_name, long num_values, Texture * const * values);

    /** Link the various modules of the shader into a single program. */
    void link();

    /** Use the shader for rendering. */
    void use();

    /** Get the OpenGL ID of the shader. */
    GLhandleARB getGLID() const { return program_id; }

  private:
    /** A value for a uniform variable. */
    struct UniformValue
    {
      float f_val;
      int i_val;
      TheaArray<float> f_array;
      TheaArray<int> i_array;
      GLTexture * texture;
    };

    /** Data related to an uniform variable. */
    struct UniformData
    {
      GLenum type;
      GLint size;
      GLint location;
      int texunit;
      bool has_value;
      UniformValue value;
      bool requires_rebind;

      /** Constructor. */
      UniformData() : has_value(false), requires_rebind(false) {}

      /** Note that the value has been changed. */
      void valueChanged() { has_value = true; requires_rebind = true; }
    };

    /** A set of uniforms read from source code. */
    typedef TheaMap<std::string, UniformData> Uniforms;

    /** Read the list of active uniforms in the shader object. */
    void readActiveUniforms();

    /** Bind the user-provided uniforms to the shader object. */
    void bindUniforms();

    /** Check if a build step (compile or link) succeeded, and throw a custom error if it did not. */
    void checkBuildStatus(GLhandleARB obj_id, GLenum status_field, std::string const & error_msg);

    GLRenderSystem * render_system;
    std::string name;
    bool complete;
    bool linked;
    bool has_vertex_module;
    bool has_fragment_module;
    GLhandleARB program_id;
    Uniforms uniforms;

}; // class GLShader

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
