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

#include "GLShader.hpp"
#include "GLCaps.hpp"
#include "../../FileSystem.hpp"
#include <cstring>
#include <fstream>

namespace Thea {
namespace Graphics {
namespace GL {

GLShader::GLShader(GLRenderSystem * render_system_, char const * name_)
: render_system(render_system_), name(name_), complete(false), linked(true), has_vertex_module(false),
  has_fragment_module(false)
{
  if (!THEA_GL_SUPPORTS(ARB_shader_objects))
    throw Error("OpenGL programmable shaders are not supported");

  program_id = glCreateProgramObjectARB();
  THEA_CHECK_GL_OK
}

GLShader::~GLShader()
{
  glDeleteObjectARB(program_id);
}

void
GLShader::attachModuleFromFile(ModuleType type, char const * path)
{
  std::ifstream in(path);
  if (!in)
    throw Error(std::string(getName()) + ": Shader module '" + std::string(path) + "' not found");

  std::string source = FileSystem::readWholeFile(path);
  attachModuleFromString(type, source.c_str());
}

void
GLShader::attachModuleFromString(ModuleType type, char const * source)
{
  GLenum gltype;
  switch (type)
  {
    case ModuleType::VERTEX:   gltype = GL_VERTEX_SHADER_ARB; break;
    case ModuleType::FRAGMENT: gltype = GL_FRAGMENT_SHADER_ARB; break;
    case ModuleType::GEOMETRY:
    {
      if (THEA_GL_SUPPORTS(EXT_geometry_shader4))
        gltype = GL_GEOMETRY_SHADER_EXT;
      else if (THEA_GL_SUPPORTS(ARB_geometry_shader4))
        gltype = GL_GEOMETRY_SHADER_ARB;  // this should have exactly the same value as the EXT version, but let's play safe
      else
        throw Error(std::string(getName()) + ": This OpenGL installation does not support geometry shaders");

      break;
    }
    default:
      throw Error(std::string(getName()) + ": Unknown module type");
  }

  GLhandleARB module_id = glCreateShaderObjectARB(gltype);
  THEA_CHECK_GL_OK

  int src_length = (int)std::strlen(source);
  GLcharARB const * src_ptr = (GLcharARB const *)source;
  glShaderSourceARB(module_id, 1, &src_ptr, &src_length);
  THEA_CHECK_GL_OK

  glCompileShaderARB(module_id);
  checkBuildStatus(module_id, GL_OBJECT_COMPILE_STATUS_ARB, "Failed to compile shader module");

  glAttachObjectARB(program_id, module_id);
  THEA_CHECK_GL_OK

  if (type == ModuleType::VERTEX)   has_vertex_module   = true;
  if (type == ModuleType::FRAGMENT) has_fragment_module = true;
  if (has_vertex_module && has_fragment_module) complete = true;

  // Update the list of active shader uniforms
  if (complete) readActiveUniforms();

  linked = false;
}

static bool
GLShader__isSamplerType(GLenum type)
{
  return type == GL_SAMPLER_1D_ARB
      || type == GL_SAMPLER_2D_ARB
      || type == GL_SAMPLER_3D_ARB
      || type == GL_SAMPLER_CUBE_ARB
      || type == GL_SAMPLER_2D_RECT_ARB
      || type == GL_SAMPLER_2D_SHADOW_ARB
      || type == GL_SAMPLER_2D_RECT_SHADOW_ARB;
}

static GLenum
GLShader__toCanonicalType(GLenum type)
{
  switch (type)
  {
    case GL_SAMPLER_1D_ARB: return GL_TEXTURE_1D;
    case GL_SAMPLER_2D_ARB: return GL_TEXTURE_2D;
    case GL_SAMPLER_3D_ARB: return GL_TEXTURE_3D;
    case GL_SAMPLER_CUBE_ARB: return GL_TEXTURE_CUBE_MAP_ARB;
    case GL_SAMPLER_2D_RECT_ARB: return GL_TEXTURE_RECTANGLE_ARB;
    case GL_SAMPLER_2D_SHADOW_ARB: return GL_TEXTURE_2D;
    case GL_SAMPLER_2D_RECT_SHADOW_ARB: return GL_TEXTURE_RECTANGLE_ARB;
    default: return type;
  }
}

void
GLShader::readActiveUniforms()
{
  alwaysAssertM(complete, std::string(getName()) + ": GL will throw an error unless the shader is complete");

  uniforms.clear();
  int last_texture_unit = 0;

  GLint max_uniform_name_length = 0, num_active_uniforms = 0;

  // On ATI cards, we are required to call glUseProgramObjectARB before glGetUniformLocationARB,
  // so we always bind the program first.

  // First, store the old program.
  GLhandleARB old_program = glGetHandleARB(GL_PROGRAM_OBJECT_ARB);

  // The shader must be linked before we can use it
  link();
  glUseProgramObjectARB(program_id);
  THEA_CHECK_GL_OK

  try
  {
    // Get the number of uniforms, and the length of the longest name.
    glGetObjectParameterivARB(program_id, GL_OBJECT_ACTIVE_UNIFORM_MAX_LENGTH_ARB, &max_uniform_name_length);
    glGetObjectParameterivARB(program_id, GL_OBJECT_ACTIVE_UNIFORMS_ARB, &num_active_uniforms);

    TheaArray<GLcharARB> name_chars(max_uniform_name_length);

    // Loop over glGetActiveUniformARB and store the results away.
    for (int i = 0; i < num_active_uniforms; ++i)
    {
      UniformData data;

      glGetActiveUniformARB(program_id, i, max_uniform_name_length, NULL, &data.size, &data.type, &name_chars[0]);
      std::string name(&name_chars[0]);

      data.location = glGetUniformLocationARB(program_id, &name_chars[0]);

      bool is_gl_builtin = (data.location == -1 || beginsWith(name, "gl_"));
      if (!is_gl_builtin)
      {
        if (GLShader__isSamplerType(data.type))
        {
          if (data.size > 1)
            throw Error(std::string(getName()) + ": Arrays of samplers are not supported");

          data.texunit = last_texture_unit++;
        }

        uniforms[name] = data;
      }
    }
  }
  catch (...)
  {
    glUseProgramObjectARB(old_program);
    throw;
  }

  glUseProgramObjectARB(old_program);
  THEA_CHECK_GL_OK
}

void
GLShader::bindUniforms()
{
  for (Uniforms::iterator ui = uniforms.begin(); ui != uniforms.end(); ++ui)
  {
    if (!ui->second.has_value)
      throw Error(std::string(getName()) + ": Uniform '" + ui->first + "' has not been set");

    if (ui->second.requires_rebind)
    {
      ui->second.requires_rebind = false;

      switch (ui->second.type)
      {
        case GL_FLOAT: glUniform1fARB(ui->second.location, ui->second.value.f_val); break;
        case GL_INT:   glUniform1iARB(ui->second.location, ui->second.value.i_val); break;
        case GL_BOOL:  glUniform1iARB(ui->second.location, ui->second.value.i_val); break;

        case GL_FLOAT_VEC2_ARB: glUniform2fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC2_ARB:   glUniform2ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC2_ARB:  glUniform2ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_VEC3_ARB: glUniform3fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC3_ARB:   glUniform3ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC3_ARB:  glUniform3ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_VEC4_ARB: glUniform4fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC4_ARB:   glUniform4ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC4_ARB:  glUniform4ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_MAT2_ARB: glUniformMatrix2fvARB(ui->second.location, ui->second.size, GL_TRUE,
                                                      &ui->second.value.f_array[0]); break;
        case GL_FLOAT_MAT3_ARB: glUniformMatrix3fvARB(ui->second.location, ui->second.size, GL_TRUE,
                                                      &ui->second.value.f_array[0]); break;
        case GL_FLOAT_MAT4_ARB: glUniformMatrix4fvARB(ui->second.location, ui->second.size, GL_TRUE,
                                                      &ui->second.value.f_array[0]); break;

        case GL_SAMPLER_1D_ARB:
        case GL_SAMPLER_2D_ARB:
        case GL_SAMPLER_3D_ARB:
        case GL_SAMPLER_CUBE_ARB:
        case GL_SAMPLER_2D_RECT_ARB:
        case GL_SAMPLER_2D_SHADOW_ARB:
        case GL_SAMPLER_2D_RECT_SHADOW_ARB:
        {
          GLenum tex_type = GLShader__toCanonicalType(ui->second.type);
          glEnable(tex_type);

          if (THEA_GL_SUPPORTS(ARB_multitexture))
            glActiveTextureARB(GL_TEXTURE0_ARB + ui->second.texunit);
          else
          {
            debugAssertM(ui->second.texunit == 0,
                         std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
          }

          glBindTexture(tex_type, ui->second.value.texture->getGLID());
          glUniform1iARB(ui->second.location, ui->second.texunit);
          ui->second.requires_rebind = true;  // samplers must always be rebound
          break;
        }

        default:
          throw Error(std::string(getName()) + ": Type of uniform '" + ui->first + "' not handled");
      }
    }
  }
}

void
GLShader::link()
{
  if (linked)
    return;

  glLinkProgramARB(program_id);
  checkBuildStatus(program_id, GL_OBJECT_LINK_STATUS_ARB, "Failed to link shader");
  linked = true;
}

void
GLShader::use()
{
  if (!complete)
    throw Error(std::string(getName()) + ": Shader is incomplete");

  link();
  glUseProgramObjectARB(program_id);
  bindUniforms();
  THEA_CHECK_GL_OK
}

void
GLShader::checkBuildStatus(GLhandleARB obj_id, GLenum status_field, std::string const & error_msg)
{
  int status;
  glGetObjectParameterivARB(obj_id, status_field, &status);
  if (!status)
  {
    int infolog_length;
    glGetObjectParameterivARB(obj_id, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infolog_length);
    if (infolog_length > 0)
    {
      TheaArray<char> infolog(infolog_length);
      GLsizei chars_written;
      glGetInfoLogARB(obj_id, (GLsizei)infolog_length, &chars_written, &infolog[0]);
      throw Error(std::string(getName()) + ": " + error_msg + " (" + &infolog[0] + ')');
    }
    else
      throw Error(std::string(getName()) + ": " + error_msg + " (unknown error)");
  }
}

#define GLShader__SET_UNIFORM_STANDARD_CHECKS(uniform_gl_type, uniform_size)                                                  \
  if (entry == uniforms.end())                                                                                                \
  {                                                                                                                           \
    THEA_WARNING << getName() << ": Uniform '" << uniform_name << "' not found";                                              \
    return;                                                                                                                   \
  }                                                                                                                           \
                                                                                                                              \
  if (GLShader__toCanonicalType(entry->second.type) != uniform_gl_type)                                                       \
    throw Error(std::string(getName()) + ": Argument does not match the declared type of uniform '"                           \
              + std::string(uniform_name) + '\'');                                                                            \
                                                                                                                              \
  if (entry->second.size != (GLint)uniform_size)                                                                              \
    throw Error(std::string(getName()) + format(": Uniform '%s' expects size %d", uniform_name, (int)(uniform_size)));

namespace GLShaderInternal {

static void toFloats(Vector2 const & v, float * f) { f[0] = v.x(); f[1] = v.y(); }
static void toFloats(Vector3 const & v, float * f) { f[0] = v.x(); f[1] = v.y(); f[2] = v.z(); }
static void toFloats(Vector4 const & v, float * f) { f[0] = v.x(); f[1] = v.y(); f[2] = v.z(); f[3] = v.w(); }

static void toFloats(ColorL const & c, float * f) { *f = c.value(); }
static void toFloats(ColorRGB const & c, float * f) { f[0] = c.r(); f[1] = c.g(); f[2] = c.b(); }
static void toFloats(ColorRGBA const & c, float * f) { f[0] = c.r(); f[1] = c.g(); f[2] = c.b(); f[3] = c.a(); }

static void toFloats(Matrix2 const & m, float * f) { m.getElementsRowMajor(f); }
static void toFloats(Matrix3 const & m, float * f) { m.getElementsRowMajor(f); }
static void toFloats(Matrix4 const & m, float * f) { m.getElementsRowMajor(f); }

} // namespace GLShaderInternal

void
GLShader::setUniform(char const * uniform_name, float value)
{
  Uniforms::iterator entry = uniforms.find(uniform_name);
  GLShader__SET_UNIFORM_STANDARD_CHECKS(GL_FLOAT, 1);
  entry->second.value.f_val = value;
  entry->second.valueChanged();
}

void
GLShader::setUniform(char const * uniform_name, int value)
{
  Uniforms::iterator entry = uniforms.find(uniform_name);
  GLShader__SET_UNIFORM_STANDARD_CHECKS(GL_INT, 1);
  entry->second.value.i_val = value;
  entry->second.valueChanged();
}

void
GLShader::setUniform(char const * uniform_name, ColorL const & value)
{
  setUniform(uniform_name, value.value());
}

void
GLShader::setUniform(char const * uniform_name, ColorL8 const & value)
{
  setUniform(uniform_name, ColorL(value).value());
}

void
GLShader::setUniform(char const * uniform_name, Texture * value)
{
  if (!value)
    throw Error(std::string(getName()) + ": Null argument passed for uniform '" + std::string(uniform_name) + '\'');

  GLTexture * gltex = dynamic_cast<GLTexture *>(value);
  debugAssertM(gltex, std::string(getName()) + ": Attempt to use a non-OpenGL texture with a GL rendersystem");

  Uniforms::iterator entry = uniforms.find(uniform_name);
  GLShader__SET_UNIFORM_STANDARD_CHECKS(gltex->getGLTarget(), 1);
  entry->second.value.texture = gltex;
  entry->second.valueChanged();
}

#define GLShader__MULTI_FLOAT_SET_UNIFORM(uniform_type, uniform_convert_type, uniform_gl_type, num_components)                \
  void                                                                                                                        \
  GLShader::setUniform(char const * uniform_name, uniform_type const & value)                                                 \
  {                                                                                                                           \
    Uniforms::iterator entry = uniforms.find(uniform_name);                                                                   \
    GLShader__SET_UNIFORM_STANDARD_CHECKS(uniform_gl_type, 1);                                                                \
    entry->second.value.f_array.resize(num_components);                                                                       \
    GLShaderInternal::toFloats(static_cast<uniform_convert_type>(value), &entry->second.value.f_array[0]);                    \
    entry->second.valueChanged();                                                                                             \
  }

GLShader__MULTI_FLOAT_SET_UNIFORM(Vector2,      Vector2,    GL_FLOAT_VEC2_ARB,  2)
GLShader__MULTI_FLOAT_SET_UNIFORM(Vector3,      Vector3,    GL_FLOAT_VEC3_ARB,  3)
GLShader__MULTI_FLOAT_SET_UNIFORM(Vector4,      Vector4,    GL_FLOAT_VEC4_ARB,  4)
GLShader__MULTI_FLOAT_SET_UNIFORM(ColorRGB,     ColorRGB,   GL_FLOAT_VEC3_ARB,  3)
GLShader__MULTI_FLOAT_SET_UNIFORM(ColorRGB8,    ColorRGB,   GL_FLOAT_VEC3_ARB,  3)
GLShader__MULTI_FLOAT_SET_UNIFORM(ColorRGBA,    ColorRGBA,  GL_FLOAT_VEC4_ARB,  4)
GLShader__MULTI_FLOAT_SET_UNIFORM(ColorRGBA8,   ColorRGBA,  GL_FLOAT_VEC4_ARB,  4)
GLShader__MULTI_FLOAT_SET_UNIFORM(Matrix2,      Matrix2,    GL_FLOAT_MAT2_ARB,  4)
GLShader__MULTI_FLOAT_SET_UNIFORM(Matrix3,      Matrix3,    GL_FLOAT_MAT3_ARB,  9)
GLShader__MULTI_FLOAT_SET_UNIFORM(Matrix4,      Matrix4,    GL_FLOAT_MAT4_ARB, 16)

void
GLShader::setUniform(char const * uniform_name, long num_values, float const * values)
{
  Uniforms::iterator entry = uniforms.find(uniform_name);
  GLShader__SET_UNIFORM_STANDARD_CHECKS(GL_FLOAT, num_values);
  if (entry->second.size < 2)
    throw Error("Attempting to set non-array uniform '" + std::string(uniform_name) + "' from array type");

  entry->second.value.f_array.resize((size_t)num_values);
  std::memcpy(&entry->second.value.f_array[0], values, (size_t)(num_values * sizeof(float)));
  entry->second.valueChanged();
}

void
GLShader::setUniform(char const * uniform_name, long num_values, int const * values)
{
  Uniforms::iterator entry = uniforms.find(uniform_name);
  GLShader__SET_UNIFORM_STANDARD_CHECKS(GL_INT, num_values);
  if (entry->second.size < 2)
    throw Error("Attempting to set non-array uniform '" + std::string(uniform_name) + "' from array type");

  entry->second.value.i_array.resize((size_t)num_values);
  std::memcpy(&entry->second.value.i_array[0], values, (size_t)(num_values * sizeof(int)));
  entry->second.valueChanged();
}

void
GLShader::setUniform(char const * uniform_name, long num_values, Texture * const * values)
{
  throw Error(std::string(getName()) + ": OpenGL texture array uniforms are not supported");
}

#define GLShader__FLOAT_ARRAY_SET_UNIFORM(uniform_type, uniform_convert_type, uniform_gl_type, num_components)                \
  void                                                                                                                        \
  GLShader::setUniform(char const * uniform_name, long num_values, uniform_type const * values)                               \
  {                                                                                                                           \
    Uniforms::iterator entry = uniforms.find(uniform_name);                                                                   \
    GLShader__SET_UNIFORM_STANDARD_CHECKS(uniform_gl_type, num_values);                                                       \
    if (entry->second.size < 2)                                                                                               \
      throw Error("Attempting to set non-array uniform '" + std::string(uniform_name) + "' from array type");                 \
                                                                                                                              \
    entry->second.value.f_array.resize((size_t)(num_components * num_values));                                          \
                                                                                                                              \
    /* Copy elements one by one to avoid packing issues */                                                                    \
    static float * array_start = &entry->second.value.f_array[0];                                                             \
    for (long i = 0; i < num_values; ++i)                                                                                     \
      GLShaderInternal::toFloats(static_cast<uniform_convert_type>(values[i]), array_start + i * num_components);             \
                                                                                                                              \
    entry->second.valueChanged();                                                                                             \
  }

GLShader__FLOAT_ARRAY_SET_UNIFORM(Vector2,      Vector2,    GL_FLOAT_VEC2_ARB,  2)
GLShader__FLOAT_ARRAY_SET_UNIFORM(Vector3,      Vector3,    GL_FLOAT_VEC3_ARB,  3)
GLShader__FLOAT_ARRAY_SET_UNIFORM(Vector4,      Vector4,    GL_FLOAT_VEC4_ARB,  4)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorL,       ColorL,     GL_FLOAT,           1)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorL8,      ColorL,     GL_FLOAT,           1)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorRGB,     ColorRGB,   GL_FLOAT_VEC3_ARB,  3)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorRGB8,    ColorRGB,   GL_FLOAT_VEC3_ARB,  3)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorRGBA,    ColorRGBA,  GL_FLOAT_VEC4_ARB,  4)
GLShader__FLOAT_ARRAY_SET_UNIFORM(ColorRGBA8,   ColorRGBA,  GL_FLOAT_VEC4_ARB,  4)
GLShader__FLOAT_ARRAY_SET_UNIFORM(Matrix2,      Matrix2,    GL_FLOAT_MAT2_ARB,  4)
GLShader__FLOAT_ARRAY_SET_UNIFORM(Matrix3,      Matrix3,    GL_FLOAT_MAT3_ARB,  9)
GLShader__FLOAT_ARRAY_SET_UNIFORM(Matrix4,      Matrix4,    GL_FLOAT_MAT4_ARB, 16)

} // namespace GL
} // namespace Graphics
} // namespace Thea
