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

#include "GlShader.hpp"
#include "GlCaps.hpp"
#include "../../FileSystem.hpp"
#include <cstring>
#include <fstream>

namespace Thea {
namespace Graphics {
namespace Gl {

GlShader::GlShader(GlRenderSystem * render_system_, char const * name_)
: render_system(render_system_), name(name_), complete(false), linked(true), has_vertex_module(false),
  has_fragment_module(false)
{
  if (!THEA_GL_SUPPORTS(ARB_shader_objects))
    throw Error("OpenGL programmable shaders are not supported");

  program_id = glCreateProgramObjectARB();
  THEA_CHECK_GL_OK
}

GlShader::~GlShader()
{
  glDeleteObjectARB(program_id);
}

int8
GlShader::attachModuleFromFile(int32 type, char const * path)
{
  std::ifstream in(path);
  if (!in) { THEA_ERROR << getName() << ": Shader module '" << path << "' not found"; return false; }

  std::string source;
  if (!FileSystem::readWholeFile(path, source))
    return false;

  return attachModuleFromString(type, source.c_str());
}

int8
GlShader::attachModuleFromString(int32 type, char const * source)
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
      { THEA_ERROR << getName() << ": This OpenGL installation does not support geometry shaders"; return false; }

      break;
    }
    default:
      THEA_ERROR << getName() << ": Unknown module type";
      return false;
  }

  GLhandleARB module_id = glCreateShaderObjectARB(gltype);
  THEA_GL_OK_OR_RETURN(0)

  int32 src_length = (int32)std::strlen(source);
  GLcharARB const * src_ptr = (GLcharARB const *)source;
  glShaderSourceARB(module_id, 1, &src_ptr, &src_length);
  THEA_GL_OK_OR_RETURN(0)

  glCompileShaderARB(module_id);
  if (!checkBuildStatus(module_id, GL_OBJECT_COMPILE_STATUS_ARB, "Failed to compile shader module")) return false;

  glAttachObjectARB(program_id, module_id);
  THEA_GL_OK_OR_RETURN(0)

  if (type == ModuleType::VERTEX)   has_vertex_module   = true;
  if (type == ModuleType::FRAGMENT) has_fragment_module = true;
  if (has_vertex_module && has_fragment_module) complete = true;

  // Update the list of active shader uniforms
  if (complete && !readActiveUniforms())
    return false;

  linked = false;

  return true;
}

namespace GlShaderInternal {

static bool
isSamplerType(GLenum type)
{
  return type == GL_SAMPLER_1D_ARB
      || type == GL_SAMPLER_2D_ARB
      || type == GL_SAMPLER_3D_ARB
      || type == GL_SAMPLER_CUBE_ARB
      || type == GL_SAMPLER_2D_RECT_ARB
      || type == GL_SAMPLER_2D_SHADOW_ARB
      || type == GL_SAMPLER_2D_RECT_SHADOW_ARB;
}

} // namespace GlShaderInternal

int8
GlShader::readActiveUniforms()
{
  if (!complete)
  { THEA_ERROR << getName() << ": GL will throw an error unless the shader is complete"; return false; }

  uniforms.clear();
  int32 last_texture_unit = 0;

  GLint max_uniform_name_length = 0, num_active_uniforms = 0;

  // On ATI cards, we are required to call glUseProgramObjectARB before glGetUniformLocationARB,
  // so we always bind the program first.

  // First, store the old program.
  GLhandleARB old_program = glGetHandleARB(GL_PROGRAM_OBJECT_ARB);

  // The shader must be linked before we can use it
  if (!link()) return false;
  glUseProgramObjectARB(program_id);
  THEA_GL_OK_OR_RETURN(0)

  try
  {
    // Get the number of uniforms, and the length of the longest name.
    glGetObjectParameterivARB(program_id, GL_OBJECT_ACTIVE_UNIFORM_MAX_LENGTH_ARB, &max_uniform_name_length);
    glGetObjectParameterivARB(program_id, GL_OBJECT_ACTIVE_UNIFORMS_ARB, &num_active_uniforms);

    Array<GLcharARB> name_chars(max_uniform_name_length);

    // Loop over glGetActiveUniformARB and store the results away.
    for (int32 i = 0; i < num_active_uniforms; ++i)
    {
      UniformData data;

      glGetActiveUniformARB(program_id, i, max_uniform_name_length, nullptr, &data.size, &data.type, &name_chars[0]);
      std::string name(&name_chars[0]);

      data.location = glGetUniformLocationARB(program_id, &name_chars[0]);

      bool is_gl_builtin = (data.location == -1 || beginsWith(name, "gl_"));
      if (!is_gl_builtin)
      {
        if (GlShaderInternal::isSamplerType(data.type))
        {
          if (data.size > 1)
            throw Error(std::string(getName()) + ": Arrays of samplers are not supported");

          data.texunit = last_texture_unit++;
        }

        uniforms[name] = data;
      }
    }
  }
  THEA_STANDARD_CATCH_BLOCKS({ glUseProgramObjectARB(old_program); return false; },
                             ERROR, "%s: Error reading active shader uniforms", getName())

  glUseProgramObjectARB(old_program);
  THEA_GL_OK_OR_RETURN(0)

  return true;
}

int8
GlShader::bindUniforms()
{
  for (Uniforms::iterator ui = uniforms.begin(); ui != uniforms.end(); ++ui)
  {
    if (!ui->second.has_value)
    { THEA_ERROR << getName() << ": Uniform '" << ui->first << "' has not been set"; return false; }

    if (ui->second.requires_rebind)
    {
      ui->second.requires_rebind = false;

      switch (ui->second.type)
      {
        case GL_FLOAT: glUniform1fARB(ui->second.location, ui->second.value.f_value); break;
        case GL_INT:   glUniform1iARB(ui->second.location, ui->second.value.i_value); break;
        case GL_BOOL:  glUniform1iARB(ui->second.location, ui->second.value.i_value); break;

        case GL_FLOAT_VEC2_ARB: glUniform2fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC2_ARB:   glUniform2ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC2_ARB:  glUniform2ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_VEC3_ARB: glUniform3fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC3_ARB:   glUniform3ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC3_ARB:  glUniform3ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_VEC4_ARB: glUniform4fvARB(ui->second.location, ui->second.size, &ui->second.value.f_array[0]); break;
        case GL_INT_VEC4_ARB:   glUniform4ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;
        case GL_BOOL_VEC4_ARB:  glUniform4ivARB(ui->second.location, ui->second.size, &ui->second.value.i_array[0]); break;

        case GL_FLOAT_MAT2_ARB: glUniformMatrix2fvARB(ui->second.location, ui->second.size, GL_FALSE,
                                                      &ui->second.value.f_array[0]); break;
        case GL_FLOAT_MAT3_ARB: glUniformMatrix3fvARB(ui->second.location, ui->second.size, GL_FALSE,
                                                      &ui->second.value.f_array[0]); break;
        case GL_FLOAT_MAT4_ARB: glUniformMatrix4fvARB(ui->second.location, ui->second.size, GL_FALSE,
                                                      &ui->second.value.f_array[0]); break;

        case GL_SAMPLER_1D_ARB:
        case GL_SAMPLER_2D_ARB:
        case GL_SAMPLER_3D_ARB:
        case GL_SAMPLER_CUBE_ARB:
        case GL_SAMPLER_2D_RECT_ARB:
        case GL_SAMPLER_2D_SHADOW_ARB:
        case GL_SAMPLER_2D_RECT_SHADOW_ARB:
        {
          GLenum tex_type = GlShaderInternal::toCanonicalType(ui->second.type);
          glEnable(tex_type);

          if (THEA_GL_SUPPORTS(ARB_multitexture))
            glActiveTextureARB(GL_TEXTURE0_ARB + ui->second.texunit);
          else if (ui->second.texunit != 0)
          { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

          glBindTexture(tex_type, ui->second.value.texture->getGlId());
          glUniform1iARB(ui->second.location, ui->second.texunit);
          ui->second.requires_rebind = true;  // samplers must always be rebound
          break;
        }

        default:
          THEA_ERROR << getName() << ": Type of uniform '" << ui->first << "' not handled";
          return false;
      }
    }
  }

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlShader::link()
{
  if (linked)
    return true;

  glLinkProgramARB(program_id);
  if (!checkBuildStatus(program_id, GL_OBJECT_LINK_STATUS_ARB, "Failed to link shader")) return false;
  linked = true;

  return true;
}

int8
GlShader::use()
{
  if (!complete)
  { THEA_ERROR << getName() << ": Shader is incomplete"; return false; }

  if (!link()) return false;
  glUseProgramObjectARB(program_id);
  if (!bindUniforms()) return false;

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlShader::checkBuildStatus(GLhandleARB obj_id, GLenum status_field, std::string const & error_msg)
{
  int32 status;
  glGetObjectParameterivARB(obj_id, status_field, &status);
  if (!status)
  {
    int32 infolog_length;
    glGetObjectParameterivARB(obj_id, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infolog_length);
    if (infolog_length > 0)
    {
      Array<char> infolog(infolog_length);
      GLsizei chars_written;
      glGetInfoLogARB(obj_id, (GLsizei)infolog_length, &chars_written, &infolog[0]);
      THEA_ERROR << getName() << ": " << error_msg << " (" << &infolog[0] << ')';
      return false;
    }
    else
    {
      THEA_ERROR << getName() << ": " << error_msg << " (unknown error)";
      return false;
    }
  }

  return true;
}

int8
GlShader::setUniform(char const * uniform_name, ITexture * value)
{
  if (!value)
  { THEA_ERROR << getName() << ": Texture for uniform '" << uniform_name << "' must be non-null"; return false; }

  GlTexture * gltex = dynamic_cast<GlTexture *>(value);
  if (!gltex)
  { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL texture with a GL rendersystem"; return false; }

  Uniforms::iterator entry = uniforms.find(uniform_name);
  THEA_GLSHADER_SET_UNIFORM_CHECKS(gltex->getGlTarget(), 1);
  entry->second.value.texture = gltex;
  entry->second.valueChanged();

  return true;
}

int8
GlShader::setUniform(char const * uniform_name, int64 num_values, ITexture * const * values)
{
  THEA_ERROR << getName() << ": OpenGL texture array uniforms are not supported";
  return false;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
