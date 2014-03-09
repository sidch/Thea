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

#include "GLFramebuffer.hpp"
#include "GLCaps.hpp"
#include <algorithm>

namespace Thea {
namespace Graphics {
namespace GL {

static GLenum
GLFramebuffer__apToGLenum(Framebuffer::AttachmentPoint ap)
{
  static GLenum const to_gl_enum[] = {
    GL_COLOR_ATTACHMENT0_EXT,
    GL_COLOR_ATTACHMENT1_EXT,
    GL_COLOR_ATTACHMENT2_EXT,
    GL_COLOR_ATTACHMENT3_EXT,
    GL_COLOR_ATTACHMENT4_EXT,
    GL_COLOR_ATTACHMENT5_EXT,
    GL_COLOR_ATTACHMENT6_EXT,
    GL_COLOR_ATTACHMENT7_EXT,
    GL_COLOR_ATTACHMENT8_EXT,
    GL_COLOR_ATTACHMENT9_EXT,
    GL_COLOR_ATTACHMENT10_EXT,
    GL_COLOR_ATTACHMENT11_EXT,
    GL_COLOR_ATTACHMENT12_EXT,
    GL_COLOR_ATTACHMENT13_EXT,
    GL_COLOR_ATTACHMENT14_EXT,
    GL_COLOR_ATTACHMENT15_EXT,
    GL_DEPTH_ATTACHMENT_EXT,
    GL_STENCIL_ATTACHMENT_EXT
  };

  return to_gl_enum[ap];
}

GLFramebuffer::GLFramebuffer(GLRenderSystem * render_system_, char const * name_)
: render_system(render_system_), name(name_), gl_fbid(0), num_attachments(0), width(0), height(0)
{
  if (!THEA_GL_SUPPORTS(EXT_framebuffer_object))
    throw Error("OpenGL framebuffer objects are not supported");

  glGenFramebuffersEXT(1, &gl_fbid);
  THEA_CHECK_GL_OK

  for (int i = 0; i < AttachmentPoint::MAX_ATTACHMENTS; ++i)
    attachment_table[i] = NULL;
}

GLFramebuffer::~GLFramebuffer()
{
  if (gl_fbid > 0)
    glDeleteFramebuffersEXT(1, &gl_fbid);
}

static bool
GLFramebuffer__isDrawBuffer(Framebuffer::AttachmentPoint ap)
{
  return (ap != Framebuffer::AttachmentPoint::DEPTH && ap != Framebuffer::AttachmentPoint::STENCIL);
}

void
GLFramebuffer::attach(AttachmentPoint ap, Texture * texture, Texture::Face face, int z_offset)
{
  GLTexture * gl_texture = dynamic_cast<GLTexture *>(texture);
  debugAssertM((texture && gl_texture) || (!texture && !gl_texture),
               std::string(getName()) + ": Texture is not a valid OpenGL texture");

  if (gl_texture)
  {
    if (gl_texture->getDimension() == Texture::Dimension::DIM_3D && z_offset >= gl_texture->getDepth())
      throw Error(std::string(getName()) + ": Z-offset lies outside depth range of 3D texture");

    if (num_attachments > 0 && (gl_texture->getWidth() != width || gl_texture->getHeight() != height))
      throw Error(std::string(getName()) + ": Texture to attach does not have same size as existing attachments");
  }

  GLTexture * current_attachment = attachment_table[ap];
  if (gl_texture == current_attachment)
    return;

  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try
  {
    // If we aren't already bound, bind us now
    if (orig_fb != gl_fbid)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
      THEA_CHECK_GL_OK
    }

    if (gl_texture)
    {
      // Get the GL attachment point
      GLenum gl_ap = GLFramebuffer__apToGLenum(ap);

      // Attach the texture
      switch (gl_texture->getDimension())
      {
        case Texture::Dimension::DIM_1D:
          glFramebufferTexture1DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGLTarget(), (GLuint)gl_texture->getGLID(),
                                    0);
          break;

        case Texture::Dimension::DIM_2D:
        case Texture::Dimension::DIM_RECTANGLE:
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGLTarget(), (GLuint)gl_texture->getGLID(),
                                    0);
          break;

        case Texture::Dimension::DIM_3D:
          glFramebufferTexture3DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGLTarget(), (GLuint)gl_texture->getGLID(),
                                    0, z_offset);
          break;

        case Texture::Dimension::DIM_CUBE_MAP:
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, GLTexture::toGLCubeMapFace(face), (GLuint)gl_texture->getGLID(),
                                    0);
        break;

        default: throw Error(std::string(getName()) + ": Unsupported texture dimensionality");
      }

      THEA_CHECK_GL_OK
      num_attachments++;

      if (num_attachments == 1)  // this is the first attachment
      {
        width  = gl_texture->getWidth();
        height = gl_texture->getHeight();
      }

      if (GLFramebuffer__isDrawBuffer(ap))
      {
        if (std::find(gl_draw_buffers.begin(), gl_draw_buffers.end(), gl_ap) == gl_draw_buffers.end())
          gl_draw_buffers.push_back(gl_ap);
      }
    }
    else
    {
      // Remove any current attachment at this point
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GLFramebuffer__apToGLenum(ap), GL_TEXTURE_2D, 0, 0);
      THEA_CHECK_GL_OK
      num_attachments--;
    }

    attachment_table[ap] = gl_texture;
  }
  catch (...)
  {
    // Restore the original framebuffer
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    throw;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_CHECK_GL_OK
  }
}

void
GLFramebuffer::detach(AttachmentPoint ap)
{
  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try
  {
    // If we aren't already bound, bind us now
    if (orig_fb != gl_fbid)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
      THEA_CHECK_GL_OK
    }

    if (attachment_table[ap])
    {
      // Get the GL attachment point
      GLenum gl_ap = GLFramebuffer__apToGLenum(ap);

      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, GL_TEXTURE_2D, 0, 0);
      THEA_CHECK_GL_OK

      attachment_table[ap] = NULL;
      num_attachments--;

      if (GLFramebuffer__isDrawBuffer(ap))
      {
        TheaArray<GLenum>::iterator loc = std::find(gl_draw_buffers.begin(), gl_draw_buffers.end(), gl_ap);
        if (loc != gl_draw_buffers.end())  // should always be true
          gl_draw_buffers.erase(loc);
      }
    }
  }
  catch (...)
  {
    // Restore the original framebuffer
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    throw;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_CHECK_GL_OK
  }
}

void
GLFramebuffer::detachAll()
{
  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try
  {
    // If we aren't already bound, bind us now
    if (orig_fb != gl_fbid)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
      THEA_CHECK_GL_OK
    }

    // Remove all current attachments
    for (int ap = 0; ap < AttachmentPoint::MAX_ATTACHMENTS; ++ap)
    {
      if (attachment_table[ap])
      {
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GLFramebuffer__apToGLenum((AttachmentPoint)ap), GL_TEXTURE_2D, 0, 0);
        THEA_CHECK_GL_OK

        attachment_table[ap] = NULL;
        num_attachments--;
      }
    }

    assert(num_attachments == 0);

    gl_draw_buffers.clear();
  }
  catch (...)
  {
    // Restore the original framebuffer
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    throw;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_CHECK_GL_OK
  }
}

void
GLFramebuffer::use()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch (status)
  {
    case GL_FRAMEBUFFER_COMPLETE_EXT: break;

    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
      throw Error(std::string(getName()) + ": Framebuffer is not attachment-complete");

    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
      throw Error(std::string(getName()) + ": Framebuffer attachments aren't of the same size");

    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
      throw Error(std::string(getName()) + ": Color attachments do not all have the same format");

    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
      throw Error(std::string(getName()) + ": Draw buffer doesn't have a valid bound texture");

    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
      throw Error(std::string(getName()) + ": Read buffer doesn't have a valid bound texture");

    default: throw Error(std::string(getName()) + ": Unknown/implementation-dependent framebuffer error");
  }

  glDrawBuffersARB((GLsizei)gl_draw_buffers.size(), &gl_draw_buffers[0]);
  glViewport(0, 0, width, height);
  THEA_CHECK_GL_OK
}

} // namespace GL
} // namespace Graphics
} // namespace Thea
