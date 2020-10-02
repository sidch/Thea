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

#include "GlFramebuffer.hpp"
#include "GlCaps.hpp"
#include <algorithm>

namespace Thea {
namespace Graphics {
namespace Gl {

namespace GlFramebufferInternal {

static GLenum
apToGLenum(int32 ap)
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

  if (ap < 0 || ap >= (intx)(sizeof(to_gl_enum) / sizeof(GLenum)))
    throw Error("GlFramebuffer: Invalid attachment point");

  return to_gl_enum[ap];
}

static int8
isDrawBuffer(int32 ap)
{
  return (ap != IFramebuffer::AttachmentPoint::DEPTH && ap != IFramebuffer::AttachmentPoint::STENCIL);
}

} // namespace GlFramebufferInternal

GlFramebuffer::GlFramebuffer(GlRenderSystem * render_system_, char const * name_)
: render_system(render_system_), name(name_), gl_fbid(0), num_attachments(0), width(0), height(0)
{
  if (!THEA_GL_SUPPORTS(EXT_framebuffer_object))
    throw Error("OpenGL framebuffer objects are not supported");

  glGenFramebuffersEXT(1, &gl_fbid);
  THEA_CHECK_GL_OK

  for (int32 i = 0; i < AttachmentPoint::NUM; ++i)
    attachment_table[i] = nullptr;
}

GlFramebuffer::~GlFramebuffer()
{
  if (gl_fbid > 0)
    glDeleteFramebuffersEXT(1, &gl_fbid);
}

int8
GlFramebuffer::attach(int32 ap, ITexture * texture, int32 face, int64 z_offset)
{
  if (ap < 0 || ap >= AttachmentPoint::NUM)
  { THEA_ERROR << getName() << ": Invalid attachment point"; return false; }

  GlTexture * gl_texture = dynamic_cast<GlTexture *>(texture);
  if (!(texture && gl_texture) && !(!texture && !gl_texture))
  { THEA_ERROR << getName() << ": Texture is not a valid OpenGL texture"; return false; }

  if (gl_texture)
  {
    if (gl_texture->getDimension() == ITexture::Dimension::DIM_3D && z_offset >= gl_texture->getDepth())
    { THEA_ERROR << getName() << ": Z-offset lies outside depth range of 3D texture"; return false; }

    if (num_attachments > 0 && (gl_texture->getWidth() != width || gl_texture->getHeight() != height))
    { THEA_ERROR << getName() << ": Texture to attach does not have same size as existing attachments"; return false; }
  }

  GlTexture * current_attachment = attachment_table[ap];
  if (gl_texture == current_attachment)
    return true;

  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try  // we need to restore the current framebuffer on failure
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
      GLenum gl_ap = GlFramebufferInternal::apToGLenum(ap);

      // Attach the texture
      switch (gl_texture->getDimension())
      {
        case ITexture::Dimension::DIM_1D:
          glFramebufferTexture1DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGlTarget(), (GLuint)gl_texture->getGlId(),
                                    0);
          break;

        case ITexture::Dimension::DIM_2D:
        case ITexture::Dimension::DIM_RECTANGLE:
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGlTarget(), (GLuint)gl_texture->getGlId(),
                                    0);
          break;

        case ITexture::Dimension::DIM_3D:
          glFramebufferTexture3DEXT(GL_FRAMEBUFFER_EXT, gl_ap, (GLenum)gl_texture->getGlTarget(), (GLuint)gl_texture->getGlId(),
                                    0, z_offset);
          break;

        case ITexture::Dimension::DIM_CUBE_MAP:
          glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, GlTexture::toGlCubeMapFace(face), (GLuint)gl_texture->getGlId(),
                                    0);
        break;

        default: throw Error("Unsupported texture dimensionality");
      }

      THEA_CHECK_GL_OK
      num_attachments++;

      if (num_attachments == 1)  // this is the first attachment
      {
        width  = gl_texture->getWidth();
        height = gl_texture->getHeight();
      }

      if (GlFramebufferInternal::isDrawBuffer(ap))
      {
        if (std::find(gl_draw_buffers.begin(), gl_draw_buffers.end(), gl_ap) == gl_draw_buffers.end())
          gl_draw_buffers.push_back(gl_ap);
      }
    }
    else
    {
      // Remove any current attachment at this point
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GlFramebufferInternal::apToGLenum(ap), GL_TEXTURE_2D, 0, 0);
      THEA_CHECK_GL_OK
      num_attachments--;
    }

    attachment_table[ap] = gl_texture;
  }
  catch (std::exception const & e)
  {
    THEA_ERROR << getName() << ": Could not attach texture to framebuffer (" << e.what() << ')';
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }
  catch (...)
  {
    THEA_ERROR << getName() << ": Could not attach texture to framebuffer";
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_GL_OK_OR_RETURN(0)
  }

  return true;
}

int8
GlFramebuffer::detach(int32 ap)
{
  if (ap < 0 || ap >= AttachmentPoint::NUM)
  { THEA_ERROR << getName() << ": Invalid attachment point"; return false; }

  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try  // we need to restore the current framebuffer on failure
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
      GLenum gl_ap = GlFramebufferInternal::apToGLenum(ap);

      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, gl_ap, GL_TEXTURE_2D, 0, 0);
      THEA_CHECK_GL_OK

      attachment_table[ap] = nullptr;
      num_attachments--;

      if (GlFramebufferInternal::isDrawBuffer(ap))
      {
        Array<GLenum>::iterator loc = std::find(gl_draw_buffers.begin(), gl_draw_buffers.end(), gl_ap);
        if (loc != gl_draw_buffers.end())  // should always be true
          gl_draw_buffers.erase(loc);
      }
    }
  }
  catch (std::exception const & e)
  {
    THEA_ERROR << getName() << ": Could not detach texture from framebuffer (" << e.what() << ')';
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }
  catch (...)
  {
    THEA_ERROR << getName() << ": Could not detach texture from framebuffer";
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_GL_OK_OR_RETURN(0)
  }

  return true;
}

int8
GlFramebuffer::detachAll()
{
  // Get current framebuffer
  GLuint orig_fb = (GLuint)glGetInteger(GL_FRAMEBUFFER_BINDING_EXT);

  try  // we need to restore the current framebuffer on failure
  {
    // If we aren't already bound, bind us now
    if (orig_fb != gl_fbid)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
      THEA_CHECK_GL_OK
    }

    // Remove all current attachments
    for (int32 ap = 0; ap < AttachmentPoint::NUM; ++ap)
    {
      if (attachment_table[ap])
      {
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GlFramebufferInternal::apToGLenum((AttachmentPoint)ap), GL_TEXTURE_2D,
                                  0, 0);
        THEA_CHECK_GL_OK

        attachment_table[ap] = nullptr;
        num_attachments--;
      }
    }

    assert(num_attachments == 0);

    gl_draw_buffers.clear();
  }
  catch (std::exception const & e)
  {
    THEA_ERROR << getName() << ": Could not detach all textures from framebuffer (" << e.what() << ')';
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }
  catch (...)
  {
    THEA_ERROR << getName() << ": Could not detach all textures from framebuffer";
    if (orig_fb != gl_fbid) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);  // restore the original framebuffer
    return false;
  }

  // Restore the original framebuffer
  if (orig_fb != gl_fbid)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, orig_fb);
    THEA_GL_OK_OR_RETURN(0)
  }

  return false;
}

int8
GlFramebuffer::use()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, gl_fbid);
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch (status)
  {
    case GL_FRAMEBUFFER_COMPLETE_EXT: break;

    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
    { THEA_ERROR << getName() << ": Framebuffer is not attachment-complete"; return false; }

    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    { THEA_ERROR << getName() << ": Framebuffer attachments aren't of the same size"; return false; }

    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
    { THEA_ERROR << getName() << ": Color attachments do not all have the same format"; return false; }

    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
    { THEA_ERROR << getName() << ": Draw buffer doesn't have a valid bound texture"; return false; }

    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
    { THEA_ERROR << getName() << ": Read buffer doesn't have a valid bound texture"; return false; }

    default: THEA_ERROR << getName() << ": Unknown/implementation-dependent framebuffer error"; return false;
  }

  glDrawBuffersARB((GLsizei)gl_draw_buffers.size(), &gl_draw_buffers[0]);
  glViewport(0, 0, (GLsizei)width, (GLsizei)height);

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
