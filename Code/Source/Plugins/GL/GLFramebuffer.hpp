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

#ifndef __Thea_Graphics_GL_GLFramebuffer_hpp__
#define __Thea_Graphics_GL_GLFramebuffer_hpp__

#include "../../Graphics/Framebuffer.hpp"
#include "../../Array.hpp"
#include "GLCommon.hpp"
#include "GLHeaders.hpp"
#include "GLTexture.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Forward declarations
class GLRenderSystem;

/** An OpenGL framebuffer. */
class THEA_GL_DLL_LOCAL GLFramebuffer : public Framebuffer
{
  public:
    /** Constructor. */
    GLFramebuffer(GLRenderSystem * render_system_, char const * name_);

    /** Destructor. */
    ~GLFramebuffer();

    /** Get the parent rendersystem. */
    GLRenderSystem * getRenderSystem() const { return render_system; }

    char const * getName() const { return name.c_str(); }

    void attach(AttachmentPoint ap, Texture * texture, Texture::Face face = Texture::Face::POS_X, int64 z_offset = 0);
    void detach(AttachmentPoint ap);
    void detachAll();

    int64 getWidth() const { return width; }
    int64 getHeight() const { return height; }

    /** Get the OpenGL ID of the framebuffer object. */
    GLuint getGLID() const { return gl_fbid; }

    /** Use the framebuffer for rendering. */
    void use();

  private:
    GLRenderSystem * render_system;
    std::string name;
    GLuint gl_fbid;
    GLTexture * attachment_table[AttachmentPoint::MAX_ATTACHMENTS];
    int32 num_attachments;
    Array<GLenum> gl_draw_buffers;
    int64 width;
    int64 height;

}; // class GLFramebuffer

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
