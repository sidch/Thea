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

#ifndef __Thea_Graphics_Gl_GlFramebuffer_hpp__
#define __Thea_Graphics_Gl_GlFramebuffer_hpp__

#include "../../Graphics/IFramebuffer.hpp"
#include "../../Array.hpp"
#include "GlCommon.hpp"
#include "GlHeaders.hpp"
#include "GlTexture.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

// Forward declarations
class GlRenderSystem;

/** An OpenGL framebuffer. */
class THEA_GL_DLL_LOCAL GlFramebuffer : public virtual IFramebuffer
{
  public:
    /** Constructor. */
    GlFramebuffer(GlRenderSystem * render_system_, char const * name_);

    /** Destructor. */
    ~GlFramebuffer();

    /** Get the parent rendersystem. */
    GlRenderSystem * getRenderSystem() const { return render_system; }

    char const * THEA_ICALL getName() const { return name.c_str(); }
    int8 THEA_ICALL setName(char const * s) { return false;  /* name is read-only */ }

    int8 THEA_ICALL attach(int32 ap, ITexture * texture, int32 face = ITexture::Face::POS_X, int64 z_offset = 0);
    int8 THEA_ICALL detach(int32 ap);
    int8 THEA_ICALL detachAll();
    int64 THEA_ICALL getWidth() const { return width; }
    int64 THEA_ICALL getHeight() const { return height; }

    /** Get the OpenGL ID of the framebuffer object. */
    GLuint getGlId() const { return gl_fbid; }

    /** Use the framebuffer for rendering. */
    int8 use();

  private:
    GlRenderSystem * render_system;
    std::string name;
    GLuint gl_fbid;
    GlTexture * attachment_table[AttachmentPoint::NUM];
    int32 num_attachments;
    Array<GLenum> gl_draw_buffers;
    int64 width;
    int64 height;

}; // class GlFramebuffer

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
