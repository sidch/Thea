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

/** An OpenGL framebuffer. */
class THEA_GL_DLL_LOCAL GLFramebuffer : public Framebuffer
{
  public:
    /** Constructor. */
    GLFramebuffer(char const * name_);

    /** Destructor. */
    ~GLFramebuffer();

    char const * getName() const { return name.c_str(); }

    void attach(AttachmentPoint ap, Texture * texture, Texture::Face face = Texture::Face::POS_X, int z_offset = 0);
    void detach(AttachmentPoint ap);
    void detachAll();

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    /** Get the OpenGL ID of the framebuffer object. */
    GLuint getGLID() const { return gl_fbid; }

    /** Use the framebuffer for rendering. */
    void use();

  private:
    std::string name;
    GLuint gl_fbid;
    GLTexture * attachment_table[AttachmentPoint::MAX_ATTACHMENTS];
    int num_attachments;
    TheaArray<GLenum> gl_draw_buffers;
    int width;
    int height;

}; // class GLFramebuffer

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
