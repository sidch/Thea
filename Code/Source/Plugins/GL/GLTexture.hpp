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

#ifndef __Thea_Graphics_GL_GLTexture_hpp__
#define __Thea_Graphics_GL_GLTexture_hpp__

#include "../../Graphics/Texture.hpp"
#include "GLCommon.hpp"
#include "GLHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Forward declarations
class GLRenderSystem;

/** An OpenGL texture. */
class THEA_GL_DLL_LOCAL GLTexture : public Texture
{
  public:
    /** Constructs an empty texture of the specified format and size. */
    GLTexture(GLRenderSystem * render_system_, char const * name_, int width_, int height_, int depth_,
              Format const * desired_format, Dimension dimension, Options const & options);

    /** Constructs a texture from a pixel buffer. The dimension argument <em>cannot</em> be DIM_CUBE_MAP. */
    GLTexture(GLRenderSystem * render_system_, char const * name_, AbstractImage const & image, Format const * desired_format,
              Dimension dimension, Options const & options);

    /** Constructs a cube-map from six pixel buffers, representing 2D images of identical format and size. */
    GLTexture(GLRenderSystem * render_system_, char const * name_, AbstractImage const * images[6],
              Format const * desired_format, Options const & options);

    /** Destructor. */
    ~GLTexture();

    /** Get the parent rendersystem. */
    GLRenderSystem * getRenderSystem() const { return render_system; }

    char const * getName() const { return name.c_str(); }

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    int getDepth() const { return depth; }
    Format const * getFormat() const { return format; }
    Dimension getDimension() const { return dimension; }

    void updateImage(AbstractImage const & image, Face face = Face::POS_X);
    void updateSubImage(AbstractImage const & image, int src_x, int src_y, int src_width, int src_height, int dst_x, int dst_y,
                        int dst_z = 0, Face face = Face::POS_X);

    void getImage(AbstractImage & image, Face face = Face::POS_X) const;
    void getSubImage(AbstractImage & image, int x, int y, int z, int subimage_width, int subimage_height, int subimage_depth,
                     Face face = Face::POS_X) const;

    /** Get the OpenGL target to which this texture is bound (e.g. GL_TEXTURE_2D). */
    GLenum getGLTarget() const { return gl_target; }

    /** Get the OpenGL ID of the texture. */
    GLuint getGLID() const { return gl_id; }

    /** Convert the label of a texture face to the corresponding GL enum. */
    static GLenum toGLCubeMapFace(Texture::Face face);

  private:
    /** Determine the internal storage format of the texture. */
    void setInternalFormat(Format const * bytes_format, Format const * desired_format);

    /** Do a series of checks to detect invalid parameters. */
    void doSanityChecks();

    /** Set texture parameters from user-specified options. */
    void setOptions(Options const & options);

    /** A quick selection of the appropriate glTexImage... call based on current state. */
    void glTexImage(void const * bytes, Format const * bytes_format, Face face);

    /** Updates the texture image and optionally sets user-specified options while doing so. */
    void _updateImage(AbstractImage const & image, Face face, Options const * options);

    GLRenderSystem * render_system;
    std::string name;
    int width;
    int height;
    int depth;
    Format const * format;
    Dimension dimension;
    GLenum gl_target;
    GLuint gl_id;

}; // class GLTexture

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
