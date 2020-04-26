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
    GLTexture(GLRenderSystem * render_system_, char const * name_, int64 width_, int64 height_, int64 depth_,
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

    int64 getWidth() const { return width; }
    int64 getHeight() const { return height; }
    int64 getDepth() const { return depth; }
    Format const * getFormat() const { return format; }
    Dimension getDimension() const { return dimension; }

    void updateImage(AbstractImage const & image, Face face = Face::POS_X);
    void updateSubImage(AbstractImage const & image,
                        int64 src_x, int64 src_y, int64 src_z, int64 src_width, int64 src_height, int64 src_depth,
                        int64 dst_x, int64 dst_y, int64 dst_z, Face face = Face::POS_X);

    void getImage(AbstractImage & image, Face face = Face::POS_X) const;
    void getSubImage(AbstractImage & image, int64 x, int64 y, int64 z,
                     int64 subimage_width, int64 subimage_height, int64 subimage_depth, Face face = Face::POS_X) const;

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
    int64 width;
    int64 height;
    int64 depth;
    Format const * format;
    Dimension dimension;
    GLenum gl_target;
    GLuint gl_id;

}; // class GLTexture

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
