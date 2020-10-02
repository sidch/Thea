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

#ifndef __Thea_Graphics_Gl_GlTexture_hpp__
#define __Thea_Graphics_Gl_GlTexture_hpp__

#include "../../Graphics/ITexture.hpp"
#include "GlCommon.hpp"
#include "GlHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

// Forward declarations
class GlRenderSystem;

/** An OpenGL texture. */
class THEA_GL_DLL_LOCAL GlTexture : public ITexture
{
  public:
    /** Constructs an empty texture of the specified format and size. */
    GlTexture(GlRenderSystem * render_system_, char const * name_, int64 width_, int64 height_, int64 depth_,
              Format const * desired_format, int32 dimension, Options const * options);

    /** Constructs a texture from a pixel buffer. The dimension argument <em>cannot</em> be DIM_CUBE_MAP. */
    GlTexture(GlRenderSystem * render_system_, char const * name_, IImage const * image, Format const * desired_format,
              int32 dimension, Options const * options);

    /** Constructs a cube-map from six pixel buffers, representing 2D images of identical format and size. */
    GlTexture(GlRenderSystem * render_system_, char const * name_, IImage const * images[6], Format const * desired_format,
              Options const * options);

    /** Destructor. */
    ~GlTexture();

    /** Get the parent rendersystem. */
    GlRenderSystem * getRenderSystem() const { return render_system; }

    char const * THEA_ICALL getName() const { return name.c_str(); }

    int64 THEA_ICALL getWidth() const { return width; }
    int64 THEA_ICALL getHeight() const { return height; }
    int64 THEA_ICALL getDepth() const { return depth; }
    Format const * THEA_ICALL getFormat() const { return format; }
    int32 THEA_ICALL getDimension() const { return dimension; }

    int8 THEA_ICALL updateImage(IImage const * image, int32 face = Face::POS_X);
    int8 THEA_ICALL updateSubImage(IImage const * image, int64 dst_x, int64 dst_y, int64 dst_z = 0, int32 face = Face::POS_X);
    int8 THEA_ICALL updateSubImage(IImage const * image,
                                   int64 src_x, int64 src_y, int64 src_z, int64 src_width, int64 src_height, int64 src_depth,
                                   int64 dst_x, int64 dst_y, int64 dst_z, int32 face = Face::POS_X);

    int8 THEA_ICALL getImage(IImage * image, int32 face = Face::POS_X) const;
    int8 THEA_ICALL getSubImage(IImage * image, int64 x, int64 y, int64 z,
                                int64 subimage_width, int64 subimage_height, int64 subimage_depth, int32 face = Face::POS_X)
                                const;

    /** Get the OpenGL target to which this texture is bound (e.g. GL_TEXTURE_2D). */
    GLenum getGlTarget() const { return gl_target; }

    /** Get the OpenGL ID of the texture. */
    GLuint getGlId() const { return gl_id; }

    /** Convert the label of a texture face to the corresponding GL enum. */
    static GLenum toGlCubeMapFace(int32 face);

  private:
    /** Determine the internal storage format of the texture. */
    int8 setInternalFormat(Format const * bytes_format, Format const * desired_format);

    /** Do a series of checks to detect invalid parameters. */
    int8 doSanityChecks();

    /** Set texture parameters from user-specified options. */
    int8 setOptions(Options const * options);

    /** A quick selection of the appropriate glTexImage... call based on current state. */
    int8 glTexImage(void const * bytes, Format const * bytes_format, int32 face);

    /** Updates the texture image and optionally sets user-specified options while doing so. */
    int8 _updateImage(IImage const * image, int32 face, Options const * options);

    GlRenderSystem * render_system;
    std::string name;
    int64 width;
    int64 height;
    int64 depth;
    Format const * format;
    int32 dimension;
    GLenum gl_target;
    GLuint gl_id;

}; // class GlTexture

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
