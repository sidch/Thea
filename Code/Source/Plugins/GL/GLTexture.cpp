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

#include "GLTexture.hpp"
#include "GLCaps.hpp"
#include "../../Math.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

static GLenum
GLTexture__dimensionToGLTarget(Texture::Dimension dimension)
{
  switch (dimension)
  {
    case Texture::Dimension::DIM_1D: return GL_TEXTURE_1D;
    case Texture::Dimension::DIM_2D: return GL_TEXTURE_2D;
    case Texture::Dimension::DIM_3D: return GL_TEXTURE_3D;
    case Texture::Dimension::DIM_RECTANGLE: return GL_TEXTURE_RECTANGLE_ARB;
    case Texture::Dimension::DIM_CUBE_MAP: return GL_TEXTURE_CUBE_MAP_ARB;
    default: throw Error("GLTexture: Unsupported dimension");
  }
}

GLenum
GLTexture::toGLCubeMapFace(Texture::Face face)
{
  switch (face)
  {
    case Texture::Face::POS_X: return GL_TEXTURE_CUBE_MAP_POSITIVE_X_ARB;
    case Texture::Face::NEG_X: return GL_TEXTURE_CUBE_MAP_NEGATIVE_X_ARB;
    case Texture::Face::POS_Y: return GL_TEXTURE_CUBE_MAP_POSITIVE_Y_ARB;
    case Texture::Face::NEG_Y: return GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_ARB;
    case Texture::Face::POS_Z: return GL_TEXTURE_CUBE_MAP_POSITIVE_Z_ARB;
    default:                   return GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_ARB;
  }
}

static void
GLTexture__setDefaultUnpackingOptions(int row_alignment)
{
  debugAssertM(row_alignment >= 1, "GLTexture: Row alignment must be positive");

  // GL's default values for everything except alignment
  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, row_alignment);
}

static void
GLTexture__setDefaultPackingOptions(int row_alignment)
{
  debugAssertM(row_alignment >= 1, "GLTexture: Row alignment must be positive");

  // GL's default values for everything except alignment
  glPixelStorei(GL_PACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_PACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_PACK_ALIGNMENT, row_alignment);
}

GLTexture::GLTexture(GLRenderSystem * render_system_, char const * name_, int width_, int height_, int depth_,
                     Format const * desired_format, Dimension dimension_, Options const & options)
: render_system(render_system_), name(name_), width(width_), height(height_), depth(depth_), dimension(dimension_)
{
  setInternalFormat(NULL, desired_format);
  doSanityChecks();

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  { GLScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GLClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    gl_target = GLTexture__dimensionToGLTarget(dimension);
    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    setOptions(options);

    switch (dimension)
    {
      case Dimension::DIM_CUBE_MAP:
        glTexImage(NULL, format, Face::NEG_X);
        glTexImage(NULL, format, Face::POS_Y);
        glTexImage(NULL, format, Face::NEG_Y);
        glTexImage(NULL, format, Face::POS_Z);
        glTexImage(NULL, format, Face::NEG_Z);
        // Fall-through...

      default:
        glTexImage(NULL, format, Face::POS_X);  // face argument is ignored if not cube map
    }
  }
}

GLTexture::GLTexture(GLRenderSystem * render_system_, char const * name_, AbstractImage const & image,
                     Format const * desired_format, Dimension dimension_, Options const & options)
: render_system(render_system_), name(name_), dimension(dimension_), gl_target(GLTexture__dimensionToGLTarget(dimension))
{
  if (dimension == Dimension::DIM_CUBE_MAP)
    throw Error(std::string(getName()) + ": This constructor cannot be used to create a cube map");

  Format const * bytes_format = TextureFormat::fromImageType(image.getType());
  setInternalFormat(bytes_format, desired_format);

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  _updateImage(image, Face::POS_X, &options);
}

GLTexture::GLTexture(GLRenderSystem * render_system_, char const * name_, AbstractImage const * images[6],
                     Format const * desired_format, Options const & options)
: render_system(render_system_), name(name_), dimension(Dimension::DIM_CUBE_MAP),
  gl_target(GLTexture__dimensionToGLTarget(dimension))
{
  if (!images[0] || !images[0]->isValid())
    throw Error(std::string(getName()) + ": All source images must be valid");

  if (images[0]->getDepth() != 1)
    throw Error(std::string(getName()) + ": Cube-mapped textures cannot be 3D");

  AbstractImage::Type type = images[0]->getType();
  width  = images[0]->getWidth();
  height = images[0]->getHeight();
  depth  = images[0]->getDepth();

  for (int i = 1; i < 6; ++i)
  {
    if (!images[i] || !images[i]->isValid())
      throw Error(std::string(getName()) + ": All source images must be valid");

    if (images[i]->getType() != type
     || images[i]->getWidth() != width || images[i]->getHeight() != height || images[i]->getDepth() != depth)
    {
      throw Error(std::string(getName()) + ": All source images must have identical type and dimensions");
    }
  }

  Format const * bytes_format = TextureFormat::fromImageType(type);
  setInternalFormat(bytes_format, desired_format);

  doSanityChecks();

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  { GLScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GLClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    gl_target = GLTexture__dimensionToGLTarget(dimension);
    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    setOptions(options);

    GLTexture__setDefaultUnpackingOptions(images[0]->getRowAlignment());
    glTexImage(images[0]->getData(), bytes_format, Face::POS_X);

    glPixelStorei(GL_PACK_ALIGNMENT, images[1]->getRowAlignment());
    glTexImage(images[1]->getData(), bytes_format, Face::NEG_X);

    glPixelStorei(GL_PACK_ALIGNMENT, images[2]->getRowAlignment());
    glTexImage(images[2]->getData(), bytes_format, Face::POS_Y);

    glPixelStorei(GL_PACK_ALIGNMENT, images[3]->getRowAlignment());
    glTexImage(images[3]->getData(), bytes_format, Face::NEG_Y);

    glPixelStorei(GL_PACK_ALIGNMENT, images[4]->getRowAlignment());
    glTexImage(images[4]->getData(), bytes_format, Face::POS_Z);

    glPixelStorei(GL_PACK_ALIGNMENT, images[5]->getRowAlignment());
    glTexImage(images[5]->getData(), bytes_format, Face::NEG_Z);
  }
}

GLTexture::~GLTexture()
{
  glDeleteTextures(1, &gl_id);
}

void
GLTexture::glTexImage(void const * bytes, Format const * bytes_format, Face face)
{
  switch (gl_target)
  {
    case GL_TEXTURE_1D:
      glTexImage1D(gl_target, 0, format->openGLFormat, width, 0, bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat,
                   bytes);
      break;

    case GL_TEXTURE_2D:
    case GL_TEXTURE_RECTANGLE_ARB:
      glTexImage2D(gl_target, 0, format->openGLFormat, width, height, 0, bytes_format->openGLBaseFormat,
                   bytes_format->openGLDataFormat, bytes);
      break;

    case GL_TEXTURE_3D:
#ifdef GL_VERSION_1_2  // Apple omits the EXT extension after moving 3D textures to core
      glTexImage3D
#else
      glTexImage3DEXT
#endif
        (gl_target, 0, format->openGLFormat, width, height, depth, 0, bytes_format->openGLBaseFormat,
         bytes_format->openGLDataFormat, bytes);
      break;

    default:  // GL_TEXTURE_CUBE_MAP_ARB
      glTexImage2D(toGLCubeMapFace(face), 0, format->openGLFormat, width, height, 0, bytes_format->openGLBaseFormat,
                   bytes_format->openGLDataFormat, bytes);
  }
  THEA_CHECK_GL_OK
}

void
GLTexture::setInternalFormat(Format const * bytes_format, Format const * desired_format)
{
  if (desired_format == Format::AUTO())
  {
    if (!bytes_format)
      throw Error(std::string(getName()) + ": Internal format cannot be automatically determined");

    format = bytes_format;
  }
  else
    format = desired_format;

  if (GLCaps::hasBug_redBlueMipmapSwap() && format == Format::RGB8())
    format = Format::RGBA8();

  if (!GLCaps::supportsTexture(format))
    throw Error(std::string(getName()) + ": Texture format not supported");
}

void
GLTexture::doSanityChecks()
{
  if (dimension == Dimension::DIM_CUBE_MAP && !THEA_GL_SUPPORTS(ARB_texture_cube_map))
    throw Error(std::string(getName()) + ": Cube map textures are not supported");

  if (dimension == Dimension::DIM_RECTANGLE && !THEA_GL_SUPPORTS(ARB_texture_rectangle))
    throw Error(std::string(getName()) + ": Rectangular textures are not supported");

  if (width < 1 || height < 1 || depth < 1)
    throw Error(std::string(getName()) + ": Texture must be at least one pixel wide in each dimension");

  if (depth > 1 && dimension != Dimension::DIM_3D)
    throw Error(std::string(getName()) + ": Only a 3D texture can have depth greater than one pixel");

  if (dimension == Dimension::DIM_1D && (height > 1 || depth > 1))
    throw Error(std::string(getName()) + ": A 1D texture cannot have height or depth greater than one pixel");

  if (!(Math::isPowerOf2(width) && Math::isPowerOf2(height) && Math::isPowerOf2(depth)) && dimension != Dimension::DIM_RECTANGLE
   && !THEA_GL_SUPPORTS(ARB_texture_non_power_of_two))  // rectangular textures can be npot by the spec
    throw Error(std::string(getName()) + ": Non-power-of-two textures are not supported");
}

void
GLTexture::setOptions(Options const & options)
{
  if (dimension == Dimension::DIM_RECTANGLE && options.wrapMode == WrapMode::TILE)
    throw Error(std::string(getName())
              + ": Tiling is not supported for rectangular textures");  // see GL_ARB_texture_rectangle spec

  GLenum wrap = GL_REPEAT;
  switch (options.wrapMode)
  {
    case WrapMode::CLAMP: wrap = THEA_GL_SUPPORTS(EXT_texture_edge_clamp) ? GL_CLAMP_TO_EDGE : GL_CLAMP; break;
    case WrapMode::TILE: wrap = GL_REPEAT; break;
    case WrapMode::ZERO:
    {
      wrap = THEA_GL_SUPPORTS(ARB_texture_border_clamp) ? GL_CLAMP_TO_BORDER_ARB : GL_CLAMP;
      GLfloat border_color[4] = {0, 0, 0, 1};
      glTexParameterfv(gl_target, GL_TEXTURE_BORDER_COLOR, border_color);
      break;
    }
    default: throw Error(std::string(getName()) + ": Unsupported wrap mode");
  }

  glTexParameteri(gl_target, GL_TEXTURE_WRAP_S, wrap);
  glTexParameteri(gl_target, GL_TEXTURE_WRAP_T, wrap);
  glTexParameteri(gl_target, GL_TEXTURE_WRAP_R, wrap);

  THEA_CHECK_GL_OK

  bool has_mipmaps = (options.interpolateMode == InterpolateMode::NEAREST_MIPMAP
                   || options.interpolateMode == InterpolateMode::BILINEAR_MIPMAP
                   || options.interpolateMode == InterpolateMode::TRILINEAR);

  if (dimension == Dimension::DIM_RECTANGLE && has_mipmaps)
    throw Error(std::string(getName())
              + ": Mipmapping is not supported for rectangular textures");  // see GL_ARB_texture_rectangle spec

  switch (options.interpolateMode)
  {
    case InterpolateMode::NEAREST_NO_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      break;

    case InterpolateMode::NEAREST_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
      break;

    case InterpolateMode::BILINEAR_NO_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      break;

    case InterpolateMode::BILINEAR_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      break;

    case InterpolateMode::TRILINEAR:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      break;

    default: throw Error(std::string(getName()) + ": Unsupported texture interpolation mode");
  }

  THEA_CHECK_GL_OK

  if (has_mipmaps)
  {
    if (THEA_GL_SUPPORTS(SGIS_generate_mipmap))
      glTexParameteri(gl_target, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
    else if (GL_VERSION_1_4)  // moved to core in 1.4
      glTexParameteri(gl_target, GL_GENERATE_MIPMAP, GL_TRUE);
    else
    {
      throw Error(std::string(getName()) + ": Automatic mipmap generation not supported");
      // TODO: Implement manual fallback (or replace everything with glGenerateMipmap for new versions of GL).
    }
  }

  THEA_CHECK_GL_OK

  if (THEA_GL_SUPPORTS(ARB_shadow))
  {
    if (options.depthReadMode == DepthReadMode::NORMAL)
    {
      glTexParameteri(gl_target, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_MODE_ARB, GL_NONE);
    }
    else
    {
      glTexParameteri(gl_target, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_FUNC_ARB,
                      (options.depthReadMode == DepthReadMode::LEQUAL) ? GL_LEQUAL : GL_GEQUAL);
    }
  }
  else if (options.depthReadMode != DepthReadMode::NORMAL)
    throw Error(std::string(getName()) + ": Comparison-based depth read modes are not supported");

  THEA_CHECK_GL_OK
}

void
GLTexture::updateImage(AbstractImage const & image, Face face)
{
  _updateImage(image, face, NULL);
}

void
GLTexture::_updateImage(AbstractImage const & image, Face face, Options const * options)
{
  if (!image.isValid())
    throw Error(std::string(getName()) + ": Cannot update texture from invalid image");

  { GLScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GLClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    Format const * bytes_format = TextureFormat::fromImageType(image.getType());
    width  = image.getWidth();
    height = image.getHeight();
    depth  = image.getDepth();

    doSanityChecks();

    if (options) setOptions(*options);
    GLTexture__setDefaultUnpackingOptions(image.getRowAlignment());

    glTexImage(image.getData(), bytes_format, face);
  }
}

void
GLTexture::updateSubImage(AbstractImage const & image, int src_x, int src_y, int src_z, int src_width, int src_height,
                          int src_depth, int dst_x, int dst_y, int dst_z, Face face)
{
  if (!image.isValid())
    throw Error(std::string(getName()) + ": Cannot update texture from invalid image");

  alwaysAssertM(src_x >= 0 && src_y >= 0 && src_z >= 0
             && src_x + src_width <= image.getWidth()
             && src_y + src_height <= image.getHeight()
             && src_z + src_depth <= image.getDepth(),
                std::string(getName()) + ": All or part of subimage lies outside source image boundaries");
  alwaysAssertM(dst_x >= 0 && dst_y >= 0 && dst_z >= 0
             && dst_x + src_width <= width && dst_y + src_height <= height && dst_z + src_depth <= depth,
                std::string(getName()) + ": All or part of subimage lies outside texture boundaries");

  Format const * bytes_format = TextureFormat::fromImageType(image.getType());

  { GLScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GLClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    GLTexture__setDefaultUnpackingOptions(image.getRowAlignment());
    glPixelStorei(GL_UNPACK_ROW_LENGTH, image.getWidth());
    glPixelStorei(GL_UNPACK_SKIP_ROWS, src_y);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, src_x);

    switch (gl_target)
    {
      case GL_TEXTURE_1D:
        glTexSubImage1D(gl_target, 0, dst_x, (GLsizei)src_width, bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat,
                        image.getData());
        break;

      case GL_TEXTURE_2D:
      case GL_TEXTURE_RECTANGLE_ARB:
        glTexSubImage2D(gl_target, 0, dst_x, dst_y, (GLsizei)src_width, (GLsizei)src_height, bytes_format->openGLBaseFormat,
                        bytes_format->openGLDataFormat, image.getData());
        break;

      case GL_TEXTURE_3D:
#ifdef GL_VERSION_1_2  // Apple omits the EXT extension after moving 3D textures to core
        glTexSubImage3D
#else
        glTexSubImage3DEXT
#endif
          (gl_target, 0, dst_x, dst_y, dst_z, (GLsizei)src_width, (GLsizei)src_height, (GLsizei)src_depth,
           bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat, image.getData());
        break;

      default:  // GL_TEXTURE_CUBE_MAP_ARB
        glTexSubImage2D(toGLCubeMapFace(face), 0, dst_x, dst_y, (GLsizei)src_width, (GLsizei)src_height,
                        bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat, image.getData());
    }
    THEA_CHECK_GL_OK
  }
}

void
GLTexture::getImage(AbstractImage & image, Face face) const
{
  { GLScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GLClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    if (depth > 1) throw Error(std::string(getName()) + ": 3D images are not currently supported");
    image.resize(image.getType(), width, height);

    Format const * bytes_format = TextureFormat::fromImageType(image.getType(), format->isDepth());
    GLTexture__setDefaultPackingOptions(image.getRowAlignment());

    if (gl_target == GL_TEXTURE_CUBE_MAP_ARB)
      glGetTexImage(toGLCubeMapFace(face), 0, bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat, image.getData());
    else
      glGetTexImage(gl_target, 0, bytes_format->openGLBaseFormat, bytes_format->openGLDataFormat, image.getData());

    THEA_CHECK_GL_OK
  }
}

void
GLTexture::getSubImage(AbstractImage & image, int x, int y, int z, int subimage_width, int subimage_height, int subimage_depth,
                       Face face) const
{
  // Until GL gets a GetTexSubImage function...
  throw Error(std::string(getName()) + ": Reading texture subimages is not supported");
}

} // namespace GL
} // namespace Graphics
} // namespace Thea
