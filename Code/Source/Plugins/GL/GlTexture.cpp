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

#include "GlTexture.hpp"
#include "GlCaps.hpp"
#include "../../Math.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

namespace GlTextureInternal {

#define THEA_GL_TEXTURE_ERROR { throw Error(std::string(getName()) + ": GlTexture: Error constructing texture"); }

static int8
setDefaultUnpackingOptions(int64 row_alignment)
{
  if (row_alignment < 1) { THEA_ERROR << "GlTexture: Row alignment must be positive"; return false; }

  // GL's default values for everything except alignment
  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, row_alignment);

  return true;
}

static int8
setDefaultPackingOptions(int64 row_alignment)
{
  if (row_alignment < 1) { THEA_ERROR << "GlTexture: Row alignment must be positive"; return false; }

  // GL's default values for everything except alignment
  glPixelStorei(GL_PACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_PACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_PACK_ALIGNMENT, row_alignment);

  return true;
}

static GLenum
dimensionToGlTarget(int32 dimension)
{
  switch (dimension)
  {
    case ITexture::Dimension::DIM_1D: return GL_TEXTURE_1D;
    case ITexture::Dimension::DIM_2D: return GL_TEXTURE_2D;
    case ITexture::Dimension::DIM_3D: return GL_TEXTURE_3D;
    case ITexture::Dimension::DIM_RECTANGLE: return GL_TEXTURE_RECTANGLE_ARB;
    case ITexture::Dimension::DIM_CUBE_MAP: return GL_TEXTURE_CUBE_MAP_ARB;
    default: return GL_INVALID_ENUM;
  }
}

} // namespace GlTextureInternal

GLenum
GlTexture::toGlCubeMapFace(int32 face)
{
  switch (face)
  {
    case ITexture::Face::POS_X: return GL_TEXTURE_CUBE_MAP_POSITIVE_X_ARB;
    case ITexture::Face::NEG_X: return GL_TEXTURE_CUBE_MAP_NEGATIVE_X_ARB;
    case ITexture::Face::POS_Y: return GL_TEXTURE_CUBE_MAP_POSITIVE_Y_ARB;
    case ITexture::Face::NEG_Y: return GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_ARB;
    case ITexture::Face::POS_Z: return GL_TEXTURE_CUBE_MAP_POSITIVE_Z_ARB;
    case ITexture::Face::NEG_Z: return GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_ARB;
    default: return GL_INVALID_ENUM;
  }
}

GlTexture::GlTexture(GlRenderSystem * render_system_, char const * name_, int64 width_, int64 height_, int64 depth_,
                     Format const * desired_format, int32 dimension_, Options const * options)
: render_system(render_system_), name(name_), width(width_), height(height_), depth(depth_), dimension(dimension_)
{
  if (!setInternalFormat(nullptr, desired_format)) THEA_GL_TEXTURE_ERROR
  if (!doSanityChecks()) THEA_GL_TEXTURE_ERROR

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  { GlScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GlClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    gl_target = GlTextureInternal::dimensionToGlTarget(dimension);
    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    if (!setOptions(options)) THEA_GL_TEXTURE_ERROR

    switch (dimension)
    {
      case Dimension::DIM_CUBE_MAP:
        if (!glTexImage(nullptr, format, Face::NEG_X)) THEA_GL_TEXTURE_ERROR
        if (!glTexImage(nullptr, format, Face::POS_Y)) THEA_GL_TEXTURE_ERROR
        if (!glTexImage(nullptr, format, Face::NEG_Y)) THEA_GL_TEXTURE_ERROR
        if (!glTexImage(nullptr, format, Face::POS_Z)) THEA_GL_TEXTURE_ERROR
        if (!glTexImage(nullptr, format, Face::NEG_Z)) THEA_GL_TEXTURE_ERROR
        // Fall-through...

      default:
        if (!glTexImage(nullptr, format, Face::POS_X)) THEA_GL_TEXTURE_ERROR  // face argument is ignored if not cube map
    }
  }
}

GlTexture::GlTexture(GlRenderSystem * render_system_, char const * name_, IImage const * image,
                     Format const * desired_format, int32 dimension_, Options const * options)
: render_system(render_system_), name(name_), dimension(dimension_), gl_target(GlTextureInternal::dimensionToGlTarget(dimension))
{
  if (!image || !image->isValid()) throw Error(std::string(getName()) + ": Source image must be valid");
  if (gl_target == GL_INVALID_ENUM) THEA_GL_TEXTURE_ERROR
  if (dimension == Dimension::DIM_CUBE_MAP)
    throw Error(std::string(getName()) + ": This constructor cannot be used to create a cube map");

  Format const * bytes_format = TextureFormat::fromImageType(IImage::Type(image->getType()));
  if (!setInternalFormat(bytes_format, desired_format)) THEA_GL_TEXTURE_ERROR

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  if (!_updateImage(image, Face::POS_X, options)) THEA_GL_TEXTURE_ERROR
}

GlTexture::GlTexture(GlRenderSystem * render_system_, char const * name_, IImage const * images[6],
                     Format const * desired_format, Options const * options)
: render_system(render_system_), name(name_), dimension(Dimension::DIM_CUBE_MAP),
  gl_target(GlTextureInternal::dimensionToGlTarget(dimension))
{
  if (!images[0] || !images[0]->isValid())
    throw Error(std::string(getName()) + ": All source images must be valid");

  if (images[0]->getDepth() != 1)
    throw Error(std::string(getName()) + ": Cube-mapped textures cannot be 3D");

  IImage::Type type = IImage::Type(images[0]->getType());
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
  if (!setInternalFormat(bytes_format, desired_format)) THEA_GL_TEXTURE_ERROR

  if (!doSanityChecks()) THEA_GL_TEXTURE_ERROR

  glGenTextures(1, &gl_id);
  THEA_CHECK_GL_OK

  { GlScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GlClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_CHECK_GL_OK

    if (!setOptions(options)) THEA_GL_TEXTURE_ERROR

    if (!GlTextureInternal::setDefaultUnpackingOptions(images[0]->getRowAlignment())) THEA_GL_TEXTURE_ERROR
    if (!glTexImage(images[0]->getData(), bytes_format, Face::POS_X)) THEA_GL_TEXTURE_ERROR

    glPixelStorei(GL_PACK_ALIGNMENT, images[1]->getRowAlignment());
    if (!glTexImage(images[1]->getData(), bytes_format, Face::NEG_X)) THEA_GL_TEXTURE_ERROR

    glPixelStorei(GL_PACK_ALIGNMENT, images[2]->getRowAlignment());
    if (!glTexImage(images[2]->getData(), bytes_format, Face::POS_Y)) THEA_GL_TEXTURE_ERROR

    glPixelStorei(GL_PACK_ALIGNMENT, images[3]->getRowAlignment());
    if (!glTexImage(images[3]->getData(), bytes_format, Face::NEG_Y)) THEA_GL_TEXTURE_ERROR

    glPixelStorei(GL_PACK_ALIGNMENT, images[4]->getRowAlignment());
    if (!glTexImage(images[4]->getData(), bytes_format, Face::POS_Z)) THEA_GL_TEXTURE_ERROR

    glPixelStorei(GL_PACK_ALIGNMENT, images[5]->getRowAlignment());
    if (!glTexImage(images[5]->getData(), bytes_format, Face::NEG_Z)) THEA_GL_TEXTURE_ERROR
  }
}

GlTexture::~GlTexture()
{
  glDeleteTextures(1, &gl_id);
}

int8
GlTexture::glTexImage(void const * bytes, Format const * bytes_format, int32 face)
{
  switch (gl_target)
  {
    case GL_TEXTURE_1D:
      glTexImage1D(gl_target, 0, format->openGlFormat(), (GLsizei)width, 0, bytes_format->openGlBaseFormat(),
                   bytes_format->openGlDataFormat(), bytes);
      break;

    case GL_TEXTURE_2D:
    case GL_TEXTURE_RECTANGLE_ARB:
      glTexImage2D(gl_target, 0, format->openGlFormat(), (GLsizei)width, (GLsizei)height, 0, bytes_format->openGlBaseFormat(),
                   bytes_format->openGlDataFormat(), bytes);
      break;

    case GL_TEXTURE_3D:
#ifdef GL_VERSION_1_2  // Apple omits the EXT extension after moving 3D textures to core
      glTexImage3D
#else
      glTexImage3DEXT
#endif
        (gl_target, 0, format->openGlFormat(), (GLsizei)width, (GLsizei)height, (GLsizei)depth, 0,
         bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), bytes);
      break;

    default:  // GL_TEXTURE_CUBE_MAP_ARB
      glTexImage2D(toGlCubeMapFace(face), 0, format->openGlFormat(), (GLsizei)width, (GLsizei)height, 0,
                   bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), bytes);
  }
  THEA_GL_OK_OR_RETURN(0)

  return true;
}

int8
GlTexture::setInternalFormat(Format const * bytes_format, Format const * desired_format)
{
  if (desired_format == TextureFormat::AUTO())
  {
    if (!bytes_format)
    {
      THEA_ERROR << getName() << ": Internal format cannot be automatically determined";
      return false;
    }

    format = bytes_format;
  }
  else
    format = desired_format;

  if (GlCaps::hasBug_redBlueMipmapSwap() && format == TextureFormat::RGB8())
    format = TextureFormat::RGBA8();

  if (!GlCaps::supportsTexture(format))
  {
    THEA_ERROR << getName() << ": Texture format not supported";
    return false;
  }

  return true;
}

int8
GlTexture::doSanityChecks()
{
  if (dimension < 0 || dimension >= Dimension::NUM)
  { THEA_ERROR << getName() << ": Invalid dimensionality"; return false; }

  if (dimension == Dimension::DIM_CUBE_MAP && !THEA_GL_SUPPORTS(ARB_texture_cube_map))
  { THEA_ERROR << getName() << ": Cube map textures are not supported"; return false; }

  if (dimension == Dimension::DIM_RECTANGLE && !THEA_GL_SUPPORTS(ARB_texture_rectangle))
  { THEA_ERROR << getName() << ": Rectangular textures are not supported"; return false; }

  if (width < 1 || height < 1 || depth < 1)
  { THEA_ERROR << getName() << ": Texture must be at least one pixel wide in each dimension"; return false; }

  if (depth > 1 && dimension != Dimension::DIM_3D)
  { THEA_ERROR << getName() << ": Only a 3D texture can have depth greater than one pixel"; return false; }

  if (dimension == Dimension::DIM_1D && (height > 1 || depth > 1))
  { THEA_ERROR << getName() << ": A 1D texture cannot have height or depth greater than one pixel"; return false; }

  if (!(Math::isPowerOf2(width) && Math::isPowerOf2(height) && Math::isPowerOf2(depth)) && dimension != Dimension::DIM_RECTANGLE
   && !THEA_GL_SUPPORTS(ARB_texture_non_power_of_two))  // rectangular textures can be npot by the spec
  { THEA_ERROR << getName() << ": Non-power-of-two textures are not supported"; return false; }

  return true;
}

int8
GlTexture::setOptions(Options const * options)
{
  if (!options) options = TextureOptions::defaults();

  if (options->wrapMode() < 0 || options->wrapMode() >= Options::WrapMode::NUM)
  { THEA_ERROR << getName() << ": Invalid wrap mode"; return false; }

  if (options->interpolateMode() < 0 || options->interpolateMode() >= Options::InterpolateMode::NUM)
  { THEA_ERROR << getName() << ": Invalid interpolation mode"; return false; }

  if (options->depthReadMode() < 0 || options->depthReadMode() >= Options::DepthReadMode::NUM)
  { THEA_ERROR << getName() << ": Invalid depth read mode"; return false; }

  if (dimension == Dimension::DIM_RECTANGLE && options->wrapMode() == Options::WrapMode::TILE)  // see GL_ARB_texture_rectangle spec
  { THEA_ERROR << getName() << ": Tiling is not supported for rectangular textures"; return false; }

  GLenum wrap = GL_REPEAT;
  switch (options->wrapMode())
  {
    case Options::WrapMode::CLAMP: wrap = THEA_GL_SUPPORTS(EXT_texture_edge_clamp) ? GL_CLAMP_TO_EDGE : GL_CLAMP; break;
    case Options::WrapMode::TILE: wrap = GL_REPEAT; break;
    case Options::WrapMode::ZERO:
    {
      wrap = THEA_GL_SUPPORTS(ARB_texture_border_clamp) ? GL_CLAMP_TO_BORDER_ARB : GL_CLAMP;
      GLfloat border_color[4] = {0, 0, 0, 1};
      glTexParameterfv(gl_target, GL_TEXTURE_BORDER_COLOR, border_color);
      break;
    }
    default: THEA_ERROR << getName() << ": Unsupported wrap mode"; return false;
  }

  glTexParameteri(gl_target, GL_TEXTURE_WRAP_S, wrap);
  glTexParameteri(gl_target, GL_TEXTURE_WRAP_T, wrap);
  glTexParameteri(gl_target, GL_TEXTURE_WRAP_R, wrap);

  THEA_GL_OK_OR_RETURN(0)

  int8 has_mipmaps = (options->interpolateMode() == Options::InterpolateMode::NEAREST_MIPMAP
                   || options->interpolateMode() == Options::InterpolateMode::BILINEAR_MIPMAP
                   || options->interpolateMode() == Options::InterpolateMode::TRILINEAR);

  if (dimension == Dimension::DIM_RECTANGLE && has_mipmaps)  // see GL_ARB_texture_rectangle spec
  { THEA_ERROR << getName() << ": Mipmapping is not supported for rectangular textures"; return false; }

  switch (options->interpolateMode())
  {
    case Options::InterpolateMode::NEAREST_NO_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      break;

    case Options::InterpolateMode::NEAREST_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
      break;

    case Options::InterpolateMode::BILINEAR_NO_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      break;

    case Options::InterpolateMode::BILINEAR_MIPMAP:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
      break;

    case Options::InterpolateMode::TRILINEAR:
      glTexParameteri(gl_target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(gl_target, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      break;

    default: THEA_ERROR << getName() << ": Unsupported texture interpolation mode"; return false;
  }

  THEA_GL_OK_OR_RETURN(0)

  if (has_mipmaps)
  {
    if (THEA_GL_SUPPORTS(SGIS_generate_mipmap))
      glTexParameteri(gl_target, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
    else if (GL_VERSION_1_4)  // moved to core in 1.4
      glTexParameteri(gl_target, GL_GENERATE_MIPMAP, GL_TRUE);
    else
    {
      THEA_ERROR << getName() << ": Automatic mipmap generation not supported";
      return false;
      // TODO: Implement manual fallback (or replace everything with glGenerateMipmap for new versions of GL).
    }
  }

  THEA_GL_OK_OR_RETURN(0)

  if (THEA_GL_SUPPORTS(ARB_shadow))
  {
    if (options->depthReadMode() == Options::DepthReadMode::NORMAL)
    {
      glTexParameteri(gl_target, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_MODE_ARB, GL_NONE);
    }
    else
    {
      glTexParameteri(gl_target, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_MODE_ARB, GL_COMPARE_R_TO_TEXTURE_ARB);
      glTexParameteri(gl_target, GL_TEXTURE_COMPARE_FUNC_ARB,
                      (options->depthReadMode() == Options::DepthReadMode::LEQUAL) ? GL_LEQUAL : GL_GEQUAL);
    }
  }
  else if (options->depthReadMode() != Options::DepthReadMode::NORMAL)
  { THEA_ERROR << getName() << ": Comparison-based depth read modes are not supported"; return false; }

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlTexture::updateImage(IImage const * image, int32 face)
{
  return _updateImage(image, face, nullptr);
}

int8
GlTexture::_updateImage(IImage const * image, int32 face, Options const * options)
{
  if (!image || !image->isValid())
  { THEA_ERROR << getName() << ": Cannot update texture from null or invalid image"; return false; }

  { GlScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GlClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);

    THEA_GL_OK_OR_RETURN(0)

    Format const * bytes_format = TextureFormat::fromImageType(IImage::Type(image->getType()));
    width  = image->getWidth();
    height = image->getHeight();
    depth  = image->getDepth();

    if (!doSanityChecks()) return false;

    if (options && !setOptions(options)) return false;
    if (!GlTextureInternal::setDefaultUnpackingOptions(image->getRowAlignment())) return false;

    return glTexImage(image->getData(), bytes_format, face);
  }
}

int8
GlTexture::updateSubImage(IImage const * image, int64 dst_x, int64 dst_y, int64 dst_z, int32 face)
{
  return updateSubImage(image, 0, 0, 0, image->getWidth(), image->getHeight(), image->getDepth(), dst_x, dst_y, dst_z, face);
}

int8
GlTexture::updateSubImage(IImage const * image, int64 src_x, int64 src_y, int64 src_z, int64 src_width, int64 src_height,
                          int64 src_depth, int64 dst_x, int64 dst_y, int64 dst_z, int32 face)
{
  if (!image || !image->isValid())
  { THEA_ERROR << getName() << ": Cannot update texture from null or invalid image"; return false; }

  if (src_x < 0 || src_y < 0 || src_z < 0
   || src_x + src_width  > image->getWidth()
   || src_y + src_height > image->getHeight()
   || src_z + src_depth  > image->getDepth())
  { THEA_ERROR << getName() << ": All or part of subimage lies outside source image boundaries"; return false; }

  if (dst_x < 0 || dst_y < 0 || dst_z < 0
   || dst_x + src_width > width || dst_y + src_height > height || dst_z + src_depth > depth)
  { THEA_ERROR << getName() << ": All or part of subimage lies outside texture boundaries"; return false; }

  Format const * bytes_format = TextureFormat::fromImageType(IImage::Type(image->getType()));

  { GlScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GlClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_GL_OK_OR_RETURN(0)

    if (!GlTextureInternal::setDefaultUnpackingOptions(image->getRowAlignment())) return false;
    glPixelStorei(GL_UNPACK_ROW_LENGTH, image->getWidth());
    glPixelStorei(GL_UNPACK_SKIP_ROWS, src_y);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, src_x);

    switch (gl_target)
    {
      case GL_TEXTURE_1D:
        glTexSubImage1D(gl_target, 0, dst_x, (GLsizei)src_width, bytes_format->openGlBaseFormat(),
                        bytes_format->openGlDataFormat(), image->getData());
        break;

      case GL_TEXTURE_2D:
      case GL_TEXTURE_RECTANGLE_ARB:
        glTexSubImage2D(gl_target, 0, dst_x, dst_y, (GLsizei)src_width, (GLsizei)src_height, bytes_format->openGlBaseFormat(),
                        bytes_format->openGlDataFormat(), image->getData());
        break;

      case GL_TEXTURE_3D:
#ifdef GL_VERSION_1_2  // Apple omits the EXT extension after moving 3D textures to core
        glTexSubImage3D
#else
        glTexSubImage3DEXT
#endif
          (gl_target, 0, dst_x, dst_y, dst_z, (GLsizei)src_width, (GLsizei)src_height, (GLsizei)src_depth,
           bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), image->getData());
        break;

      default:  // GL_TEXTURE_CUBE_MAP_ARB
        glTexSubImage2D(toGlCubeMapFace(face), 0, dst_x, dst_y, (GLsizei)src_width, (GLsizei)src_height,
                        bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), image->getData());
    }
    THEA_GL_OK_OR_RETURN(0)
  }

  return true;
}

int8
GlTexture::getImage(IImage * image, int32 face) const
{
  if (!image)
  { THEA_ERROR << getName() << ": Cannot read texture to null image"; return false; }

  { GlScope scope(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
    GlClientScope client_scope(GL_CLIENT_PIXEL_STORE_BIT);

    glEnable(gl_target);
    glBindTexture(gl_target, gl_id);
    THEA_GL_OK_OR_RETURN(0)

    if (depth > 1) { THEA_ERROR << getName() << ": 3D images are not currently supported"; return false; }
    if (!image->resize(image->getType(), width, height)) return false;

    Format const * bytes_format = TextureFormat::fromImageType(IImage::Type(image->getType()), format->isDepth());
    if (!GlTextureInternal::setDefaultPackingOptions(image->getRowAlignment())) return false;

    if (gl_target == GL_TEXTURE_CUBE_MAP_ARB)
    {
      GLenum gl_face = toGlCubeMapFace(face);
      if (gl_face == GL_INVALID_ENUM) return false;

      glGetTexImage(gl_face, 0, bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), image->getData());
    }
    else
      glGetTexImage(gl_target, 0, bytes_format->openGlBaseFormat(), bytes_format->openGlDataFormat(), image->getData());

    THEA_GL_OK_OR_RETURN(0)
  }

  return true;
}

int8
GlTexture::getSubImage(IImage * image, int64 x, int64 y, int64 z, int64 subimage_width, int64 subimage_height,
                       int64 subimage_depth, int32 face) const
{
  // Until GL gets a GetTexSubImage function...
  THEA_ERROR << getName() << ": Reading texture subimages is not supported";
  return false;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
