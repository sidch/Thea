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

#ifndef __Thea_Graphics_Texture_hpp__
#define __Thea_Graphics_Texture_hpp__

#include "../Common.hpp"
#include "../AbstractImage.hpp"
#include "../NamedObject.hpp"
#include "TextureFormat.hpp"

namespace Thea {
namespace Graphics {

/**
 * An interface for a texture.
 *
 * @todo Make this safe for passing across shared library boundaries.
 */
class THEA_API Texture : public AbstractNamedObject
{
  public:
    /** %Texture dimensionality (enum class). */
    struct THEA_API Dimension
    {
      /** Supported values. */
      enum Value
      {
        DIM_1D,         ///< 1D texture.
        DIM_2D,         ///< 2D texture.
        DIM_3D,         ///< 3D texture.
        DIM_RECTANGLE,  ///< Rectangular texture (deprecated, retained only for compatibility with OpenGL).
        DIM_CUBE_MAP    ///< Cube-mapped texture.
      };

      THEA_ENUM_CLASS_BODY(Dimension)
    };

    /** Wrap modes (enum class). */
    struct THEA_API WrapMode
    {
      /** Supported values. */
      enum Value
      {
        CLAMP,  ///< Clamp out-of-range pixels to the border color.
        TILE,   ///< Tile the image to assign colors to out-of-range pixels.
        ZERO    ///< Out-of-range pixels are all set to zero.
      };

      THEA_ENUM_CLASS_BODY(WrapMode)
    };

    /** Interpolation modes (enum class). */
    struct THEA_API InterpolateMode
    {
      /** Supported values. */
      enum Value
      {
        NEAREST_NO_MIPMAP,   ///< No mipmapping or filtering.
        NEAREST_MIPMAP,      ///< Mipmapping, but no filtering during lookup.
        BILINEAR_NO_MIPMAP,  ///< Bilinear filtering, no mipmapping.
        BILINEAR_MIPMAP,     ///< Bilinear filtering with mipmapping.
        TRILINEAR            ///< Trilinear filtering (bilinear + interpolation between mipmaps)
      };

      THEA_ENUM_CLASS_BODY(InterpolateMode)
    };

    /** Depth read modes, useful for shadow mapping (enum class).*/
    struct THEA_API DepthReadMode
    {
      /** Supported values. */
      enum Value
      {
        NORMAL,  ///< See documentation of GL_TEXTURE_COMPARE_MODE and GL_TEXTURE_COMPARE_FUNC.
        LEQUAL,  ///< See documentation of GL_TEXTURE_COMPARE_MODE and GL_TEXTURE_COMPARE_FUNC.
        GEQUAL   ///< See documentation of GL_TEXTURE_COMPARE_MODE and GL_TEXTURE_COMPARE_FUNC.
      };

      THEA_ENUM_CLASS_BODY(DepthReadMode)
    };

    /** Cube map faces (enum class).*/
    struct THEA_API Face
    {
      /** Supported values. */
      enum Value
      {
        POS_X = 0,  ///< Face with outward normal along positive X.
        NEG_X,      ///< Face with outward normal along negative X.
        POS_Y,      ///< Face with outward normal along negative Y.
        NEG_Y,      ///< Face with outward normal along negative Y.
        POS_Z,      ///< Face with outward normal along negative Z.
        NEG_Z       ///< Face with outward normal along negative Z.
      };

      THEA_ENUM_CLASS_BODY(Face)
    };

    /** %Texture options. */
    struct THEA_API Options
    {
      WrapMode         wrapMode;         ///< Texture wrapping mode.
      InterpolateMode  interpolateMode;  ///< Texture interpolation mode.
      DepthReadMode    depthReadMode;    ///< Depth read mode.

      /** Get the set of default options. */
      static Options const & defaults()
      {
        static Options const ops = { WrapMode::CLAMP, InterpolateMode::BILINEAR_NO_MIPMAP, DepthReadMode::NORMAL };
        return ops;
      }
    };

    /** %Texture storage format. */
    typedef TextureFormat Format;

    /** Destructor. */
    virtual ~Texture() {}

    /** Get the width of the texture in pixels. */
    virtual int64 getWidth() const = 0;

    /** Get the height of the texture in pixels. */
    virtual int64 getHeight() const = 0;

    /** Get the depth of the texture in pixels. */
    virtual int64 getDepth() const = 0;

    /** Get the storage format of the texture. */
    virtual Format const * getFormat() const = 0;

    /** Get the dimensionality of the texture. */
    virtual Dimension getDimension() const = 0;

    /** Update (a face of) the texture from a pixel buffer. The face argument is ignored for non-cube map textures. */
    virtual void updateImage(AbstractImage const & image, Face face = Face::POS_X) = 0;

    /** Update a part of (a face of) the texture from a pixel buffer. The face argument is ignored for non-cube map textures. */
    virtual void updateSubImage(AbstractImage const & image, int64 dst_x, int64 dst_y, int64 dst_z = 0,
                                Face face = Face::POS_X)
    {
      updateSubImage(image, 0, 0, 0, image.getWidth(), image.getHeight(), image.getDepth(), dst_x, dst_y, dst_z, face);
    }

    /**
     * Update a part of (a face of) the texture from a portion of a pixel buffer. The block of the source image with corner
     * (\a src_x, \a src_y, \a src_z) and size \a src_width x \a src_height x \a src_depth is copied to the corresponding block
     * of the texture with corner (\a dst_x, \a dst_y, \a dst_z). The face argument is ignored for non-cube map textures.
     */
    virtual void updateSubImage(AbstractImage const & image,
                                int64 src_x, int64 src_y, int64 src_z, int64 src_width, int64 src_height, int64 src_depth,
                                int64 dst_x, int64 dst_y, int64 dst_z, Face face = Face::POS_X) = 0;

    /** Copy (a face of) the texture into a pixel buffer. The face argument is ignored for non-cube map textures. */
    virtual void getImage(AbstractImage & image, Face face = Face::POS_X) const = 0;

    /** Copy a part of (a face of) the texture into a pixel buffer. The face argument is ignored for non-cube map textures. */
    virtual void getSubImage(AbstractImage & image, int64 x, int64 y, int64 z,
                             int64 subimage_width, int64 subimage_height, int64 subimage_depth,
                             Face face = Face::POS_X) const = 0;

}; // class Texture

} // namespace Graphics
} // namespace Thea

#endif
