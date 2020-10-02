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

#ifndef __Thea_Graphics_ITexture_hpp__
#define __Thea_Graphics_ITexture_hpp__

#include "../Common.hpp"
#include "../IImage.hpp"
#include "../NamedObject.hpp"
#include "TextureFormat.hpp"

namespace Thea {
namespace Graphics {

/** Interface for texture options. */
class THEA_API ITextureOptions
{
  public:
    /** Wrap modes (enum class). */
    struct THEA_API WrapMode
    {
      /** Supported values. */
      enum Value
      {
        CLAMP,  ///< Clamp out-of-range pixels to the border color.
        TILE,   ///< Tile the image to assign colors to out-of-range pixels.
        ZERO,   ///< Out-of-range pixels are all set to zero.
        NUM     ///< [Internal] Number of allowed wrap modes.
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
        TRILINEAR,           ///< Trilinear filtering (bilinear + interpolation between mipmaps)
        NUM                  ///< [Internal] Number of allowed interpolation modes.
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
        GEQUAL,  ///< See documentation of GL_TEXTURE_COMPARE_MODE and GL_TEXTURE_COMPARE_FUNC.
        NUM      ///< [Internal] Number of allowed depth read modes.
      };

      THEA_ENUM_CLASS_BODY(DepthReadMode)
    };

    /** Destructor. */
    virtual ~ITextureOptions() = 0;

    /** Get the texture wrapping mode, as a value from the WrapMode enum. */
    virtual int32 THEA_ICALL wrapMode() const = 0;

    /** Set the texture wrapping mode, as a value from the WrapMode enum. */
    virtual ITextureOptions & THEA_ICALL setWrapMode(int32 mode) = 0;

    /** Get the texture interpolation mode, as a value from the InterpolateMode enum. */
    virtual int32 THEA_ICALL interpolateMode() const = 0;

    /** Set the texture interpolation mode, as a value from the InterpolateMode enum. */
    virtual ITextureOptions & THEA_ICALL setInterpolateMode(int32 mode) = 0;

    /** Get the depth read mode, as a value from the DepthReadMode enum. */
    virtual int32 THEA_ICALL depthReadMode() const = 0;

    /** Set the depth read mode, as a value from the DepthReadMode enum. */
    virtual ITextureOptions & THEA_ICALL setDepthReadMode(int32 mode) = 0;

}; // class ITextureOptions

inline ITextureOptions::~ITextureOptions() {}

/** Texture options. */
class THEA_API TextureOptions : public ITextureOptions
{
  public:
    /** Constructor. Initializes with default options. */
    TextureOptions()
    : wrap_mode(WrapMode::CLAMP),
      interpolate_mode(InterpolateMode::BILINEAR_NO_MIPMAP),
      depth_read_mode(DepthReadMode::NORMAL)
    {}

    int32 THEA_ICALL wrapMode() const { return wrap_mode; }
    ITextureOptions & THEA_ICALL setWrapMode(int32 mode) { wrap_mode = mode; return *this; }
    int32 THEA_ICALL interpolateMode() const { return interpolate_mode; }
    ITextureOptions & THEA_ICALL setInterpolateMode(int32 mode) { interpolate_mode = mode; return *this; }
    int32 THEA_ICALL depthReadMode() const { return depth_read_mode; }
    ITextureOptions & THEA_ICALL setDepthReadMode(int32 mode) { depth_read_mode = mode; return *this; }

    /** Get the set of default options. */
    static TextureOptions const * defaults() { static TextureOptions const def; return &def; }

  private:
    int32 wrap_mode;
    int32 interpolate_mode;
    int32 depth_read_mode;

}; // class TextureOptions

/** Interface for a texture. */
class THEA_API ITexture : public INamedObject
{
  public:
    /** Texture dimensionality (enum class). */
    struct THEA_API Dimension
    {
      /** Supported values. */
      enum Value
      {
        DIM_1D,         ///< 1D texture.
        DIM_2D,         ///< 2D texture.
        DIM_3D,         ///< 3D texture.
        DIM_RECTANGLE,  ///< Rectangular texture (deprecated, retained only for compatibility with OpenGL).
        DIM_CUBE_MAP,   ///< Cube-mapped texture.
        NUM             ///< [Internal] Number of allowed dimensionality values.
      };

      THEA_ENUM_CLASS_BODY(Dimension)
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
        NEG_Z,      ///< Face with outward normal along negative Z.
        NUM         ///< [Internal] Number of allowed faces.
      };

      THEA_ENUM_CLASS_BODY(Face)
    };

    /** Texture storage format. */
    typedef ITextureFormat Format;

    /** Texture options. */
    typedef ITextureOptions Options;

    /** Destructor. */
    virtual ~ITexture() = 0;

    /** Get the width of the texture in pixels. */
    virtual int64 THEA_ICALL getWidth() const = 0;

    /** Get the height of the texture in pixels. */
    virtual int64 THEA_ICALL getHeight() const = 0;

    /** Get the depth of the texture in pixels. */
    virtual int64 THEA_ICALL getDepth() const = 0;

    /** Get the storage format of the texture. */
    virtual Format const * THEA_ICALL getFormat() const = 0;

    /** Get the dimensionality of the texture, as a value from the Dimension enum. */
    virtual int32 THEA_ICALL getDimension() const = 0;

    /**
     * Update (a face of) the texture from a pixel buffer. The face argument is ignored for non-cube map textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL updateImage(IImage const * image, int32 face = Face::POS_X) = 0;

    /**
     * Update a part of (a face of) the texture from a pixel buffer. The face argument is ignored for non-cube map textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL updateSubImage(IImage const * image, int64 dst_x, int64 dst_y, int64 dst_z = 0, int32 face = Face::POS_X) = 0;

    /**
     * Update a part of (a face of) the texture from a portion of a pixel buffer. The block of the source image with corner
     * (\a src_x, \a src_y, \a src_z) and size \a src_width x \a src_height x \a src_depth is copied to the corresponding block
     * of the texture with corner (\a dst_x, \a dst_y, \a dst_z). The face argument is ignored for non-cube map textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL updateSubImage(IImage const * image,
                                int64 src_x, int64 src_y, int64 src_z, int64 src_width, int64 src_height, int64 src_depth,
                                int64 dst_x, int64 dst_y, int64 dst_z, int32 face = Face::POS_X) = 0;

    /**
     * Copy (a face of) the texture into a pixel buffer. The face argument is ignored for non-cube map textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getImage(IImage * image, int32 face = Face::POS_X) const = 0;

    /**
     * Copy a part of (a face of) the texture into a pixel buffer. The face argument is ignored for non-cube map textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getSubImage(IImage * image, int64 x, int64 y, int64 z,
                             int64 subimage_width, int64 subimage_height, int64 subimage_depth,
                             int32 face = Face::POS_X) const = 0;

}; // class ITexture

inline ITexture::~ITexture() {}

} // namespace Graphics
} // namespace Thea

#endif
