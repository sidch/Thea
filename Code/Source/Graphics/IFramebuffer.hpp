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

#ifndef __Thea_Graphics_IFramebuffer_hpp__
#define __Thea_Graphics_IFramebuffer_hpp__

#include "../Common.hpp"
#include "../NamedObject.hpp"
#include "ITexture.hpp"

namespace Thea {
namespace Graphics {

/** Interface for a framebuffer. */
class THEA_API IFramebuffer : public INamedObject
{
  public:
    /** Attachment points (enum class). These are guaranteed to be consecutive integers starting from 0. */
    struct THEA_API AttachmentPoint
    {
      /** Supported values. */
      enum Value
      {
        COLOR_0 = 0,  ///< Color buffer attachment point 0.
        COLOR_1,      ///< Color buffer attachment point 1.
        COLOR_2,      ///< Color buffer attachment point 2.
        COLOR_3,      ///< Color buffer attachment point 3.
        COLOR_4,      ///< Color buffer attachment point 4.
        COLOR_5,      ///< Color buffer attachment point 5.
        COLOR_6,      ///< Color buffer attachment point 6.
        COLOR_7,      ///< Color buffer attachment point 7.
        COLOR_8,      ///< Color buffer attachment point 8.
        COLOR_9,      ///< Color buffer attachment point 9.
        COLOR_10,     ///< Color buffer attachment point 10.
        COLOR_11,     ///< Color buffer attachment point 11.
        COLOR_12,     ///< Color buffer attachment point 12.
        COLOR_13,     ///< Color buffer attachment point 13.
        COLOR_14,     ///< Color buffer attachment point 14.
        COLOR_15,     ///< Color buffer attachment point 15.
        DEPTH,        ///< Depth buffer attachment point.
        STENCIL,      ///< Stencil buffer attachment point.
        NUM           ///< [Internal] Number of allowed attachment points.
      };

      THEA_ENUM_CLASS_BODY(AttachmentPoint)
    };

    /** Destructor. */
    virtual ~IFramebuffer() = 0;

    /**
     * Attach a render-texture to an attachment point. Specifying a null texture will cause any existing attachment to be
     * removed.
     *
     * @param attachment_point A value from the AttachmentPoint enum.
     * @param texture The texture to attach.
     * @param face The face of the texture to attach, as a value from the ITexture::Face enum.
     * @param z_offset The depth slice of the texture to write to, for 3D textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL attach(int32 attachment_point, ITexture * texture, int32 face = ITexture::Face::POS_X, int64 z_offset = 0) = 0;

    /**
     * Detach the current attachment at an attachment point, specified as a value from the AttachmentPoint enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL detach(int32 attachment_point) = 0;

    /**
     * Remove all current attachments.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL detachAll() = 0;

    /** Get the width of the framebuffer in pixels. */
    virtual int64 THEA_ICALL getWidth() const = 0;

    /** Get the height of the framebuffer in pixels. */
    virtual int64 THEA_ICALL getHeight() const = 0;

}; // class IFramebuffer

inline IFramebuffer::~IFramebuffer() {}

} // namespace Graphics
} // namespace Thea

#endif
