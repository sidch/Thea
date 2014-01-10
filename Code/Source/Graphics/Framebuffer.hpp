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

#ifndef __Thea_Graphics_Framebuffer_hpp__
#define __Thea_Graphics_Framebuffer_hpp__

#include "../Common.hpp"
#include "../NamedObject.hpp"
#include "Texture.hpp"

namespace Thea {
namespace Graphics {

/** An interface for a framebuffer. */
class THEA_API Framebuffer : public AbstractNamedObject
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
        MAX_ATTACHMENTS  // only for reading the number of elements in the enum
      };

      THEA_ENUM_CLASS_BODY(AttachmentPoint)
    };

    /** Destructor. */
    virtual ~Framebuffer() {}

    /**
     * Attach a render-texture to an attachment point. Specifying a null texture will cause any existing attachment to be
     * removed.
     */
    virtual void attach(AttachmentPoint ap, Texture * texture, Texture::Face face = Texture::Face::POS_X, int z_offset = 0) = 0;

    /** Detach the current attachment at an attachment point. */
    virtual void detach(AttachmentPoint ap) = 0;

    /** Remove all current attachments. */
    virtual void detachAll() = 0;

    /** Get the width of the framebuffer in pixels. */
    virtual int getWidth() const = 0;

    /** Get the height of the framebuffer in pixels. */
    virtual int getHeight() const = 0;

}; // class Framebuffer

} // namespace Graphics
} // namespace Thea

#endif
