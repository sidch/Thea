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

#ifndef __Thea_Graphics_DrawableObject_hpp__
#define __Thea_Graphics_DrawableObject_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBox3.hpp"
#include "RenderOptions.hpp"
#include "RenderSystem.hpp"

namespace Thea {

/** %Graphics functionality. */
namespace Graphics {

/** An object that can be displayed onscreen. */
class THEA_API DrawableObject
{
  public:
    THEA_DEF_POINTER_TYPES(DrawableObject, shared_ptr, weak_ptr)

    /** Destructor. */
    virtual ~DrawableObject() {}

    /** Upload GPU resources to the graphics system. */
    virtual void uploadToGraphicsSystem(RenderSystem & render_system) = 0;

    /** Draw the object using the specified options via the specified rendering system. */
    virtual void draw(RenderSystem & render_system, RenderOptions const & options = RenderOptions::defaults()) const = 0;

    /**
     * Recompute and cache the bounding box for the object. Use this function for objects where the recomputation takes
     * non-trivial time, so that subsequent calls to getBounds() are efficient (and accurate as long as the object does not
     * change).
     */
    virtual void updateBounds() = 0;

    /**
     * Get the bounding box for the object. This should return quickly, using a cached value computed by updateBounds() if
     * necessary.
     */
    virtual AxisAlignedBox3 const & getBounds() const = 0;

}; // class DrawableObject

} // namespace Graphics
} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Graphics::DrawableObject)

#endif
