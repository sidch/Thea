//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_GraphicsWidget_hpp__
#define __Browse3D_GraphicsWidget_hpp__

#include "Common.hpp"
#include "../../Graphics/DrawableObject.hpp"

namespace Thea {
namespace Graphics {

class Shader;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

/** A drawable widget. */
class GraphicsWidget : public Graphics::DrawableObject
{
  public:
    THEA_DEF_POINTER_TYPES(GraphicsWidget, shared_ptr, weak_ptr)

    /** Get the bounding box of the model. */
    AxisAlignedBox3 const & getBounds() const { static AxisAlignedBox3 const dummy; return dummy; }

    /** Update the bounding box of the part. */
    void updateBounds() {}

    /** Upload any data needed to draw the part to the GPU. */
    void uploadToGraphicsSystem(Graphics::RenderSystem & render_system) {}

    /** Select a Phong shader for rendering. */
    static void setPhongShader(Graphics::RenderSystem & render_system);

    /** Get the shader currently being used. */
    static Graphics::Shader * getShader();

    /** Set the lighting parameters. */
    static void setLight(Vector3 const & dir, ColorRGB const & color, ColorRGB const & ambient_color_);

    /** Set two-sided lighting on/off. */
    static void setTwoSided(bool value);

    /** Get the direction of incident light. */
    static Vector3 const & getLightDirection() { return light_dir; }

    /** Get the color of incident light. */
    static ColorRGB const & getLightColor() { return light_color; }

    /** Get the color of ambient light. */
    static ColorRGB const & getAmbientColor() { return ambient_color; }

    /** Check if two-sided lighting is on or off. */
    static bool isTwoSided() { return two_sided; }

  private:
    /** Set shader uniforms related to lighting. */
    static void setLightingUniforms(Graphics::Shader * s = NULL);

    /** Get the wrapped Phong shader. */
    static Graphics::Shader * getPhongShader(Graphics::RenderSystem & render_system);

    static Graphics::Shader * shader;
    static Vector3 light_dir;
    static ColorRGB light_color;
    static ColorRGB ambient_color;
    static bool two_sided;

}; // class GraphicsWidget

} // namespace Browse3D

#endif
