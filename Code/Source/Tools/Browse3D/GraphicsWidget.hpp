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
// First version: 2011
//
//============================================================================

#ifndef __Browse3D_GraphicsWidget_hpp__
#define __Browse3D_GraphicsWidget_hpp__

#include "Common.hpp"
#include "../../AxisAlignedBox3.hpp"
#include "../../Graphics/IDrawable.hpp"

namespace Thea {
namespace Graphics {

class IShader;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

/** A drawable widget. */
class GraphicsWidget : public Graphics::IDrawable
{
  public:
    THEA_DECL_SMART_POINTERS(GraphicsWidget)

    /** Get the bounding box of the model. */
    virtual AxisAlignedBox3 const & getBounds() const { static AxisAlignedBox3 const dummy; return dummy; }

    /** Update the bounding box of the part. */
    virtual void updateBounds() {}

    /** Select an appropriate shader for surface rendering. Prefers matcap over Phong if both are available. */
    static bool setSurfaceShader(Graphics::IRenderSystem & render_system);

    /** Select the Phong shader for rendering. */
    static bool setPhongShader(Graphics::IRenderSystem & render_system);

    /** Select the matcap shader for rendering. */
    static bool setMatcapShader(Graphics::IRenderSystem & render_system);

    /** Get the shader currently being used. */
    static Graphics::IShader * getShader();

    /** Set the lighting parameters for shaders that support these settings. */
    static void setLight(Vector3 const & dir, ColorRgb const & color, ColorRgb const & ambient_color_);

    /** Set two-sided lighting on/off, for shaders that support these settings. */
    static void setTwoSided(bool value);

    /** Set flat shading on/off. */
    static void setFlatShaded(bool value);

    /** Get the direction of incident light. */
    static Vector3 const & getLightDirection() { return light_dir; }

    /** Get the color of incident light. */
    static ColorRgb const & getLightColor() { return light_color; }

    /** Get the color of ambient light. */
    static ColorRgb const & getAmbientColor() { return ambient_color; }

    /** Check if two-sided lighting is on or off. */
    static bool isTwoSided() { return two_sided; }

    /** Check if flat shading is on or off. */
    static bool isFlatShaded() { return flat_shaded; }

  private:
    /** Get the wrapped Phong shader, or a null pointer if it could not be initialized. */
    static Graphics::IShader * getPhongShader(Graphics::IRenderSystem & render_system);

    /** Get the wrapped matcap shader, or a null pointer if it could not be initialized. */
    static Graphics::IShader * getMatcapShader(Graphics::IRenderSystem & render_system);

    /** Set shader uniforms related to Phong shading (including the surface material). Returns null on error. */
    static Graphics::IShader * setPhongUniforms(Graphics::IRenderSystem & render_system);

    /** Set shader uniforms related to matcap materials (including the matcap texture itself). Returns null on error. */
    static Graphics::IShader * setMatcapUniforms(Graphics::IRenderSystem & render_system);

    static Graphics::IShader * shader;
    static Vector3 light_dir;
    static ColorRgb light_color;
    static ColorRgb ambient_color;
    static bool two_sided;
    static bool flat_shaded;

}; // class GraphicsWidget

} // namespace Browse3D

#endif
