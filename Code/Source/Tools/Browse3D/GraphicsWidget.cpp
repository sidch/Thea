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

#include "GraphicsWidget.hpp"
#include "App.hpp"
#include "../../Graphics/IShader.hpp"
#include "../../Application.hpp"
#include "../../MatrixWrapper.hpp"

namespace Browse3D {

Graphics::IShader * GraphicsWidget::shader = nullptr;
Graphics::IShader * phong_shader = nullptr;

Vector3 GraphicsWidget::light_dir       =  Vector3(-1, -1, -2);
ColorRgb GraphicsWidget::light_color    =  ColorRgb(1, 1, 1);
ColorRgb GraphicsWidget::ambient_color  =  ColorRgb(1, 1, 1);
bool GraphicsWidget::two_sided          =  true;

Graphics::IShader *
GraphicsWidget::getPhongShader(Graphics::IRenderSystem & render_system)
{
  using namespace Graphics;

  if (!phong_shader)
  {
    phong_shader = render_system.createShader("Phong shader");

    phong_shader->attachModuleFromFile(IShader::ModuleType::VERTEX,
                                       Application::getResourcePath("Materials/PhongVert.glsl").c_str());
    phong_shader->attachModuleFromFile(IShader::ModuleType::FRAGMENT,
                                       Application::getResourcePath("Materials/PhongFrag.glsl").c_str());

    setLightingUniforms(phong_shader);
    Vector4 mat = (app().options().no_shading ? Vector4(1, 0, 0, 0) : app().options().material);
    phong_shader->setUniform("material", &asLvalue(Math::wrapMatrix(mat)));
  }

  return phong_shader;
}

void
GraphicsWidget::setLightingUniforms(Graphics::IShader * s)
{
  using namespace Graphics;

  if (!s) s = shader;
  if (!s) return;

  Vector3 v_light_color(light_color.r(), light_color.g(), light_color.b());
  Vector3 v_ambient_color(ambient_color.r(), ambient_color.g(), ambient_color.b());

  s->setUniform("light_dir", &asLvalue(Math::wrapMatrix(light_dir)));
  s->setUniform("light_color", &asLvalue(Math::wrapMatrix(v_light_color)));
  s->setUniform("ambient_color", &asLvalue(Math::wrapMatrix(v_ambient_color)));
  s->setUniform("two_sided", two_sided ? (Real)1 : (Real)0);
}

void
GraphicsWidget::setPhongShader(Graphics::IRenderSystem & render_system)
{
  shader = getPhongShader(render_system);

  setLightingUniforms();  // don't worry about multiple calls, let's play safe

  render_system.setShader(shader);
}

Graphics::IShader *
GraphicsWidget::getShader()
{
  return shader;
}

void
GraphicsWidget::setLight(Vector3 const & dir, ColorRgb const & color, ColorRgb const & ambient_color_)
{
  light_dir = dir;
  light_color = color;
  ambient_color = ambient_color_;
  setLightingUniforms();
}

void
GraphicsWidget::setTwoSided(bool value)
{
  two_sided = value;
  setLightingUniforms();
}

} // namespace Browse3D
