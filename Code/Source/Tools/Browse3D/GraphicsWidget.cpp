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
#include "../../Graphics/Shader.hpp"
#include "../../Application.hpp"

namespace Browse3D {

Graphics::Shader * GraphicsWidget::shader = nullptr;
Graphics::Shader * phong_shader = nullptr;

Vector3 GraphicsWidget::light_dir       =  Vector3(-1, -1, -2);
ColorRGB GraphicsWidget::light_color    =  ColorRGB(1, 1, 1);
ColorRGB GraphicsWidget::ambient_color  =  ColorRGB(1, 1, 1);
bool GraphicsWidget::two_sided          =  true;

Graphics::Shader *
GraphicsWidget::getPhongShader(Graphics::RenderSystem & render_system)
{
  using namespace Graphics;

  if (!phong_shader)
  {
    phong_shader = render_system.createShader("Phong shader");

    phong_shader->attachModuleFromFile(Shader::ModuleType::VERTEX,
                                       Application::getResourcePath("Materials/PhongVert.glsl").c_str());
    phong_shader->attachModuleFromFile(Shader::ModuleType::FRAGMENT,
                                       Application::getResourcePath("Materials/PhongFrag.glsl").c_str());

    setLightingUniforms(phong_shader);
    phong_shader->setUniform("material", (app().options().no_shading ? Vector4(1, 0, 0, 0) : app().options().material));
  }

  return phong_shader;
}

void
GraphicsWidget::setLightingUniforms(Graphics::Shader * s)
{
  using namespace Graphics;

  if (!s) s = shader;
  if (!s) return;

  s->setUniform("light_dir", light_dir);
  s->setUniform("light_color", light_color);
  s->setUniform("ambient_color", ambient_color);
  s->setUniform("two_sided", two_sided ? 1.0f : 0.0f);
}

void
GraphicsWidget::setPhongShader(Graphics::RenderSystem & render_system)
{
  shader = getPhongShader(render_system);

  setLightingUniforms();  // don't worry about multiple calls, let's play safe

  render_system.setShader(shader);
}

Graphics::Shader *
GraphicsWidget::getShader()
{
  return shader;
}

void
GraphicsWidget::setLight(Vector3 const & dir, ColorRGB const & color, ColorRGB const & ambient_color_)
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
