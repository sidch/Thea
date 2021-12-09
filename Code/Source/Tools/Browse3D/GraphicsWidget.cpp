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
#include "../../Graphics/ITexture.hpp"
#include "../../Application.hpp"
#include "../../FileSystem.hpp"
#include "../../Image.hpp"
#include "../../MatrixWrapper.hpp"

namespace Browse3D {

Graphics::IShader  * GraphicsWidget::shader = nullptr;
Graphics::IShader  * phong_shader = nullptr;
Graphics::IShader  * matcap_shader = nullptr;
Graphics::ITexture * matcap_tex = nullptr;

Vector3 GraphicsWidget::light_dir       =  Vector3(-1, -1, -2);
ColorRgb GraphicsWidget::light_color    =  ColorRgb(1, 1, 1);
ColorRgb GraphicsWidget::ambient_color  =  ColorRgb(1, 1, 1);
bool GraphicsWidget::two_sided          =  true;
bool GraphicsWidget::flat_shaded        =  false;

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
  }

  return setPhongUniforms(render_system);
}

Graphics::IShader *
GraphicsWidget::getMatcapShader(Graphics::IRenderSystem & render_system)
{
  using namespace Graphics;

  if (!matcap_shader
   && !app().options().matcap.empty()  // faster than the next filesystem test
   && FileSystem::fileExists(app().options().matcap))
  {
    matcap_shader = render_system.createShader("Matcap shader");

    matcap_shader->attachModuleFromFile(IShader::ModuleType::VERTEX,
                                       Application::getResourcePath("Materials/PhongVert.glsl").c_str());
    matcap_shader->attachModuleFromFile(IShader::ModuleType::FRAGMENT,
                                       Application::getResourcePath("Materials/MatcapFrag.glsl").c_str());
  }

  return setMatcapUniforms(render_system);
}

Graphics::IShader *
GraphicsWidget::setPhongUniforms(Graphics::IRenderSystem & render_system)
{
  if (!phong_shader) { return nullptr; }

  Vector4 mat = (app().options().no_shading ? Vector4(1, 0, 0, 0) : app().options().material);
  Vector3 v_light_color(light_color.r(), light_color.g(), light_color.b());
  Vector3 v_ambient_color(ambient_color.r(), ambient_color.g(), ambient_color.b());

  phong_shader->setUniform("material", &asLvalue(Math::wrapMatrix(mat)));
  phong_shader->setUniform("light_dir", &asLvalue(Math::wrapMatrix(light_dir)));
  phong_shader->setUniform("light_color", &asLvalue(Math::wrapMatrix(v_light_color)));
  phong_shader->setUniform("ambient_color", &asLvalue(Math::wrapMatrix(v_ambient_color)));
  phong_shader->setUniform("two_sided", two_sided ? (Real)1 : (Real)0);
  phong_shader->setUniform("flat_shaded", flat_shaded ? (Real)1 : (Real)0);

  return phong_shader;
}

Graphics::IShader *
GraphicsWidget::setMatcapUniforms(Graphics::IRenderSystem & render_system)
{
  if (!matcap_shader) { return nullptr; }

  if (!matcap_tex)
  {
    try
    {
      Image matcap_img(app().options().matcap);
      matcap_tex = render_system.createTexture("Matcap", &matcap_img, Graphics::TextureFormat::AUTO(),
                                               Graphics::ITexture::Dimension::DIM_2D);
      if (!matcap_tex) { return nullptr; }
    }
    THEA_CATCH(return nullptr;, ERROR, "Could not load matcap image from file '%s'", app().options().matcap.c_str())
  }

  matcap_shader->setUniform("matcap_tex", matcap_tex);
  matcap_shader->setUniform("two_sided", two_sided ? (Real)1 : (Real)0);
  matcap_shader->setUniform("flat_shaded", flat_shaded ? (Real)1 : (Real)0);

  return matcap_shader;
}

bool
GraphicsWidget::setSurfaceShader(Graphics::IRenderSystem & render_system)
{
  if (!app().options().matcap.empty() && setMatcapShader(render_system))
    return true;

  return setPhongShader(render_system);
}

bool
GraphicsWidget::setPhongShader(Graphics::IRenderSystem & render_system)
{
  shader = getPhongShader(render_system);  // also rebinds uniforms to current values
  if (!shader) { return false; }

  render_system.setShader(shader);

  return true;
}

bool
GraphicsWidget::setMatcapShader(Graphics::IRenderSystem & render_system)
{
  shader = getMatcapShader(render_system);  // also rebinds uniforms to current values
  if (!shader) { return false; }

  render_system.setShader(shader);

  return true;
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
}

void
GraphicsWidget::setTwoSided(bool value)
{
  two_sided = value;
}

void
GraphicsWidget::setFlatShaded(bool value)
{
  flat_shaded = value;
}

} // namespace Browse3D
