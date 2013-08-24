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

#include "GraphicsWidget.hpp"
#include "../../Graphics/Shader.hpp"
#include "../../Application.hpp"

namespace Browse3D {

Graphics::Shader * GraphicsWidget::shader;
Vector3 GraphicsWidget::light_dir;
ColorRGB GraphicsWidget::light_color;
ColorRGB GraphicsWidget::ambient_color;

namespace GraphicsWidgetInternal {

Graphics::Shader * phong_shader = NULL;

void
setDefaultPhongUniforms(Graphics::Shader & shader)
{
  shader.setUniform("ambient_color", ColorRGB(1, 0.8f, 0.7f));
  shader.setUniform("light_dir", Vector3(-1, -1, -2));
  shader.setUniform("light_color", ColorRGB(1, 1, 1));

  shader.setUniform("material", Vector4(0.2f, 0.6f, 0.2f, 25));
  if (shader.hasUniform("metallicity"))
    shader.setUniform("metallicity", 0.5f);
}

Graphics::Shader *
getPhongShader(Graphics::RenderSystem & render_system)
{
  using namespace Graphics;

  if (!phong_shader)
  {
    phong_shader = render_system.createShader("Phong shader");

    phong_shader->attachModuleFromFile(Shader::ModuleType::VERTEX,
                                       Application::getFullResourcePath("Materials/PhongVert.glsl"));
    phong_shader->attachModuleFromFile(Shader::ModuleType::FRAGMENT,
                                       Application::getFullResourcePath("Materials/PhongNoSpecularFrag.glsl"));

    setDefaultPhongUniforms(*phong_shader);
  }

  return phong_shader;
}

} // namespace GraphicsWidgetInternal

void
GraphicsWidget::setLightingUniforms()
{
  using namespace Graphics;

  if (!shader) return;

  shader->setUniform("light_dir", light_dir);
  shader->setUniform("light_color", light_color);
  shader->setUniform("ambient_color", ambient_color);
}

void
GraphicsWidget::setPhongShader(Graphics::RenderSystem & render_system)
{
  shader = GraphicsWidgetInternal::getPhongShader(render_system);

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

} // namespace Browse3D
