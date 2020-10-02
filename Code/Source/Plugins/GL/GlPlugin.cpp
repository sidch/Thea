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

#include "GlPlugin.hpp"
#include "GlRenderSystem.hpp"

namespace Thea {

static Graphics::Gl::GlPlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_GL_API IPlugin *
dllStartPlugin(IFactoryRegistry * registry_)
{
  plugin = new Graphics::Gl::GlPlugin(registry_);
  return plugin;
}

/** DLL stop routine. Uninstalls plugin. */
extern "C" THEA_GL_API void
dllStopPlugin()
{
  delete plugin;
}

namespace Graphics {
namespace Gl {

static char const * GL_PLUGIN_NAME        =  "OpenGL RenderSystem";
static char const * GL_RENDERSYSTEM_NAME  =  "OpenGL";

GlPlugin::GlPlugin(IFactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(GL_PLUGIN_NAME) + ": Factory registry must be non-null");
}

GlPlugin::~GlPlugin()
{
  uninstall();
  delete factory;
}

char const *
GlPlugin::getName() const
{
  return GL_PLUGIN_NAME;
}

void
GlPlugin::install()
{}

void
GlPlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new GlRenderSystemFactory;

    registry->addRenderSystemFactory(GL_RENDERSYSTEM_NAME, factory);
    started = true;
  }
}

void
GlPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllRenderSystems();

    registry->removeRenderSystemFactory(GL_RENDERSYSTEM_NAME);
    started = false;
  }
}

void
GlPlugin::uninstall()
{
  shutdown();  // not currently dependent on presence of other plugins

  // Don't delete the factory, we don't want anyone to be able to create a second GL rendersystem by installing the plugin
  // object again (though this shouldn't normally happen). Delete it only when the plugin is properly destroyed.
  //
  // ^^^ This is expected to change when we figure out how to properly cleanup GL state after a rendersystem is destroyed.
}

} // namespace Gl
} // namespace Graphics

} // namespace Thea
