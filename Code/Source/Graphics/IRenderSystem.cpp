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

#include "IRenderSystem.hpp"

namespace Thea {
namespace Graphics {

IRenderSystemFactory *
RenderSystemManager::getFactory(std::string const & type)
{
  FactoryMap::const_iterator installed = installed_factories.find(toLower(type));
  if (installed == installed_factories.end())
    throw Error("No factory for rendersystems of type '" + type + "' is installed");

  return installed->second;
}

bool
RenderSystemManager::installFactory(std::string const & type, IRenderSystemFactory * factory)
{
  debugAssertM(factory, "RenderSystemManager: Null factory cannot be installed");

  std::string type_lc = toLower(type);
  FactoryMap::const_iterator installed = installed_factories.find(type_lc);
  if (installed == installed_factories.end())
  {
    installed_factories[type_lc] = factory;
    return true;
  }
  else
    return false;
}

void
RenderSystemManager::uninstallFactory(std::string const & type)
{
  installed_factories.erase(type);
}

} // namespace Graphics
} // namespace Thea
