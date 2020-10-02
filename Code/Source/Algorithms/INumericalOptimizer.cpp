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

#include "INumericalOptimizer.hpp"

namespace Thea {
namespace Algorithms {

INumericalOptimizerFactory *
NumericalOptimizerManager::getFactory(std::string const & type)
{
  FactoryMap::const_iterator installed = installed_factories.find(toLower(type));
  if (installed == installed_factories.end())
    throw Error("No factory for numerical optimizers of type '" + type + "' is installed");

  return installed->second;
}

bool
NumericalOptimizerManager::installFactory(std::string const & type, INumericalOptimizerFactory * factory)
{
  debugAssertM(factory, "NumericalOptimizerManager: Null factory cannot be installed");

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
NumericalOptimizerManager::uninstallFactory(std::string const & type)
{
  installed_factories.erase(type);
}

} // namespace Algorithms
} // namespace Thea
