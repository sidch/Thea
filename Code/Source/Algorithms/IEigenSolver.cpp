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

#include "IEigenSolver.hpp"

namespace Thea {
namespace Algorithms {

IEigenSolverFactory *
EigenSolverManager::getFactory(std::string const & type)
{
  FactoryMap::const_iterator installed = installed_factories.find(toLower(type));
  if (installed == installed_factories.end())
    throw Error("No factory for eigensolvers of type '" + type + "' is installed");

  return installed->second;
}

bool
EigenSolverManager::installFactory(std::string const & type, IEigenSolverFactory * factory)
{
  debugAssertM(factory, "EigenSolverManager: Null factory cannot be installed");

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
EigenSolverManager::uninstallFactory(std::string const & type)
{
  installed_factories.erase(type);
}

} // namespace Algorithms
} // namespace Thea
