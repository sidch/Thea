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

#include "IPlugin.hpp"
#include "Application.hpp"
#include "DynLib.hpp"
#include "Algorithms/IEigenSolver.hpp"
#include "Algorithms/ILinearSolver.hpp"
#include "Algorithms/INumericalOptimizer.hpp"
#include "Graphics/IRenderSystem.hpp"

namespace Thea {

typedef IPlugin * (*THEA_DLL_START_PLUGIN)(IFactoryRegistry *);
typedef void      (*THEA_DLL_STOP_PLUGIN) (void);

PluginManager::~PluginManager()
{
  unloadAllPlugins();
}

IPlugin *
PluginManager::load(std::string const & path)
{
  if (dynlibs.find(path) != dynlibs.end()) return nullptr;

  DynLib * dynlib = Application::getDynLibManager().load(path);

  THEA_DLL_START_PLUGIN start_func = (THEA_DLL_START_PLUGIN)dynlib->getSymbol("dllStartPlugin");
  IPlugin * plugin = start_func(this);

  if (!plugin)
  {
    Application::getDynLibManager().unload(dynlib);

    THEA_ERROR << "PluginManager: Could not initialize plugin '" << path << '\'';
    return nullptr;
  }

  install(plugin, dynlib, path);

  return plugin;
}

void
PluginManager::unload(std::string const & path)
{
  unload(dynlibs.find(path));
}

void
PluginManager::unload(DynLibMap::iterator lib)
{
  if (lib != dynlibs.end())
  {
    PluginDynLib pd = lib->second;  // the iterator won't be valid after the uninstall() below

    uninstall(pd.plugin);

    THEA_DLL_STOP_PLUGIN stop_func = (THEA_DLL_STOP_PLUGIN)pd.dynlib->getSymbol("dllStopPlugin");
    stop_func();

    Application::getDynLibManager().unload(pd.dynlib);
  }
}

void
PluginManager::startupAllPlugins()
{
  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
    (*pi)->startup();
}

void
PluginManager::shutdownAllPlugins()
{
  for (PluginList::reverse_iterator pi = plugins.rbegin(); pi != plugins.rend(); ++pi)
    (*pi)->shutdown();
}

void
PluginManager::unloadAllPlugins()
{
  shutdownAllPlugins();

  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
  {
    THEA_LOG << "PluginManager: Uninstalling plugin '" << (*pi)->getName() << '\'';
    (*pi)->uninstall();
  }

  plugins.clear();

  for (DynLibMap::iterator di = dynlibs.begin(); di != dynlibs.end(); ++di)
  {
    THEA_DLL_STOP_PLUGIN stop_func = (THEA_DLL_STOP_PLUGIN)di->second.dynlib->getSymbol("dllStopPlugin");
    stop_func();

    Application::getDynLibManager().unload(di->second.dynlib);
  }

  dynlibs.clear();
}

void
PluginManager::unloadAllDylibs()
{
  for (DynLibMap::iterator di = dynlibs.begin(); di != dynlibs.end(); )
  {
    DynLibMap::iterator curr = di++;
    unload(curr);
  }

  dynlibs.clear();  // just in case, though the unload()'s should have taken care of this
}

void
PluginManager::install(IPlugin * plugin)
{
  install(plugin, nullptr, "");
}

void
PluginManager::install(IPlugin * plugin, DynLib * dynlib, std::string const & path)
{
  alwaysAssertM(plugin, "PluginManager: Plugin pointer cannot be null");
  alwaysAssertM(!dynlib || !path.empty(), "PluginManager: Plugin loaded from a dynamic library must have a non-empty path");

  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
    if (*pi == plugin)  // already installed
      return;

  THEA_LOG << "PluginManager: Installing plugin '" << plugin->getName() << '\'';

  plugins.push_back(plugin);
  if (dynlib)
    dynlibs[path] = PluginDynLib(plugin, dynlib);

  plugin->install();
}

void
PluginManager::uninstall(IPlugin * plugin)
{
  if (!plugin) return;

  THEA_LOG << "PluginManager: Uninstalling plugin '" << plugin->getName() << '\'';
  plugin->uninstall();

  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
    if (*pi == plugin)
    {
      plugins.erase(pi);
      break;
    }

  for (DynLibMap::iterator di = dynlibs.begin(); di != dynlibs.end(); ++di)
    if (di->second.plugin == plugin)
    {
      dynlibs.erase(di);
      break;
    }
}

void
PluginManager::addEigenSolverFactory(char const * name, Algorithms::IEigenSolverFactory * factory)
{
  Application::getEigenSolverManager().installFactory(name, factory);
}

void
PluginManager::removeEigenSolverFactory(char const * name)
{
  Application::getEigenSolverManager().uninstallFactory(name);
}

void
PluginManager::addLinearSolverFactory(char const * name, Algorithms::ILinearSolverFactory * factory)
{
  Application::getLinearSolverManager().installFactory(name, factory);
}

void
PluginManager::removeLinearSolverFactory(char const * name)
{
  Application::getLinearSolverManager().uninstallFactory(name);
}

void
PluginManager::addRenderSystemFactory(char const * name, Graphics::IRenderSystemFactory * factory)
{
  Application::getRenderSystemManager().installFactory(name, factory);
}

void
PluginManager::removeRenderSystemFactory(char const * name)
{
  Application::getRenderSystemManager().uninstallFactory(name);
}

void
PluginManager::addNumericalOptimizerFactory(char const * name, Algorithms::INumericalOptimizerFactory * factory)
{
  Application::getNumericalOptimizerManager().installFactory(name, factory);
}

void
PluginManager::removeNumericalOptimizerFactory(char const * name)
{
  Application::getNumericalOptimizerManager().uninstallFactory(name);
}

} // namespace Thea
