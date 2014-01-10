//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#include "Plugin.hpp"
#include "Application.hpp"
#include "DynLib.hpp"
#include "Algorithms/EigenSolver.hpp"
#include "Algorithms/LinearSolver.hpp"
#include "Algorithms/NumericalOptimizer.hpp"
#include "Graphics/RenderSystem.hpp"

namespace Thea {

typedef Plugin * (*THEA_DLL_START_PLUGIN)(FactoryRegistry *);
typedef void     (*THEA_DLL_STOP_PLUGIN) (void);

PluginManager::~PluginManager()
{
  unloadAllPlugins();
}

Plugin *
PluginManager::load(std::string const & path)
{
  if (dynlibs.find(path) != dynlibs.end()) return NULL;

  DynLib * dynlib = Application::getDynLibManager().load(path);

  THEA_DLL_START_PLUGIN start_func = (THEA_DLL_START_PLUGIN)dynlib->getSymbol("dllStartPlugin");
  Plugin * plugin = start_func(this);

  if (!plugin)
  {
    Application::getDynLibManager().unload(dynlib);

    THEA_ERROR << "PluginManager: Could not initialize plugin '" << path << '\'';
    return NULL;
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
PluginManager::install(Plugin * plugin)
{
  install(plugin, NULL, "");
}

void
PluginManager::install(Plugin * plugin, DynLib * dynlib, std::string const & path)
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
PluginManager::uninstall(Plugin * plugin)
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
PluginManager::addEigenSolverFactory(char const * name, Algorithms::EigenSolverFactory * factory)
{
  Application::getEigenSolverManager().installFactory(name, factory);
}

void
PluginManager::removeEigenSolverFactory(char const * name)
{
  Application::getEigenSolverManager().uninstallFactory(name);
}

void
PluginManager::addLinearSolverFactory(char const * name, Algorithms::LinearSolverFactory * factory)
{
  Application::getLinearSolverManager().installFactory(name, factory);
}

void
PluginManager::removeLinearSolverFactory(char const * name)
{
  Application::getLinearSolverManager().uninstallFactory(name);
}

void
PluginManager::addRenderSystemFactory(char const * name, Graphics::RenderSystemFactory * factory)
{
  Application::getRenderSystemManager().installFactory(name, factory);
}

void
PluginManager::removeRenderSystemFactory(char const * name)
{
  Application::getRenderSystemManager().uninstallFactory(name);
}

void
PluginManager::addNumericalOptimizerFactory(char const * name, Algorithms::NumericalOptimizerFactory * factory)
{
  Application::getNumericalOptimizerManager().installFactory(name, factory);
}

void
PluginManager::removeNumericalOptimizerFactory(char const * name)
{
  Application::getNumericalOptimizerManager().uninstallFactory(name);
}

} // namespace Thea
