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
#include "DynLib.hpp"
#include "Algorithms/LinearSolver.hpp"
#include "Algorithms/EigenSolver.hpp"
#include "Graphics/RenderSystem.hpp"

namespace Thea {

typedef Plugin * (*THEA_DLL_START_PLUGIN)(void);
typedef void     (*THEA_DLL_STOP_PLUGIN) (void);

// Static variables
bool                       PluginManager::initialized = false;
PluginManager::DynLibMap   PluginManager::dynlibs;
PluginManager::PluginList  PluginManager::plugins;

void
PluginManager::init()
{
  if (initialized) return;

  // Initialize all managers
  Algorithms::LinearSolverManager::_init();
  Algorithms::EigenSolverManager::_init();
  Graphics::RenderSystemManager::_init();

  initialized = true;
}

void
PluginManager::finish()
{
  if (!initialized) return;

  unloadAllPlugins();

  // Shutdown all managers
  Graphics::RenderSystemManager::_finish();
  Algorithms::EigenSolverManager::_finish();
  Algorithms::LinearSolverManager::_finish();

  initialized = false;
}

Plugin *
PluginManager::load(std::string const & path)
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  if (dynlibs.find(path) != dynlibs.end()) return NULL;

  DynLib * dynlib = DynLibManager::load(path);
  dynlibs[path] = dynlib;

  THEA_DLL_START_PLUGIN start_func = (THEA_DLL_START_PLUGIN)dynlib->getSymbol("dllStartPlugin");
  return start_func();  // this must call install(), which adds the plugin to the installed list
}

void
PluginManager::unload(std::string const & path)
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  DynLibMap::iterator dyn_loaded = dynlibs.find(path);
  if (dyn_loaded != dynlibs.end())
  {
    THEA_DLL_STOP_PLUGIN stop_func = (THEA_DLL_STOP_PLUGIN)dyn_loaded->second->getSymbol("dllStopPlugin");
    stop_func();  // this must call uninstall(), which removes the plugin from the installed list

    DynLibManager::unload(dyn_loaded->second);
    dynlibs.erase(dyn_loaded);
  }
}

void
PluginManager::startupAllPlugins()
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
    (*pi)->startup();
}

void
PluginManager::shutdownAllPlugins()
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  for (PluginList::reverse_iterator pi = plugins.rbegin(); pi != plugins.rend(); ++pi)
    (*pi)->shutdown();
}

void
PluginManager::unloadAllPlugins()
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  shutdownAllPlugins();

  for (DynLibMap::iterator di = dynlibs.begin(); di != dynlibs.end(); ++di)
  {
    THEA_DLL_STOP_PLUGIN stop_func = (THEA_DLL_STOP_PLUGIN)di->second->getSymbol("dllStopPlugin");
    stop_func();  // this must call uninstall(), which removes the plugin from the installed list

    DynLibManager::unload(di->second);
  }

  dynlibs.clear();

  // There should be only static libs left now
  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
  {
    THEA_LOG << "PluginManager: Uninstalling plugin '" << (*pi)->getName() << '\'';
    (*pi)->uninstall();
  }

  plugins.clear();
}

void
PluginManager::install(Plugin * plugin)
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");
  alwaysAssertM(plugin, "PluginManager: Plugin pointer cannot be null");

  THEA_LOG << "PluginManager: Installing plugin '" << plugin->getName() << '\'';
  plugins.push_back(plugin);
  plugin->install();
}

void
PluginManager::uninstall(Plugin * plugin)
{
  alwaysAssertM(initialized, "PluginManager: Manager not initialized");

  if (!plugin) return;

  THEA_LOG << "PluginManager: Uninstalling plugin '" << plugin->getName() << '\'';
  plugin->uninstall();

  for (PluginList::iterator pi = plugins.begin(); pi != plugins.end(); ++pi)
    if (*pi == plugin)
    {
      plugins.erase(pi);
      break;
    }
}

} // namespace Thea
