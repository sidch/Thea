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

#ifndef __Thea_Plugin_hpp__
#define __Thea_Plugin_hpp__

#include "Common.hpp"
#include "List.hpp"
#include "Map.hpp"

namespace Thea {

// Forward declarations
namespace Algorithms {

class EigenSolverFactory;
class LinearSolverFactory;
class NumericalOptimizerFactory;

} // namespace Graphics

namespace Graphics {

class RenderSystemFactory;

} // namespace Graphics

/** The interface for a plugin that may be statically or dynamically loaded. */
class THEA_API Plugin
{
  public:
    /** Destructor. */
    virtual ~Plugin() {}

    /** Get the name of the plugin. */
    virtual char const * getName() const = 0;

    /**
     * Installation routine for the plugin. This is called when the plugin is first registered with the plugin manager, and
     * executes initialization tasks that depend only on core functionality and not on the availability of other plugins.
     *
     * @see uninstall()
     */
    virtual void install() = 0;

    /**
     * Startup routine for the plugin. This function assumes all dependencies have been loaded and can hence interact with
     * other plugins. If called multiple times without an intervening shutdown(), the second and subsequent calls should have
     * no effect.
     *
     * @see startup()
     */
    virtual void startup() = 0;

    /**
     * Shutdown routine for the plugin. This function assumes all dependencies are still loaded and can hence interact with
     * other plugins. If called multiple times without an intervening startup(), the second and subsequent calls should have no
     * effect.
     *
     * @see shutdown()
     */
    virtual void shutdown() = 0;

    /**
     * Uninstallation routine for the plugin. This is called when the plugin is removed from the plugin manager, and executes
     * cleanup tasks that depend only on core functionality and not on the availability of other plugins.
     *
     * @see install()
     */
    virtual void uninstall() = 0;

}; // class Plugin

// Forward declaration
class DynLib;

/** Interface provided to DLL's for registering rendersystem factories etc. */
class THEA_API FactoryRegistry
{
  public:
    /** Destructor. */
    virtual ~FactoryRegistry() {}

    /** Register an eigensolver factory. */
    virtual void addEigenSolverFactory(char const * name, Algorithms::EigenSolverFactory * factory) = 0;

    /** Unregister an eigensolver factory. */
    virtual void removeEigenSolverFactory(char const * name) = 0;

    /** Register a linear solver factory. */
    virtual void addLinearSolverFactory(char const * name, Algorithms::LinearSolverFactory * factory) = 0;

    /** Unregister a linear solver factory. */
    virtual void removeLinearSolverFactory(char const * name) = 0;

    /** Register a rendersystem factory. */
    virtual void addRenderSystemFactory(char const * name, Graphics::RenderSystemFactory * factory) = 0;

    /** Unregister a rendersystem factory. */
    virtual void removeRenderSystemFactory(char const * name) = 0;

    /** Register a numerical optimizer factory. */
    virtual void addNumericalOptimizerFactory(char const * name, Algorithms::NumericalOptimizerFactory * factory) = 0;

    /** Unregister a numerical optimizer factory. */
    virtual void removeNumericalOptimizerFactory(char const * name) = 0;

}; // class FactoryRegistry

/** Manages set of installed plugins (static or dynamic). There should be exactly one active PluginManager per application. */
class THEA_API PluginManager : public FactoryRegistry
{
  public:
    /** Destructor. */
    ~PluginManager();

    /**
     * Load a dynamically linked plugin from a path. The plugin will be automatically installed via install().
     *
     * @return A pointer to the new or previously loaded plugin on success, NULL on failure.
     *
     * @see unload()
     */
    Plugin * load(std::string const & path);

    /**
     * Unload a dynamically linked plugin, using the same path specified to load() it. The plugin will be automatically
     * uninstalled via uninstall().
     *
     * \warning This does <b>not</b> call the plugin's \link Plugin::shutdown() shutdown() \endlink routine, so be sure to call
     * it explicitly before calling this function.
     *
     * @see load()
     */
    void unload(std::string const & path);

    /**
     * Install a new plugin. The plugin is registered with the manager and its \link Plugin::install() install() \endlink
     * function called.
     *
     * @note Do not use this function to install plugins from dynamic libraries: use load() instead.
     *
     * @param plugin A pointer to the plugin to be installed, which should not be null.
     *
     * @see uninstall()
     */
    void install(Plugin * plugin);

    /**
     * Uninstall a plugin. The plugin is deregistered from the manager and its \link Plugin::uninstall() uninstall() \endlink
     * function called.
     *
     * @note Do not use this function to uninstall plugins from dynamic libraries: use unload() instead.
     *
     * @param plugin A pointer to the plugin to be uninstalled.
     *
     * @see install()
     */
    void uninstall(Plugin * plugin);

    /**
     * Execute the \link Plugin::startup() startup() \endlink routine of all plugins, in the order in which they were installed.
     */
    void startupAllPlugins();

    /**
     * Execute the \link Plugin::shutdown() shutdown() \endlink routine of all plugins, in the reverse of the order in which
     * they were installed.
     */
    void shutdownAllPlugins();

    /**
     * Unload/uninstall plugins. This will initially call shutdownAllPlugins(), so if you want to shutdown the plugins in an
     * order other than the reverse of the installation order, you should do it manually before calling this function.
     */
    void unloadAllPlugins();

    void addEigenSolverFactory(char const * name, Algorithms::EigenSolverFactory * factory);
    void addLinearSolverFactory(char const * name, Algorithms::LinearSolverFactory * factory);
    void addRenderSystemFactory(char const * name, Graphics::RenderSystemFactory * factory);
    void addNumericalOptimizerFactory(char const * name, Algorithms::NumericalOptimizerFactory * factory);

    void removeEigenSolverFactory(char const * name);
    void removeLinearSolverFactory(char const * name);
    void removeRenderSystemFactory(char const * name);
    void removeNumericalOptimizerFactory(char const * name);

  private:
    /** Holds a plugin and its associated dynamic library. */
    struct PluginDynLib
    {
      Plugin * plugin;
      DynLib * dynlib;

      PluginDynLib(Plugin * p = NULL, DynLib * d = NULL) : plugin(p), dynlib(d) {}

    }; // struct PluginDynLib

    typedef TheaMap<std::string, PluginDynLib>  DynLibMap;   ///< Maps paths to the corresponding dynamic libraries.
    typedef TheaList<Plugin *>                  PluginList;  ///< List of plugins.

    friend class DynLibManager;

    /**
     * Install a new plugin. The plugin is registered with the manager and its \link Plugin::install() install() \endlink
     * function called.
     *
     * @param plugin A pointer to the plugin to be installed, which should not be null.
     * @param dynlib A pointer to the dynamic library from which the plugin was loaded, if any.
     * @param path The path to the dynamic library from which the plugin was loaded, if any.
     */
    void install(Plugin * plugin, DynLib * dynlib, std::string const & path);

    /** Unload a plugin, given a handle to it in the dynlib map. */
    void unload(DynLibMap::iterator lib);

    /**
     * Unload all plugins that were loaded from dynamic libraries.
     *
     * @see unloadAllPlugins()
     */
    void unloadAllDylibs();

    DynLibMap   dynlibs;      ///< Set of dynamically loaded libraries.
    PluginList  plugins;      ///< Set of installed plugins.

}; // class PluginManager

} // namespace Thea

#endif
