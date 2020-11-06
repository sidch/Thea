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

#ifndef __Thea_IPlugin_hpp__
#define __Thea_IPlugin_hpp__

#include "Common.hpp"
#include "List.hpp"
#include "Map.hpp"

namespace Thea {

// Forward declarations
namespace Algorithms {

class IEigenSolverFactory;
class ILinearSolverFactory;
class INumericalOptimizerFactory;

} // namespace Graphics

namespace Graphics {

class IRenderSystemFactory;

} // namespace Graphics

/** Interface for a plugin that may be statically or dynamically loaded. */
class THEA_API IPlugin
{
  public:
    /** Destructor. */
    virtual ~IPlugin() = 0;

    /** Get the name of the plugin. */
    virtual char const * THEA_ICALL getName() const = 0;

    /**
     * Installation routine for the plugin. This is called when the plugin is first registered with the plugin manager, and
     * executes initialization tasks that depend only on core functionality and not on the availability of other plugins.
     *
     * @see uninstall()
     */
    virtual void THEA_ICALL install() = 0;

    /**
     * Startup routine for the plugin. This function assumes all dependencies have been loaded and can hence interact with
     * other plugins. If called multiple times without an intervening shutdown(), the second and subsequent calls should have
     * no effect.
     *
     * @see startup()
     */
    virtual void THEA_ICALL startup() = 0;

    /**
     * Shutdown routine for the plugin. This function assumes all dependencies are still loaded and can hence interact with
     * other plugins. If called multiple times without an intervening startup(), the second and subsequent calls should have no
     * effect.
     *
     * @see shutdown()
     */
    virtual void THEA_ICALL shutdown() = 0;

    /**
     * Uninstallation routine for the plugin. This is called when the plugin is removed from the plugin manager, and executes
     * cleanup tasks that depend only on core functionality and not on the availability of other plugins.
     *
     * @see install()
     */
    virtual void THEA_ICALL uninstall() = 0;

}; // class IPlugin

inline IPlugin::~IPlugin() {}

// Forward declaration
class DynLib;

/** Interface provided to DLL's for registering rendersystem factories etc. */
class THEA_API IFactoryRegistry
{
  public:
    /** Destructor. */
    virtual ~IFactoryRegistry() = 0;

    /** Register an eigensolver factory. */
    virtual void THEA_ICALL addEigenSolverFactory(char const * name, Algorithms::IEigenSolverFactory * factory) = 0;

    /** Unregister an eigensolver factory. */
    virtual void THEA_ICALL removeEigenSolverFactory(char const * name) = 0;

    /** Register a linear solver factory. */
    virtual void THEA_ICALL addLinearSolverFactory(char const * name, Algorithms::ILinearSolverFactory * factory) = 0;

    /** Unregister a linear solver factory. */
    virtual void THEA_ICALL removeLinearSolverFactory(char const * name) = 0;

    /** Register a rendersystem factory. */
    virtual void THEA_ICALL addRenderSystemFactory(char const * name, Graphics::IRenderSystemFactory * factory) = 0;

    /** Unregister a rendersystem factory. */
    virtual void THEA_ICALL removeRenderSystemFactory(char const * name) = 0;

    /** Register a numerical optimizer factory. */
    virtual void THEA_ICALL addNumericalOptimizerFactory(char const * name,
                                                         Algorithms::INumericalOptimizerFactory * factory) = 0;

    /** Unregister a numerical optimizer factory. */
    virtual void THEA_ICALL removeNumericalOptimizerFactory(char const * name) = 0;

}; // class IFactoryRegistry

inline IFactoryRegistry::~IFactoryRegistry() {}

/** Manages set of installed plugins (static or dynamic). There should be exactly one active PluginManager per application. */
class THEA_API PluginManager : public virtual IFactoryRegistry
{
  public:
    /** Destructor. */
    ~PluginManager();

    /**
     * Load a dynamically linked plugin from a path. The plugin will be automatically installed via install().
     *
     * @return A pointer to the new or previously loaded plugin on success, nullptr on failure.
     *
     * @see unload()
     */
    IPlugin * load(std::string const & path);

    /**
     * Unload a dynamically linked plugin, using the same path specified to load() it. The plugin will be automatically
     * uninstalled via uninstall().
     *
     * \warning This does <b>not</b> call the plugin's \link IPlugin::shutdown() shutdown() \endlink routine, so be sure to call
     * it explicitly before calling this function.
     *
     * @see load()
     */
    void unload(std::string const & path);

    /**
     * Install a new plugin. The plugin is registered with the manager and its \link IPlugin::install() install() \endlink
     * function called.
     *
     * @note Do not use this function to install plugins from dynamic libraries: use load() instead.
     *
     * @param plugin A pointer to the plugin to be installed, which should not be null.
     *
     * @see uninstall()
     */
    void install(IPlugin * plugin);

    /**
     * Uninstall a plugin. The plugin is deregistered from the manager and its \link IPlugin::uninstall() uninstall() \endlink
     * function called.
     *
     * @note Do not use this function to uninstall plugins from dynamic libraries: use unload() instead.
     *
     * @param plugin A pointer to the plugin to be uninstalled.
     *
     * @see install()
     */
    void uninstall(IPlugin * plugin);

    /**
     * Execute the \link IPlugin::startup() startup() \endlink routine of all plugins, in the order in which they were installed.
     */
    void startupAllPlugins();

    /**
     * Execute the \link IPlugin::shutdown() shutdown() \endlink routine of all plugins, in the reverse of the order in which
     * they were installed.
     */
    void shutdownAllPlugins();

    /**
     * Unload/uninstall plugins. This will initially call shutdownAllPlugins(), so if you want to shutdown the plugins in an
     * order other than the reverse of the installation order, you should do it manually before calling this function.
     */
    void unloadAllPlugins();

    void addEigenSolverFactory(char const * name, Algorithms::IEigenSolverFactory * factory);
    void addLinearSolverFactory(char const * name, Algorithms::ILinearSolverFactory * factory);
    void addRenderSystemFactory(char const * name, Graphics::IRenderSystemFactory * factory);
    void addNumericalOptimizerFactory(char const * name, Algorithms::INumericalOptimizerFactory * factory);

    void removeEigenSolverFactory(char const * name);
    void removeLinearSolverFactory(char const * name);
    void removeRenderSystemFactory(char const * name);
    void removeNumericalOptimizerFactory(char const * name);

  private:
    /** Holds a plugin and its associated dynamic library. */
    struct PluginDynLib
    {
      IPlugin * plugin;
      DynLib * dynlib;

      PluginDynLib(IPlugin * p = nullptr, DynLib * d = nullptr) : plugin(p), dynlib(d) {}

    }; // struct PluginDynLib

    typedef Map<std::string, PluginDynLib>  DynLibMap;   ///< Maps paths to the corresponding dynamic libraries.
    typedef List<IPlugin *>                  PluginList;  ///< List of plugins.

    friend class DynLibManager;

    /**
     * Install a new plugin. The plugin is registered with the manager and its \link IPlugin::install() install() \endlink
     * function called.
     *
     * @param plugin A pointer to the plugin to be installed, which should not be null.
     * @param dynlib A pointer to the dynamic library from which the plugin was loaded, if any.
     * @param path The path to the dynamic library from which the plugin was loaded, if any.
     */
    void install(IPlugin * plugin, DynLib * dynlib, std::string const & path);

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
