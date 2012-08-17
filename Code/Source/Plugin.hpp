//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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
#include "NamedObject.hpp"

namespace Thea {

/** The interface for a plugin that may be statically or dynamically loaded. */
class THEA_API Plugin : public virtual NamedObject
{
  public:
    /** Destructor. */
    virtual ~Plugin() {}

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

/** Manages set of installed plugins (static or dynamic). */
class THEA_API PluginManager
{
  public:
    /**
     * Initialization routine. You must call this at the beginning of the program if you want to use any plugins. It's safe to
     * call this function multiple times -- a second or later call will be ignored unless it was immediately preceded by
     * finish().
     *
     * Once you call init(), you must call finish() at the end of your program to properly clean up.
     *
     * @see finish()
     */
    static void init();

    /**
     * Termination routine. You must call this at the end of the program after using any plugins. It's safe to call this
     * function multiple times -- a second or later call will be ignored unless it was immediately preceded by init(). It's also
     * safe to call this function without ever having called init() first.
     *
     * @see init()
     */
    static void finish();

    /**
     * Load a dynamically linked plugin from a path. The plugin's \link Plugin::install() install() \endlink function will be
     * automatically called.
     *
     * @return A pointer to the new plugin on success, NULL on failure or if the plugin has already been loaded.
     *
     * @see unload()
     */
    static Plugin * load(std::string const & path);

    /**
     * Unload a dynamically linked plugin, using the same path specified to load() it. The plugin's \link Plugin::uninstall()
     * uninstall() \endlink function will be automatically called.
     *
     * \warning This does <b>not</b> call the plugin's \link Plugin::shutdown() shutdown() \endlink routine, so be sure to call
     * it explicitly before calling this function.
     *
     * @see load()
     */
    static void unload(std::string const & path);

    /**
     * Execute the \link Plugin::startup() startup() \endlink routine of all plugins, in the order in which they were installed.
     */
    static void startupAllPlugins();

    /**
     * Execute the \link Plugin::shutdown() shutdown() \endlink routine of all plugins, in the reverse of the order in which
     * they were installed.
     */
    static void shutdownAllPlugins();

    /**
     * Unload/uninstall plugins. This will initially call shutdownAllPlugins(), so if you want to shutdown the plugins in an
     * order other than the reverse of the installation order, you should do it manually before calling this function.
     */
    static void unloadAllPlugins();

    /**
     * This function should be called by <code>dllStartPlugin()</code> to register the plugin. It may also be used to install
     * statically linked plugins.
     *
     * @param plugin A pointer to the plugin to be installed, which should not be null.
     *
     * @see uninstall()
     */
    static void install(Plugin * plugin);

    /**
     * This function should be called by <code>dllStopPlugin()</code> to register the plugin. It may also be used to uninstall
     * statically linked plugins.
     *
     * @param plugin A pointer to the plugin to be uninstalled.
     *
     * @see install()
     */
    static void uninstall(Plugin * plugin);

  private:
    typedef TheaMap<std::string, DynLib *>  DynLibMap;   ///< Maps paths to the corresponding dynamic libraries.
    typedef TheaList<Plugin *>              PluginList;  ///< List of plugins.

    static bool        initialized;  ///< Has the manager been initialized?
    static DynLibMap   dynlibs;      ///< Set of dynamically loaded libraries.
    static PluginList  plugins;      ///< Set of installed plugins.

}; // class PluginManager

} // namespace Thea

#endif
