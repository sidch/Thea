//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_App_hpp__
#define __Browse3D_App_hpp__

#include "Common.hpp"
#include <QObject>

namespace Thea {

class Plugin;

namespace Graphics {

class RenderSystem;
class RenderSystemFactory;

} // namespace Graphics

} // namespace Thea

namespace Browse3D {

class MainWindow;

/** Main application properties and functions. */
class App : public QObject
{
    Q_OBJECT

  public:
    /** Run-time options, typically specified on the command-line or via a configuration file. */
    struct Options
    {
      /** Constructor. */
      Options();

      QString plugin_dir;           ///< Directory containing plugins.
      QString resource_dir;         ///< Directory containing resources (shaders, textures, ...)
      QString working_dir;          ///< The application's initial working directory.
      QString model;                ///< The initial shape to load.

      bool bg_plain;                ///< Draw the background in a single plain color?
      ColorRGB bg_color;            ///< Background color.

      bool accentuate_features;     ///< Make feature distributions easier to view?
      bool show_graph;              ///< Show point adjacency graph, if available?

    }; // struct Options

    /** Constructor. */
    App();

    /**
     * Initialize the application from command-line parameters. Should be called on program startup.
     *
     * @return True if the application was successfully initialized, else false.
     */
    bool init(int argc, char * argv[]);  // return true if the app should continue, false if it should quit

    /** Get the set of program options. */
    Options const & options() const { return opts; }

    /** Get a textual representation of all the program options. */
    std::string optsToString() const;

    /** Check if the application has a properly constructed rendersystem */
    bool hasRenderSystem() const { return has_render_system.value(); }

    /** Get a handle to the application's main window. */
    MainWindow * getMainWindow() { return main_window; }

    /** Get a handle to the application's rendersystem. */
    Graphics::RenderSystem * getRenderSystem() { return render_system; }

  public slots:
    /** Cleanup data held by the application. Call right before program exit. */
    void cleanup();

  private:
    /**
     * Parse program options from the command-line, or from a configuration file specified on the command-line.
     *
     * @return True on success, false on failure.
     */
    bool parseOptions(int argc, char * argv[]);

    /** Create the application's main window. */
    void createMainWindow();

    /** Load necessary or optional plugins. */
    void loadPlugins();

    /** Create the application's rendersystem. */
    void createRenderSystem();

    Options opts;
    MainWindow * main_window;
    AtomicInt32 has_render_system;
    Plugin * gl_plugin;
    Graphics::RenderSystemFactory * render_system_factory;
    Graphics::RenderSystem * render_system;

}; // class App

/** Get the global application object, constructed on program startup. */
inline App & app()
{
  static App a;
  return a;
}

} // namespace Browse3D

#endif
