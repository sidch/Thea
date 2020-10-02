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

#ifndef __Browse3D_App_hpp__
#define __Browse3D_App_hpp__

#include "Common.hpp"
#include "../../AffineTransform3.hpp"
#include "../../MatVec.hpp"
#include <wx/app.h>
#include <vector>

namespace Thea {

class IPlugin;

namespace Graphics {

class IRenderSystem;
class IRenderSystemFactory;

} // namespace Graphics

} // namespace Thea

namespace Browse3D {

class MainWindow;

/** Main application properties and functions. */
class App : public wxApp
{
  public:
    /** Run-time options, typically specified on the command-line or via a configuration file. */
    struct Options
    {
      /** Constructor. */
      Options();

      std::string plugin_dir;                              ///< Directory containing plugins.
      std::string resource_dir;                            ///< Directory containing resources (shaders, textures, ...)
      std::string working_dir;                             ///< The application's initial working directory.
      std::string model;                                   ///< The initial shape to load.
      AffineTransform3 model_transform;                    ///< The transformation of the initial shape.
      Array<std::string> overlays;                         ///< The initial overlays to load.
      Array<AffineTransform3> overlay_transforms;          ///< The transforms of the overlays.
      std::string features;                                ///< Path to directory/file containing features to load.
      std::string elem_labels;                             ///< Path to directory/file containing face/point labels to load.

      bool accentuate_features;                             ///< Make feature distributions easier to view?
      bool color_cube_features;                             /**< Map 0-centered 3D feature sets to RGB color-cube, if
                                                                 accentuate_features == true; */
      bool show_normals;                                    ///< Draw arrows for normals?
      bool show_graph;                                      ///< Show point adjacency graph, if available?

      ColorRgba color;                                      ///< Model color.
      bool bg_plain;                                        ///< Draw the background in a single plain color?
      ColorRgba bg_color;                                   ///< Background color.
      bool two_sided;                                       ///< Use two-sided lighting?
      bool flat;                                            ///< Flat-shade all meshes?
      Vector4 material;                                     ///< Surface material coefficients (kd, ka, ks, ksp).
      bool fancy_points;                                    ///< Draw points as shaded spheres?
      bool fancy_colors;                                    ///< Color points by a function of position?
      Real point_scale;                                     ///< Scale point sizes by this factor.
      bool no_axes;                                         ///< Hide the coordinate axes?
      bool no_shading;                                      ///< No shading, just render raw colors?

    }; // struct Options

    /** Constructor. */
    App();

    /** Get the set of program options. */
    Options const & options() const { return opts; }

    /** Get a textual representation of all the program options. */
    std::string optsToString() const;

    /** Check if the application has a properly constructed rendersystem */
    bool hasRenderSystem() const { return has_render_system.value(); }

    /** Get a handle to the application's main window. */
    MainWindow * getMainWindow() { return main_window; }

    /** Get a handle to the application's rendersystem. */
    Graphics::IRenderSystem * getRenderSystem() { return render_system; }

    //=========================================================================================================================
    // GUI callbacks etc
    //=========================================================================================================================

    /** Called when application is launched. */
    bool OnInit();

    /** Called when the main program loop throws an exception. */
    bool OnExceptionInMainLoop();

    /** Cleanup data held by the application. */
    int OnExit();

  private:
    /**
     * Parse program options from the command-line, or from a configuration file specified on the command-line.
     *
     * @return True on success, false on failure.
     */
    bool parseOptions(int argc, char * argv[]);

    /**
     * Parse program options from the command-line, or from a configuration file specified on the command-line.
     *
     * @return True on success, false on failure.
     */
    bool parseOptions(std::vector<std::string> const & args);

    /** Create the application's main window. */
    void createMainWindow();

    /** Load necessary or optional plugins. */
    void loadPlugins();

    /** Create the application's rendersystem. */
    void createRenderSystem();

    Options opts;
    MainWindow * main_window;
    AtomicInt32 has_render_system;
    IPlugin * gl_plugin;
    Graphics::IRenderSystemFactory * render_system_factory;
    Graphics::IRenderSystem * render_system;

}; // class App

} // namespace Browse3D

DECLARE_APP(Browse3D::App);

namespace Browse3D {

/** Get the global application object, constructed on program startup. */
inline App & app()
{
  return wxGetApp();
}

} // namespace Browse3D

#endif
