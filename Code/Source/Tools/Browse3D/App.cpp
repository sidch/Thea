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

#include "App.hpp"
#include "MainWindow.hpp"
#include "../../Application.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../IPlugin.hpp"
#include "../../Graphics/IRenderSystem.hpp"
#include "../../ThirdParty/CLI11/CLI11.hpp"
#include <wx/cmdline.h>
#include <wx/image.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>

IMPLEMENT_APP(Browse3D::App);

namespace Browse3D {

namespace AppInternal {

static ColorRgba  const  DEFAULT_COLOR(1.0f, 0.9f, 0.8f, 1.0f);
static Vector4    const  DEFAULT_MATERIAL(0.3f, 0.7f, 0.2f, 25);

} // namespace AppInternal

App::Options::Options()
: accentuate_features(false), color_cube_features(false), show_normals(false), show_graph(false),
  color(AppInternal::DEFAULT_COLOR), bg_plain(false), bg_color(ColorRgb::black()), two_sided(true), flat(false), edges(false),
  material(AppInternal::DEFAULT_MATERIAL), fancy_points(false), fancy_colors(false), point_scale(1), no_axes(false)
{
}

App::App()
: main_window(nullptr),
  has_render_system(false),
  gl_plugin(nullptr),
  render_system_factory(nullptr),
  render_system(nullptr)
{
}

std::string
App::optsToString() const
{
  std::ostringstream oss;

  oss << "\nProgram options"
      << "\n=============== "
      << "\n  plugin-dir = " << opts.plugin_dir
      << "\n  resource-dir = " << opts.resource_dir
      << "\n  working-dir = " << opts.working_dir
      << "\n  model = " << opts.model
      << "\n  overlays = { " << stringJoin(opts.overlays, ", ") << " }"
      << "\n  features = " << opts.features
      << "\n  elem-labels = " << opts.elem_labels
      << "\n  emph-features = " << opts.accentuate_features
      << "\n  color-cube = " << opts.color_cube_features
      << "\n  normals = " << opts.show_normals
      << "\n  graph = " << opts.show_graph
      << "\n  color = " << opts.color.toString()
      << "\n  bg-plain = " << opts.bg_plain
      << "\n  bg-color = " << opts.bg_color.toString()
      << "\n  two-sided = " << opts.two_sided
      << "\n  flat = " << opts.flat
      << "\n  edges = " << opts.edges
      << "\n  material = " << toString(opts.material)
      << "\n  matcap = " << opts.matcap
      << "\n  fancy-points = " << opts.fancy_points
      << "\n  fancy-colors = " << opts.fancy_colors
      << "\n  point-scale = " << opts.point_scale
      << "\n  no-axes = " << opts.no_axes
      << "\n  no-shading = " << opts.no_shading
      << '\n';

  return oss.str();
}

bool
parseModel(std::string const & str, std::string & path, AffineTransform3 & transform)
{
  Array<std::string> fields;
  stringSplit(str, '|', fields);
  if (fields.empty() || fields.size() > 2)
  {
    THEA_ERROR << "Could not parse model: " << str;
    return false;
  }

  path = trimWhitespace(fields[0]);

  if (fields.size() >= 2)
  {
    float m[3][4];
    if (std::sscanf(fields[1].c_str(), " [ L : [ %f , %f , %f ; %f , %f , %f ; %f , %f , %f ] , T : ( %f , %f , %f ) ]",
                    &m[0][0], &m[0][1], &m[0][2],
                    &m[1][0], &m[1][1], &m[1][2],
                    &m[2][0], &m[2][1], &m[2][2],
                    &m[0][3], &m[1][3], &m[2][3]) == 12)
    {}
    else if (std::sscanf(fields[1].c_str(), " %f %f %f %f  %f %f %f %f  %f %f %f %f",
                         &m[0][0], &m[0][1], &m[0][2], &m[0][3],
                         &m[1][0], &m[1][1], &m[1][2], &m[1][3],
                         &m[2][0], &m[2][1], &m[2][2], &m[2][3]) == 12)
    {}
    else if (std::sscanf(fields[1].c_str(), " %f %f %f", &m[0][3], &m[1][3], &m[2][3]) == 3)
    {
      m[0][0] = m[1][1] = m[2][2] = 1;
      m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] = m[2][1] = 0;
    }
    else
    {
      THEA_ERROR << "Could not parse model transformation: " << fields[1];
      return false;
    }

    transform = AffineTransform3((Matrix3() << m[0][0], m[0][1], m[0][2],
                                               m[1][0], m[1][1], m[1][2],
                                               m[2][0], m[2][1], m[2][2]).finished(),
                                 Vector3(m[0][3], m[1][3], m[2][3]));
  }
  else
    transform.setIdentity();

  return true;
}

// Returns number of components found
int
parseVector(std::string const & str, Vector4 & v)
{
  std::string s = trimWhitespace(str);
  if (s.length() >= 2 && s[0] == '(' && s[s.length() - 1] == ')')
    s = s.substr(1, s.length() - 2);

  double x[4];
  size_t num_parsed = std::sscanf(s.c_str(), " %lf , %lf , %lf , %lf", &x[0], &x[1], &x[2], &x[3]);
  if (num_parsed < 2)
    num_parsed = std::sscanf(s.c_str(), " %lf %lf %lf %lf", &x[0], &x[1], &x[2], &x[3]);

  v.fill(-1);
  for (size_t i = 0; i < num_parsed; ++i)
    v[i] = (Real)x[i];

  return (int)num_parsed;
}

bool
App::parseOptions(int argc, char * argv[])
{
  Array<std::string> args;
  for (int i = 1; i < argc; ++i)  // omit the program path
    args.push_back(argv[i]);

  return parseOptions(args);
}

bool
App::parseOptions(Array<std::string> const & args)
{
  std::string app_dir(FilePath::parent(Application::programPath()));
  std::string def_resource_dir = FilePath::concat(app_dir, "../../../Resources");
#ifdef _MSC_VER
  // Visual Studio puts executables in Debug|Release subdirectory, so handle this possibility
  if (!FileSystem::directoryExists(def_resource_dir))
    def_resource_dir = FilePath::concat(app_dir, "../../../../Resources");
#endif

  CLI::App app{"Browse3D"};

  std::string s_model;
  Array<std::string> s_overlays;
  std::string s_color;
  std::string s_bg_color;
  std::string s_material;

  app.set_config("--conf", "Browse3D.conf", "Configuration file (overridden by duplicate cmdline options)");
  app.add_option("--plugin-dir", opts.plugin_dir, "Plugins directory");
  app.add_option("--resource-dir", opts.resource_dir, "Resources directory")->default_val(def_resource_dir);
  app.add_option("--working-dir", opts.working_dir, "Working directory")->default_val(".");
  app.add_option("model", s_model, "Model to load on startup, with optional transform");
  app.add_option("overlay", s_overlays, "Overlay model(s) to load on startup");
  app.add_option("-f,--features", opts.features, "Directory/file containing features to load");
  app.add_option("-l,--labels", opts.elem_labels, "Directory/file containing face/point labels to load");
  app.add_flag("-e,--emph-features", "Make feature distributions easier to view");
  app.add_flag("-3,--color-cube", "Map 0-centered 3D feature sets to RGB color-cube, if --emph-features");
  app.add_flag("-n,--normals", "Draw normals as arrows");
  app.add_flag("-g,--graph", "Show point adjacency graph");
  app.add_option("-c,--color", s_color, "Model color");
  app.add_option("-b,--bg-color", s_bg_color, "Background color");
  app.add_option("--two-sided", opts.two_sided, "Use two-sided lighting?")->default_val(true);
  app.add_flag("-0,--flat", "Flat shade meshes");
  app.add_flag("-j,--edges", "Show mesh edges");
  app.add_option("-k,--material", s_material, "Surface material coefficients (ka, kd, ks, ksp), or path to a matcap image");
  app.add_flag("--fancy-points", "Draw points as shaded spheres");
  app.add_flag("--fancy-colors", "Color points by a function of position");
  app.add_option("-p,--point-scale", opts.point_scale, "Scale point sizes by this factor")->default_val(1);
  app.add_flag("--no-axes", "Hide the coordinate axes");
  app.add_flag("--no-shading", "No shading, just render raw colors");

  try
  {
    app.parse(argc, argv);
  }
  catch (CLI::ParseError const & e)
  {
    app.exit(e);
    return false;
  }

  opts.plugin_dir    =  FileSystem::resolve(opts.plugin_dir);
  opts.resource_dir  =  FileSystem::resolve(opts.resource_dir);
  opts.working_dir   =  FileSystem::resolve(opts.working_dir);

  if (!s_model.empty() && !parseModel(s_model, opts.model, opts.model_transform))
    return false;

  opts.overlays.clear();
  opts.overlay_transforms.clear();
  for (size_t i = 0; i < s_overlays.size(); ++i)
  {
    std::string path;
    AffineTransform3 transform;
    if (!parseModel(s_overlays[i], path, transform))
      return false;

    opts.overlays.push_back(path);
    opts.overlay_transforms.push_back(transform);
  }

  opts.features             =  FileSystem::resolve(opts.features);
  opts.elem_labels          =  FileSystem::resolve(opts.elem_labels);

  opts.accentuate_features  =  (app.count("--emph-features") > 0);
  opts.color_cube_features  =  (app.count("--color-cube") > 0);
  opts.show_normals         =  (app.count("--normals") > 0);
  opts.show_graph           =  (app.count("--graph") > 0);
  opts.flat                 =  (app.count("--flat") > 0);
  opts.edges                =  (app.count("--edges") > 0);
  opts.fancy_points         =  (app.count("--fancy-points") > 0);
  opts.fancy_colors         =  (app.count("--fancy-colors") > 0);
  opts.no_axes              =  (app.count("--no-axes") > 0);
  opts.no_shading           =  (app.count("--no-shading") > 0);

  if (!s_color.empty())
  {
    std::stringstream ss; ss << std::hex << s_color;
    uint32 argb; ss >> argb;
    opts.color = ColorRgb::fromARGB(argb);  // Rgb and not Rgba to ensure it's opaque
  }

  if (!s_bg_color.empty())
  {
    std::stringstream ss; ss << std::hex << s_bg_color;
    uint32 argb; ss >> argb;
    opts.bg_plain = true;
    opts.bg_color = ColorRgb::fromARGB(argb);  // Rgb and not Rgba to ensure it's opaque
  }
  else
  {
    opts.bg_plain = false;
    opts.bg_color = ColorRgb::black();
  }

  if (!s_material.empty())
  {
    if (FileSystem::fileExists(s_material))
      opts.matcap = s_material;
    else
    {
      Vector4 m;
      int num_fields = parseVector(s_material, m);
      if (num_fields > 0)
      {
        for (int i = 0; i < num_fields; ++i)
          if (m[i] >= -0.001)
            opts.material[i] = m[i];
      }
      else
      {
        THEA_ERROR << "Invalid material: " << s_material;
        return false;
      }
    }
  }

  Application::setResourceArchive(opts.resource_dir);

  return true;
}

void
App::createMainWindow()
{
  // Create the main window, and hence a rendering context
  main_window = new MainWindow;
  main_window->Show(true);
}

void
App::loadPlugins()
{
  // Try to load the OpenGL plugin
  Array<std::string> plugin_dirs; plugin_dirs.push_back(opts.plugin_dir);
  std::string plugin_path = Application::getPluginPath("TheaPluginGL", &plugin_dirs);
  if (plugin_path.empty())
    throw Error("Could not locate OpenGL plugin 'TheaPluginGL'");

  THEA_CONSOLE << "Loading OpenGL plugin: " << plugin_path;
  gl_plugin = Application::getPluginManager().load(plugin_path);

  // Start up the plugin (a GL context should already exist in a UI widget)
  gl_plugin->startup();
}

void
App::createRenderSystem()
{
  if (!has_render_system.exchange(true))
  {
    render_system_factory = Application::getRenderSystemManager().getFactory("OpenGL");
    render_system = render_system_factory->createRenderSystem("OpenGL");

    THEA_CONSOLE << "\nRenderSystem: " << render_system->describeSystem();
  }
}

//=============================================================================================================================
// GUI callbacks etc
//=============================================================================================================================

bool
App::OnInit()
{
  if (!parseOptions(this->argc, this->argv))
    return false;

  THEA_CONSOLE << "Started Browse3D\n" << optsToString();

  wxImage::AddHandler(new wxPNGHandler);

  createMainWindow();

  // Load plugins and create a rendersystem
  loadPlugins();
  createRenderSystem();

  return true;
}

bool
App::OnExceptionInMainLoop()
{
  try
  {
    throw; // Rethrow the current exception.
  }
  THEA_CATCH(return -1;, ERROR, "%s", "An error occurred")

  // Exit the main loop and thus terminate the program.
  return false;
}

int
App::OnExit()
{
  if (render_system_factory)
    render_system_factory->destroyRenderSystem(render_system);

  Application::getPluginManager().unloadAllPlugins();

  return wxApp::OnExit();
}

} // namespace Browse3D
