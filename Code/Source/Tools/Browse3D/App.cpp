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
#include <wx/cmdline.h>
#include <wx/image.h>
#include <boost/program_options.hpp>
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
  has_render_system(0),
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
      << "\n  fancy-points = " << opts.fancy_points
      << "\n  fancy-colors = " << opts.fancy_colors
      << "\n  point-scale = " << opts.point_scale
      << "\n  no-axes = " << opts.no_axes
      << "\n  no-shading= " << opts.no_shading
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
  {
    transform = AffineTransform3::identity();
  }

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
  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i)  // omit the program path
    args.push_back(argv[i]);

  return parseOptions(args);
}

bool
App::parseOptions(std::vector<std::string> const & args)
{
  namespace po = boost::program_options;

  static std::string const usage("Usage: Browse3D [options] [model]");

  std::string conf_file;

  // po::options_description hidden;
  // hidden.add_options()("hidden-option", po::value<std::string>(&hidden_option), "");

  std::string app_dir(FilePath::parent(Application::programPath()));
  std::string def_resource_dir = FilePath::concat(app_dir, "../../../Resources");
#ifdef _MSC_VER
  // Visual Studio puts executables in Debug|Release subdirectory, so handle this possibility
  if (!FileSystem::directoryExists(def_resource_dir))
    def_resource_dir = FilePath::concat(app_dir, "../../../../Resources");
#endif

  std::string s_model;
  std::vector<std::string> s_overlays;
  std::string s_color;
  std::string s_bg_color;
  std::string s_material;

  po::options_description visible("Allowed options");
  visible.add_options()
          ("help,h",               "Print this help message")
          ("version,v",            "Print the program version")
          ("conf",                 po::value<std::string>(&conf_file)->default_value("Browse3D.conf"), "Configuration file (overridden by duplicate cmdline options)")
          ("plugin-dir",           po::value<std::string>(&opts.plugin_dir), "Plugins directory")
          ("resource-dir",         po::value<std::string>(&opts.resource_dir)->default_value(def_resource_dir), "Resources directory")
          ("working-dir",          po::value<std::string>(&opts.working_dir)->default_value("."), "Working directory")
          ("model",                po::value<std::string>(&s_model), "Model to load on startup, with optional transform")
          ("overlay",              po::value< std::vector<std::string> >(&s_overlays), "Overlay model(s) to load on startup")
          ("features,f",           po::value<std::string>(&opts.features), "Directory/file containing features to load")
          ("elem-labels,l",        po::value<std::string>(&opts.elem_labels), "Directory/file containing face/point labels to load")
          ("emph-features,e",      "Make feature distributions easier to view")
          ("color-cube,3",         "Map 0-centered 3D feature sets to RGB color-cube, if --emph-features")
          ("normals,n",            "Draw normals as arrows")
          ("graph,g",              "Show point adjacency graph")
          ("color,c",              po::value<std::string>(&s_color), "Model color")
          ("bg-color,b",           po::value<std::string>(&s_bg_color), "Background color")
          ("two-sided",            po::value<bool>(&opts.two_sided)->default_value(true), "Use two-sided lighting?")
          ("flat,0",               "Flat shade meshes")
          ("edges,j",              "Show mesh edges")
          ("material,k",           po::value<std::string>(&s_material), "Surface material coefficients (ka, kd, ks, ksp)")
          ("fancy-points",         "Draw points as shaded spheres")
          ("fancy-colors",         "Color points by a function of position")
          ("point-scale,p",        po::value<Real>(&opts.point_scale)->default_value(1), "Scale point sizes by this factor")
          ("no-axes",              "Hide the coordinate axes")
          ("no-shading",           "No shading, just render raw colors")
  ;

  po::options_description desc;
  desc.add(visible) /* .add(hidden) */ ;

  if (argc < 1)
  {
    THEA_CONSOLE << usage;
    std::cerr << visible;  // should be intercepted by out
    return 0;
  }

  po::positional_options_description pdesc;
  pdesc.add("model", 1);
  pdesc.add("overlay", -1);

  // Read cmdline options first (overrides conflicting config file values)
  po::parsed_options cmdline_parsed = po::basic_command_line_parser<char>(args).options(desc).positional(pdesc).run();
  po::variables_map vm;
  po::store(cmdline_parsed, vm);

  // Now read the config file, if it is found
  if (vm.count("conf") > 0 && FileSystem::fileExists(conf_file))
  {
    THEA_CONSOLE << "Reading options from config file:" << conf_file;

    std::ifstream conf_in(conf_file.c_str());
    po::parsed_options conf_file_parsed = po::parse_config_file(conf_in, desc);
    po::store(conf_file_parsed, vm);
  }

  po::notify(vm);

  bool quit = false;

  if (vm.count("version") > 0)
  {
    THEA_CONSOLE << "Browse3D version 2.0";
    THEA_CONSOLE << "Siddhartha Chaudhuri, 2016";
    quit = true;
  }

  if (vm.count("help") > 0)
  {
    if (quit) THEA_CONSOLE << "";
    THEA_CONSOLE << usage;
    THEA_CONSOLE << visible;  // should be intercepted by THEA_CONSOLE
    quit = true;
  }

  if (quit)
    return false;

  opts.plugin_dir    =  FileSystem::resolve(opts.plugin_dir);
  opts.resource_dir  =  FileSystem::resolve(opts.resource_dir);
  opts.working_dir   =  FileSystem::resolve(opts.working_dir);

  if (!s_model.empty())
  {
    if (!parseModel(s_model, opts.model, opts.model_transform))
      return false;
  }

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
  opts.accentuate_features  =  (vm.count("emph-features") > 0);
  opts.color_cube_features  =  (vm.count("color-cube") > 0);
  opts.show_normals         =  (vm.count("normals") > 0);
  opts.show_graph           =  (vm.count("graph") > 0);
  opts.flat                 =  (vm.count("flat") > 0);
  opts.edges                =  (vm.count("edges") > 0);
  opts.fancy_points         =  (vm.count("fancy-points") > 0);
  opts.fancy_colors         =  (vm.count("fancy-colors") > 0);
  opts.no_axes              =  (vm.count("no-axes") > 0);
  opts.no_shading           =  (vm.count("no-shading") > 0);

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
    Vector4 m;
    int num_fields = parseVector(s_material, m);
    if (num_fields <= 0)
    {
      THEA_ERROR << "Could not parse material: " << s_material;
      return false;
    }

    for (int i = 0; i < num_fields; ++i)
      if (m[i] >= -0.001)
        opts.material[i] = m[i];
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

  // Start up the plugin (a GL context should already exist in a QGLWidget)
  gl_plugin->startup();
}

void
App::createRenderSystem()
{
  render_system_factory = Application::getRenderSystemManager().getFactory("OpenGL");
  render_system = render_system_factory->createRenderSystem("OpenGL");
  has_render_system = 1;
}

//=============================================================================================================================
// GUI callbacks etc
//=============================================================================================================================

bool
App::OnInit()
{
  if (!parseOptions(this->argc, this->argv))
    return false;

  THEA_CONSOLE << "Started Browse3D\n";
  THEA_CONSOLE << optsToString();

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
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

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
