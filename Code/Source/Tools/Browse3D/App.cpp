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

#include "App.hpp"
#include "MainWindow.hpp"
#include "../../Application.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Plugin.hpp"
#include "../../Graphics/RenderSystem.hpp"
#include <wx/cmdline.h>
#include <wx/image.h>
#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>

IMPLEMENT_APP(Browse3D::App);

namespace Browse3D {

App::Options::Options()
: accentuate_features(false), color_cube_features(false), show_normals(false), show_graph(false), bg_plain(false),
  bg_color(ColorRGB::black()), two_sided(true), flat(false), fancy_points(false), fancy_colors(false)
{
}

App::App()
: main_window(NULL),
  has_render_system(0),
  gl_plugin(NULL),
  render_system_factory(NULL),
  render_system(NULL)
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
      << "\n  face-labels = " << opts.face_labels
      << "\n  emph-features = " << opts.accentuate_features
      << "\n  color-cube = " << opts.color_cube_features
      << "\n  show-normals = " << opts.show_normals
      << "\n  show-graph = " << opts.show_graph
      << "\n  bg-plain = " << opts.bg_plain
      << "\n  bg-color = " << opts.bg_color.toString()
      << "\n  two-sided = " << opts.two_sided
      << "\n  flat = " << opts.flat
      << "\n  fancy-points = " << opts.fancy_points
      << "\n  fancy-colors = " << opts.fancy_colors
      << "\n  point-scale = " << opts.point_scale
      << '\n';

  return oss.str();
}

bool
parseModel(std::string const & str, std::string & path, AffineTransform3 & transform)
{
  TheaArray<std::string> fields;
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

    transform = AffineTransform3(Matrix3(m[0][0], m[0][1], m[0][2],
                                         m[1][0], m[1][1], m[1][2],
                                         m[2][0], m[2][1], m[2][2]),
                                 Vector3(m[0][3], m[1][3], m[2][3]));
  }
  else
  {
    transform = AffineTransform3::identity();
  }

  return true;
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
#ifdef _MSC_VER
  // Visual Studio puts executables in Debug|Release subdirectory
  std::string def_resource_dir = FilePath::concat(app_dir, "../../../../Resources");
#else
  std::string def_resource_dir = FilePath::concat(app_dir, "../../../Resources");
#endif

  std::string s_model;
  std::vector<std::string> s_overlays;
  std::string s_bg_color;

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
          ("face-labels,l",        po::value<std::string>(&opts.face_labels), "Directory/file containing face labels to load")
          ("emph-features,e",      "Make feature distributions easier to view")
          ("color-cube,3",         "Map 0-centered 3D feature sets to RGB color-cube, if --emph-features")
          ("normals,n",            "Draw normals as arrows")
          ("graph,g",              "Show point adjacency graph")
          ("bg",                   po::value<std::string>(&s_bg_color), "Background color")
          ("two-sided",            po::value<bool>(&opts.two_sided)->default_value(true), "Use two-sided lighting?")
          ("flat",                 "Flat shade all meshes?")
          ("fancy-points",         "Draw points as shaded spheres?")
          ("fancy-colors,c",       "Color points by a function of position?")
          ("point-scale,s",        po::value<Real>(&opts.point_scale)->default_value(1), "Scale point sizes by this factor")
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
  opts.face_labels          =  FileSystem::resolve(opts.face_labels);
  opts.accentuate_features  =  (vm.count("emph-features") > 0);
  opts.color_cube_features  =  (vm.count("color-cube") > 0);
  opts.show_normals         =  (vm.count("normals") > 0);
  opts.show_graph           =  (vm.count("graph") > 0);
  opts.flat                 =  (vm.count("flat") > 0);
  opts.fancy_points         =  (vm.count("fancy-points") > 0);
  opts.fancy_colors         =  (vm.count("fancy-colors") > 0);

  if (!s_bg_color.empty())
  {
    std::stringstream ss;
    ss << std::hex << s_bg_color;

    uint32 argb;
    ss >> argb;

    opts.bg_plain = true;
    opts.bg_color = ColorRGB::fromARGB(argb);
  }
  else
  {
    opts.bg_plain = false;
    opts.bg_color = ColorRGB::black();
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
  TheaArray<std::string> plugin_dirs; plugin_dirs.push_back(opts.plugin_dir);
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
