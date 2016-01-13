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
#include "../../FileSystem.hpp"
#include "../../Plugin.hpp"
#include "../../Graphics/RenderSystem.hpp"
#include <QApplication>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QGLFormat>
#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>

namespace Browse3D {

App::Options::Options()
: bg_plain(false), bg_color(ColorRGB::black()), two_sided(true), accentuate_features(false), show_graph(false)
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
      << "\n  overlays = { " << opts.overlays.join(", ") << " }"
      << "\n  features = " << opts.features
      << "\n  bg-plain = " << opts.bg_plain
      << "\n  bg-color = " << opts.bg_color.toString()
      << "\n  two-sided = " << opts.two_sided
      << "\n  flat = " << opts.flat
      << "\n  accentuate-features = " << opts.accentuate_features
      << "\n  color-cube-features = " << opts.color_cube_features
      << "\n  show-graph = " << opts.show_graph
      << "\n  fancy-points = " << opts.fancy_points
      << "\n  fancy-colors = " << opts.fancy_colors
      << "\n  point-scale = " << opts.point_scale
      << '\n';

  return oss.str();
}

bool
App::init(int argc, char * argv[])
{
  if (!parseOptions(argc, argv))
    return false;

  qDebug() << "Started Browse3D";
  std::string opt_str = optsToString();
  qDebug() << opt_str;

  createMainWindow();

  // Load plugins and create a rendersystem
  loadPlugins();
  createRenderSystem();
  main_window->update();

  return true;
}

bool
parseModel(std::string const & str, QString & path, AffineTransform3 & transform)
{
  TheaArray<std::string> fields;
  stringSplit(str, '|', fields);
  if (fields.empty() || fields.size() > 2)
  {
    THEA_ERROR << "Could not parse model: " << str;
    return false;
  }

  path = toQString(trimWhitespace(fields[0]));

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
  namespace po = boost::program_options;

  static std::string const usage("Usage: Browse3D [options] [model]");

  std::string conf_file;

  // po::options_description hidden;
  // hidden.add_options()("hidden-option", po::value<std::string>(&hidden_option), "");

  std::string app_dir(QFile::encodeName(QApplication::applicationDirPath()).data());
  std::string def_plugin_dir = getFullPath(app_dir, "../lib");
#ifdef _MSC_VER
  // Visual Studio puts executables in Debug|Release subdirectory
  std::string def_resource_dir = getFullPath(app_dir, "../../../../Resources");
#else
  std::string def_resource_dir = getFullPath(app_dir, "../../../Resources");
#endif

  std::string s_plugin_dir;
  std::string s_resource_dir;
  std::string s_working_dir;
  std::string s_model;
  std::vector<std::string> s_overlays;
  std::string s_features;
  std::string s_bg_color;

  po::options_description visible("Allowed options");
  visible.add_options()
          ("help,h",               "Print this help message")
          ("version,v",            "Print the program version")
          ("conf",                 po::value<std::string>(&conf_file)->default_value("Browse3D.conf"), "Configuration file (overridden by duplicate cmdline options)")
          ("plugin-dir",           po::value<std::string>(&s_plugin_dir)->default_value(def_plugin_dir), "Plugins directory")
          ("resource-dir",         po::value<std::string>(&s_resource_dir)->default_value(def_resource_dir), "Resources directory")
          ("working-dir",          po::value<std::string>(&s_working_dir)->default_value("."), "Working directory")
          ("model",                po::value<std::string>(&s_model), "Model to load on startup")
          ("overlay",              po::value< std::vector<std::string> >(&s_overlays), "Overlay model(s) to load on startup")
          ("features,f",           po::value<std::string>(&s_features), "Directory or file containing features to load")
          ("emph-features,e",      "Make feature distributions easier to view")
          ("color-cube,3",         "Map 0-centered 3D feature sets to RGB color-cube, if --emph-features")
          ("graph,g",              "Show point adjacency graph")
          ("bg",                   po::value<std::string>(&s_bg_color), "Background color")
          ("two-sided",            po::value<bool>(&opts.two_sided)->default_value(true), "Use two-sided lighting?")
          ("flat",                 po::value<bool>(&opts.flat)->default_value(true), "Flat shade all meshes?")
          ("fancy-points",         "Draw points as shaded spheres?")
          ("fancy-colors,c",       "Color points by a function of position?")
          ("point-scale,s",        po::value<Real>(&opts.point_scale)->default_value(1), "Scale point sizes by this factor")
  ;

  po::options_description desc;
  desc.add(visible) /* .add(hidden) */ ;

  if (argc < 1)
  {
    qDebug() << usage;
    std::cerr << visible;  // should be intercepted by out
    return 0;
  }

  po::positional_options_description pdesc;
  pdesc.add("model", 1);
  pdesc.add("overlay", -1);

  // Read cmdline options first (overrides conflicting config file values)
  po::parsed_options cmdline_parsed = po::basic_command_line_parser<char>(argc, argv).options(desc).positional(pdesc).run();
  po::variables_map vm;
  po::store(cmdline_parsed, vm);

  // Now read the config file, if it is found
  if (vm.count("conf") > 0) conf_file = vm["conf"].as<std::string>();
  if (QFile::exists(QFile::decodeName(conf_file.c_str())))
  {
    qDebug() << "Reading options from config file:" << conf_file;

    std::ifstream conf_in(conf_file.c_str());
    po::parsed_options conf_file_parsed = po::parse_config_file(conf_in, desc);
    po::store(conf_file_parsed, vm);
  }

  po::notify(vm);

  bool quit = false;

  if (vm.count("version") > 0)
  {
    qDebug() << "Browse3D version 1.1";
    qDebug() << "Siddhartha Chaudhuri/Stanford University/Princeton University, 2013";
    quit = true;
  }

  if (vm.count("help") > 0)
  {
    if (quit) qDebug();
    qDebug() << usage;
    std::cout << visible;  // should be intercepted by qDebug()
    quit = true;
  }

  if (quit)
    return false;

  if (!s_plugin_dir.empty())    opts.plugin_dir    =  QDir(QFile::decodeName(s_plugin_dir.c_str())).canonicalPath();
  if (!s_resource_dir.empty())  opts.resource_dir  =  QDir(QFile::decodeName(s_resource_dir.c_str())).canonicalPath();
  if (!s_working_dir.empty())   opts.working_dir   =  QDir(QFile::decodeName(s_working_dir.c_str())).canonicalPath();

  if (!s_model.empty())
  {
    if (!parseModel(s_model, opts.model, opts.model_transform))
      return false;
  }

  opts.overlays.clear();
  opts.overlay_transforms.clear();
  for (size_t i = 0; i < s_overlays.size(); ++i)
  {
    QString path;
    AffineTransform3 transform;
    if (!parseModel(s_overlays[i], path, transform))
      return false;

    opts.overlays << path;
    opts.overlay_transforms.push_back(transform);
  }

  if (!s_features.empty()) opts.features = toQString(s_features);
  opts.accentuate_features  =  (vm.count("emph-features") > 0);
  opts.color_cube_features  =  (vm.count("color-cube") > 0);
  opts.show_graph           =  (vm.count("graph") > 0);
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

  Application::setResourceArchive(s_resource_dir);

  return true;
}

void
App::createMainWindow()
{
  // Enable antialiasing
  QGLFormat glf = QGLFormat::defaultFormat();
  glf.setSampleBuffers(true);
  glf.setSamples(4);
  QGLFormat::setDefaultFormat(glf);

  // Create the main window, and hence a rendering context
  main_window = new MainWindow;
  main_window->init();
  main_window->raise();
  main_window->activateWindow();
  main_window->show();
}

void
App::loadPlugins()
{
  // Try to load the OpenGL plugin
#ifdef THEA_DEBUG_BUILD
  std::string s_plugin_dir         =  toStdString(opts.plugin_dir);

#ifdef THEA_WINDOWS
  std::string debug_plugin_path    =  getFullPath(s_plugin_dir, "TheaPluginGLd");
  std::string release_plugin_path  =  getFullPath(s_plugin_dir, "TheaPluginGL");
#else
  std::string debug_plugin_path    =  getFullPath(s_plugin_dir, "libTheaPluginGLd");
  std::string release_plugin_path  =  getFullPath(s_plugin_dir, "libTheaPluginGL");
#endif

#ifdef THEA_WINDOWS
  std::string debug_plugin_path_ext = debug_plugin_path + ".dll";
#elif THEA_OSX
  std::string debug_plugin_path_ext = debug_plugin_path + ".dylib";
#else
  std::string debug_plugin_path_ext = debug_plugin_path + ".so";
#endif

  std::string plugin_path = FileSystem::exists(debug_plugin_path_ext) ? debug_plugin_path : release_plugin_path;
#else

#ifdef THEA_WINDOWS
  std::string plugin_path = getFullPath(toStdString(opts.plugin_dir), "TheaPluginGL");
#else
  std::string plugin_path = getFullPath(toStdString(opts.plugin_dir), "libTheaPluginGL");
#endif

#endif

  qDebug() << "Loading OpenGL plugin:" << plugin_path;
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

void
App::cleanup()
{
  // QObjects without parents must be manually deleted

  if (main_window && !main_window->parent())
  {
    delete main_window;
    main_window = NULL;
  }

  if (render_system_factory)
    render_system_factory->destroyRenderSystem(render_system);

  Application::getPluginManager().unloadAllPlugins();
}

} // namespace Browse3D
