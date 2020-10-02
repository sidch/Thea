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
// First version: 2013
//
//============================================================================

#include "Application.hpp"
#include "Algorithms/IEigenSolver.hpp"
#include "Algorithms/ILinearSolver.hpp"
#include "Algorithms/INumericalOptimizer.hpp"
#include "Graphics/IRenderSystem.hpp"
#include "DynLib.hpp"
#include "FilePath.hpp"
#include "FileSystem.hpp"
#include "IPlugin.hpp"
#include "System.hpp"

#if defined(THEA_WINDOWS)

#  include <windows.h>

#elif defined(THEA_LINUX)

#  include <unistd.h>

#elif defined(THEA_MAC)

#  include <mach-o/dyld.h>
#  include <stdlib.h>

#endif

namespace Thea {

Application::Globals Application::globals;

Application::Globals::Globals()
{
  eigen_solver_mgr = new Algorithms::EigenSolverManager;
  linear_solver_mgr = new Algorithms::LinearSolverManager;
  numerical_optimizer_mgr = new Algorithms::NumericalOptimizerManager;
  render_system_mgr = new Graphics::RenderSystemManager;
  dynlib_mgr = new DynLibManager;
  plugin_mgr = new PluginManager;
}

Application::Globals::~Globals()
{
  delete plugin_mgr;
  delete dynlib_mgr;
  delete render_system_mgr;
  delete numerical_optimizer_mgr;
  delete linear_solver_mgr;
  delete eigen_solver_mgr;
}

std::string
Application::programPath()
{
  char path[8192];

#ifdef THEA_WINDOWS
  {
    GetModuleFileNameA(nullptr, path, sizeof(path));
  }
#elif defined(THEA_MAC)
  {
    char unresolved_path[8192];
    uint32_t size = (uint32_t)sizeof(unresolved_path);
    if (_NSGetExecutablePath(unresolved_path, &size) != 0)
      throw Error("Application: Executable path is longer than buffer size");

    if (!realpath(unresolved_path, path))
      throw Error("Application: Could not resolve executable path");
  }
#else
  {
    int ret = readlink("/proc/self/exe", path, sizeof(path));

    // In case of an error, leave the handling up to the caller
    if (ret == -1)
      return "";

    debugAssertM((int)sizeof(path) > ret, "Application: String too short to store current program path");
    // Ensure proper nullptr termination
    path[ret] = 0;
  }
#endif

  return path;
}

void
Application::setResourceArchive(std::string const & path)
{
  if (!path.empty())
  {
    if (!FileSystem::directoryExists(path))
      throw Error("Resource archive '" + path + "' does not exist or is not a valid directory");

    _resourceArchive() = FileSystem::resolve(path);

    THEA_DEBUG << "Resource archive set to '" << _resourceArchive() << '\'';
  }
}

std::string
Application::getResourcePath(std::string const & resource_name)
{
  std::string path = FilePath::concat(_resourceArchive(), resource_name);
  return FileSystem::exists(path) ? path : std::string();
}

std::string
Application::getPluginPath(std::string const & plugin_name, Array<std::string> const * plugin_dirs)
{
  std::string app_dir = FilePath::parent(FileSystem::resolve(Application::programPath()));
  std::string app_plugin_dir1 = FilePath::concat(FilePath::parent(app_dir), "lib");     // in parent
  std::string app_plugin_dir2 = FilePath::concat(FilePath::parent(app_dir), "../lib");  // in grandparent

  Array<std::string> search_dirs;
  search_dirs.push_back(app_dir);             // the application directory has highest priority
  search_dirs.push_back(app_plugin_dir1);     // the parent directory is next
  search_dirs.push_back(app_plugin_dir2);     // the grandparent directory is next

  if (plugin_dirs)                            // user-supplied search directories are next
    search_dirs.insert(search_dirs.end(), plugin_dirs->begin(), plugin_dirs->end());

  search_dirs.push_back("/usr/lib");          // system directories follow
  search_dirs.push_back("/usr/local/lib");
  search_dirs.push_back("/sw/lib");           //  -- Fink
  search_dirs.push_back("/opt/local/lib");    //  -- Homebrew

  std::string dir = FilePath::parent(plugin_name);
  std::string basename = FilePath::completeBaseName(plugin_name);
  std::string ext = FilePath::extension(plugin_name);

#ifndef THEA_WINDOWS
  if (!beginsWith(basename, "lib"))
    basename = "lib" + basename;
#endif

#ifdef THEA_WINDOWS
  if (ext != "dll")
  {
    if (!ext.empty()) basename = basename + "." + ext;  // e.g. we got something like libGL.3
    ext = "dll";
  }
#elif THEA_MAC
  if (ext != "dylib")
  {
    if (!ext.empty()) basename = basename + "." + ext;
    ext = "dylib";
  }
#else
  if (ext != "so")
  {
    if (!ext.empty()) basename = basename + "." + ext;
    ext = "so";
  }
#endif

  std::string debug_filename    =  basename + "d." + ext;
  std::string release_filename  =  basename + "." + ext;

  for (size_t i = 0; i < search_dirs.size(); ++i)
  {
    if (search_dirs[i].empty())  // if you need to specify the current working directory, use "."
      continue;

    std::string debug_path    =  FilePath::concat(search_dirs[i], debug_filename);
    std::string release_path  =  FilePath::concat(search_dirs[i], release_filename);

#ifdef THEA_DEBUG_BUILD

    // Debug preferred over release
    if (FileSystem::fileExists(debug_path))
      return FileSystem::resolve(debug_path);
    else if (FileSystem::fileExists(release_path))
      return FileSystem::resolve(release_path);

#else // THEA_DEBUG_BUILD

    // Release preferred over debug
    if (FileSystem::fileExists(release_path))
      return FileSystem::resolve(release_path);
    else if (FileSystem::fileExists(debug_path))
      return FileSystem::resolve(debug_path);

#endif // THEA_DEBUG_BUILD
  }

  return std::string();
}

std::string &
Application::_resourceArchive()
{
  static std::string resource_archive = FilePath::parent(programPath());
  return resource_archive;
}

} // namespace Thea
