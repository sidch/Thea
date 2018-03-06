//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Stanford University
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

#include "Application.hpp"
#include "Algorithms/EigenSolver.hpp"
#include "Algorithms/LinearSolver.hpp"
#include "Algorithms/NumericalOptimizer.hpp"
#include "Graphics/RenderSystem.hpp"
#include "DynLib.hpp"
#include "FilePath.hpp"
#include "FileSystem.hpp"
#include "Plugin.hpp"
#include "System.hpp"

#if defined(THEA_WINDOWS)

#  include <windows.h>

#elif defined(THEA_LINUX)

#  include <unistd.h>

#elif defined(THEA_OSX)

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
    GetModuleFileNameA(NULL, path, sizeof(path));
  }
#elif defined(THEA_OSX)
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
    // Ensure proper NULL termination
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
Application::getPluginPath(std::string const & plugin_name, TheaArray<std::string> const * plugin_dirs)
{
  std::string app_dir = FilePath::parent(FileSystem::resolve(Application::programPath()));
  std::string app_plugin_dir = FilePath::concat(FilePath::parent(app_dir), "lib");

  TheaArray<std::string> search_dirs;
  search_dirs.push_back(app_dir);             // the application directory has highest priority
  search_dirs.push_back(app_plugin_dir);      // the relative plugins directory is next

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
#elif THEA_OSX
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
