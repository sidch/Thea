//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Application_hpp__
#define __Thea_Application_hpp__

#include "Common.hpp"

namespace Thea {

// Forward declarations
class DynLibManager;
class PluginManager;

namespace Algorithms {

class EigenSolverManager;
class LinearSolverManager;
class NumericalOptimizerManager;

} // namespace Graphics

namespace Graphics {

class RenderSystemManager;

} // namespace Graphics

/**
 * Access and modify global properties of the current application.
 */
class THEA_API Application
{
  public:
    /** Get the path to the current program. */
    static std::string programPath();

    /**
     * Get the archive in which the application will look for resources (by default the directory containing the executable).
     */
    static std::string const & getResourceArchive() { return _resourceArchive(); }

    /**
     * Set the archive in which the application will look for resources.
     */
    static void setResourceArchive(std::string const & path);

    /** Get the fully qualified path to a resource, or the empty string if it cannot be found. */
    static std::string getResourcePath(std::string const & resource_name);

    /**
     * Get the fully qualified path to a plugin, or the empty string if it cannot be found. The search order is:
     * - Directory containing executable
     * - "lib" folder in parent of above directory
     * - User-supplied \a plugin_dirs, if not null
     * - %System directories, e.g. <code>/usr/lib</code> and <code>/usr/local/lib</code>
     *
     * @param plugin_name Name of the plugin, with or without extension, relative path, leading "lib", etc.
     * @param plugin_dirs Optional additional set of search directories.
     */
    static std::string getPluginPath(std::string const & plugin_name, TheaArray<std::string> const * plugin_dirs = NULL);

    /** Get the global plugin manager. */
    static PluginManager & getPluginManager() { return *globals.plugin_mgr; }

    /** Get the global dynamic library manager. */
    static DynLibManager & getDynLibManager() { return *globals.dynlib_mgr; }

    /** Get the global eigensolver manager. */
    static Algorithms::EigenSolverManager & getEigenSolverManager() { return *globals.eigen_solver_mgr; }

    /** Get the global linear solver manager. */
    static Algorithms::LinearSolverManager & getLinearSolverManager() { return *globals.linear_solver_mgr; }

    /** Get the global numerical optimizer manager. */
    static Algorithms::NumericalOptimizerManager & getNumericalOptimizerManager() { return *globals.numerical_optimizer_mgr; }

    /** Get the global rendersystem manager. */
    static Graphics::RenderSystemManager & getRenderSystemManager() { return *globals.render_system_mgr; }

  private:
    /** Encapsulates global objects. */
    struct Globals
    {
      PluginManager                          *  plugin_mgr;
      DynLibManager                          *  dynlib_mgr;
      Algorithms::EigenSolverManager         *  eigen_solver_mgr;
      Algorithms::LinearSolverManager        *  linear_solver_mgr;
      Algorithms::NumericalOptimizerManager  *  numerical_optimizer_mgr;
      Graphics::RenderSystemManager          *  render_system_mgr;

      /** Constructor. */
      Globals();

      /** Destructor. */
      ~Globals();

    }; // struct Globals

    /**
     * Wraps the resource archive variable to ensure that a) we don't need a separate implementation file and b) the variable is
     * initialized to the default value on first call and not earlier.
     */
    static std::string & _resourceArchive();

    static Globals globals;  ///< Global objects.

}; // class Application

} // namespace Thea

#endif
