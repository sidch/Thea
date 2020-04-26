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
// First version: 2009
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
    static std::string getPluginPath(std::string const & plugin_name, Array<std::string> const * plugin_dirs = nullptr);

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
