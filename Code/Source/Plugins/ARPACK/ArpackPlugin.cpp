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

#include "ArpackPlugin.hpp"
#include "ArpackEigenSolver.hpp"

namespace Thea {

static Algorithms::ArpackPlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_ARPACK_API IPlugin *
dllStartPlugin(IFactoryRegistry * registry_)
{
  plugin = new Algorithms::ArpackPlugin(registry_);
  return plugin;
}

/** DLL stop routine. Uninstalls plugin. */
extern "C" THEA_ARPACK_API void
dllStopPlugin()
{
  delete plugin;
}

namespace Algorithms {

static char const * ARPACK_PLUGIN_NAME       =  "ARPACK EigenSolver";
static char const * ARPACK_EIGENSOLVER_NAME  =  "ARPACK";

ArpackPlugin::ArpackPlugin(IFactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(ARPACK_PLUGIN_NAME) + ": Factory registry must be non-null");
}

ArpackPlugin::~ArpackPlugin()
{
  uninstall();
}

char const *
ArpackPlugin::getName() const
{
  return ARPACK_PLUGIN_NAME;
}

void
ArpackPlugin::install()
{}

void
ArpackPlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new ArpackEigenSolverFactory;

    registry->addEigenSolverFactory(ARPACK_EIGENSOLVER_NAME, factory);
    started = true;
  }
}

void
ArpackPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllEigenSolvers();

    registry->removeEigenSolverFactory(ARPACK_EIGENSOLVER_NAME);
    started = false;
  }
}

void
ArpackPlugin::uninstall()
{
  shutdown();  // not currently dependent on presence of other plugins

  if (factory)
  {
    delete factory;
    factory = nullptr;
  }
}

} // namespace Algorithms

} // namespace Thea
