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

#include "ARPACKPlugin.hpp"
#include "ARPACKEigenSolver.hpp"

namespace Thea {

static Algorithms::ARPACKPlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_ARPACK_API Plugin *
dllStartPlugin(FactoryRegistry * registry_)
{
  plugin = new Algorithms::ARPACKPlugin(registry_);
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

ARPACKPlugin::ARPACKPlugin(FactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(ARPACK_PLUGIN_NAME) + ": Factory registry must be non-null");
}

ARPACKPlugin::~ARPACKPlugin()
{
  uninstall();
}

char const *
ARPACKPlugin::getName() const
{
  return ARPACK_PLUGIN_NAME;
}

void
ARPACKPlugin::install()
{}

void
ARPACKPlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new ARPACKEigenSolverFactory;

    registry->addEigenSolverFactory(ARPACK_EIGENSOLVER_NAME, factory);
    started = true;
  }
}

void
ARPACKPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllEigenSolvers();

    registry->removeEigenSolverFactory(ARPACK_EIGENSOLVER_NAME);
    started = false;
  }
}

void
ARPACKPlugin::uninstall()
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
