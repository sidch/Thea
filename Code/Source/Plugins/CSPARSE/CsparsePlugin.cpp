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

#include "CsparsePlugin.hpp"
#include "CsparseLinearSolver.hpp"

namespace Thea {

static Algorithms::CsparsePlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_CSPARSE_API IPlugin *
dllStartPlugin(IFactoryRegistry * registry_)
{
  plugin = new Algorithms::CsparsePlugin(registry_);
  return plugin;
}

/** DLL stop routine. Uninstalls plugin. */
extern "C" THEA_CSPARSE_API void
dllStopPlugin()
{
  delete plugin;
}

namespace Algorithms {

static char const * CSPARSE_PLUGIN_NAME        =  "CSPARSE LinearSolver";
static char const * CSPARSE_LINEARSOLVER_NAME  =  "CSPARSE";

CsparsePlugin::CsparsePlugin(IFactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(CSPARSE_PLUGIN_NAME) + ": Factory registry must be non-null");
}

CsparsePlugin::~CsparsePlugin()
{
  uninstall();
}

char const *
CsparsePlugin::getName() const
{
  return CSPARSE_PLUGIN_NAME;
}

void
CsparsePlugin::install()
{}

void
CsparsePlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new CsparseLinearSolverFactory;

    registry->addLinearSolverFactory(CSPARSE_LINEARSOLVER_NAME, factory);
    started = true;
  }
}

void
CsparsePlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllLinearSolvers();

    registry->removeLinearSolverFactory(CSPARSE_LINEARSOLVER_NAME);
    started = false;
  }
}

void
CsparsePlugin::uninstall()
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
