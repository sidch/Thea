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

#include "CSPARSEPlugin.hpp"
#include "CSPARSELinearSolver.hpp"

namespace Thea {

static Algorithms::CSPARSEPlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_CSPARSE_API Plugin *
dllStartPlugin(FactoryRegistry * registry_)
{
  plugin = new Algorithms::CSPARSEPlugin(registry_);
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

CSPARSEPlugin::CSPARSEPlugin(FactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(CSPARSE_PLUGIN_NAME) + ": Factory registry must be non-null");
}

CSPARSEPlugin::~CSPARSEPlugin()
{
  uninstall();
}

char const *
CSPARSEPlugin::getName() const
{
  return CSPARSE_PLUGIN_NAME;
}

void
CSPARSEPlugin::install()
{}

void
CSPARSEPlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new CSPARSELinearSolverFactory;

    registry->addLinearSolverFactory(CSPARSE_LINEARSOLVER_NAME, factory);
    started = true;
  }
}

void
CSPARSEPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllLinearSolvers();

    registry->removeLinearSolverFactory(CSPARSE_LINEARSOLVER_NAME);
    started = false;
  }
}

void
CSPARSEPlugin::uninstall()
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
