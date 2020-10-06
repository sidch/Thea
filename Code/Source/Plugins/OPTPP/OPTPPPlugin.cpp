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

#include "OPTPPPlugin.hpp"
#include "OPTPPNumericalOptimizer.hpp"

namespace Thea {

static Algorithms::OPTPPPlugin * plugin = nullptr;

/** DLL start routine. Installs plugin. */
extern "C" THEA_OPTPP_API IPlugin *
dllStartPlugin(IFactoryRegistry * registry_)
{
  plugin = new Algorithms::OPTPPPlugin(registry_);
  return plugin;
}

/** DLL stop routine. Uninstalls plugin. */
extern "C" THEA_OPTPP_API void
dllStopPlugin()
{
  delete plugin;
}

namespace Algorithms {

static char const * OPTPP_PLUGIN_NAME              =  "OPT++ INumericalOptimizer";
static char const * OPTPP_NUMERICALOPTIMIZER_NAME  =  "OPT++";

OPTPPPlugin::OPTPPPlugin(IFactoryRegistry * registry_)
: registry(registry_), factory(nullptr), started(false)
{
  alwaysAssertM(registry, std::string(OPTPP_PLUGIN_NAME) + ": Factory registry must be non-null");
}

OPTPPPlugin::~OPTPPPlugin()
{
  uninstall();
}

char const *
OPTPPPlugin::getName() const
{
  return OPTPP_PLUGIN_NAME;
}

void
OPTPPPlugin::install()
{}

void
OPTPPPlugin::startup()
{
  if (!started)
  {
    if (!factory)
      factory = new OPTPPNumericalOptimizerFactory;

    registry->addNumericalOptimizerFactory(OPTPP_NUMERICALOPTIMIZER_NAME, factory);
    started = true;
  }
}

void
OPTPPPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllNumericalOptimizers();

    registry->removeNumericalOptimizerFactory(OPTPP_NUMERICALOPTIMIZER_NAME);
    started = false;
  }
}

void
OPTPPPlugin::uninstall()
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
