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

#include "OPTPPPlugin.hpp"
#include "OPTPPNumericalOptimizer.hpp"

namespace Thea {

static Algorithms::OPTPPPlugin * plugin = NULL;

/** DLL start routine. Installs plugin. */
extern "C" THEA_OPTPP_API Plugin *
dllStartPlugin()
{
  plugin = new Algorithms::OPTPPPlugin();
  PluginManager::install(plugin);
  return plugin;
}

/** DLL stop routine. Uninstalls plugin. */
extern "C" THEA_OPTPP_API void
dllStopPlugin()
{
  PluginManager::uninstall(plugin);
  delete plugin;
}

namespace Algorithms {

static std::string const OPTPP_PLUGIN_NAME              =  "OPT++ NumericalOptimizer";
static std::string const OPTPP_NUMERICALOPTIMIZER_NAME  =  "OPT++";

OPTPPPlugin::OPTPPPlugin()
: NamedObject(OPTPP_PLUGIN_NAME), factory(NULL), started(false)
{}

OPTPPPlugin::~OPTPPPlugin()
{
  uninstall();
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

    NumericalOptimizerManager::installFactory(OPTPP_NUMERICALOPTIMIZER_NAME, factory);
    started = true;
  }
}

void
OPTPPPlugin::shutdown()
{
  if (started)
  {
    factory->destroyAllNumericalOptimizers();
    NumericalOptimizerManager::uninstallFactory(OPTPP_NUMERICALOPTIMIZER_NAME);
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
    factory = NULL;
  }
}

} // namespace Algorithms

} // namespace Thea
