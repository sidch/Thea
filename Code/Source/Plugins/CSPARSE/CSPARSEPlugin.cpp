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

#include "CSPARSEPlugin.hpp"
#include "CSPARSELinearSolver.hpp"

namespace Thea {

static Algorithms::CSPARSEPlugin * plugin = NULL;

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
: registry(registry_), factory(NULL), started(false)
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
    factory = NULL;
  }
}

} // namespace Algorithms

} // namespace Thea
