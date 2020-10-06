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

#include "OPTPPNumericalOptimizer.hpp"
#include "OPTPPConstraints.hpp"
#include <NLF.h>
#include <CompoundConstraint.h>

namespace Thea {
namespace Algorithms {

OPTPPNumericalOptimizer::OPTPPNumericalOptimizer(std::string const & name_)
: NamedObject(name_)
{}

bool
OPTPPNumericalOptimizer::minimize(IScalarFunction const & objective, double const * hint, Options const & options)
{
  // TODO

  has_solution = false;
  solution.clear();

  return has_solution;
}

OPTPPNumericalOptimizerFactory::~OPTPPNumericalOptimizerFactory()
{
  destroyAllNumericalOptimizers();
}

INumericalOptimizer *
OPTPPNumericalOptimizerFactory::createNumericalOptimizer(char const * name)
{
  OPTPPNumericalOptimizer * ls = new OPTPPNumericalOptimizer(name);
  optimizers.insert(ls);
  return ls;
}

void
OPTPPNumericalOptimizerFactory::destroyNumericalOptimizer(INumericalOptimizer * optimizer)
{
  optimizers.erase(optimizer);
  delete optimizer;
}

void
OPTPPNumericalOptimizerFactory::destroyAllNumericalOptimizers()
{
  for (NumericalOptimizerSet::iterator li = optimizers.begin(); li != optimizers.end(); ++li)
    delete *li;

  optimizers.clear();
}

} // namespace Algorithms
} // namespace Thea
