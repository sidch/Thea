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

#ifndef __Thea_Algorithms_OPTPPConstraints_hpp__
#define __Thea_Algorithms_OPTPPConstraints_hpp__

#include "OPTPPCommon.hpp"
#include "../../Algorithms/LinearConstraint.hpp"
#include "../../Algorithms/NonlinearConstraint.hpp"
#include <LinearConstraint.h>
#include <NonLinearConstraint.h>

namespace Thea {
namespace Algorithms {

// Convert a linear constraint to OPT++ form.
OPTPP::LinearConstraint * toOPTPPConstraint(LinearConstraint const & constraint);

// Convert a nonlinear constraint to OPT++ form.
OPTPP::NonLinearConstraint * toOPTPPConstraint(NonlinearConstraint const & constraint);

} // namespace Algorithms
} // namespace Thea

#endif
