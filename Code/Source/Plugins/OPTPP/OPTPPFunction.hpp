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

#ifndef __Thea_Algorithms_OPTPPFunction_hpp__
#define __Thea_Algorithms_OPTPPFunction_hpp__

#include "OPTPPCommon.hpp"
#include "../../Algorithms/IScalarFunction.hpp"
#include <NLP.h>

namespace Thea {
namespace Algorithms {

// Convert a scalar function to OPT++ form.
OPTPP::NLPBase * toOPTPPFunction(IScalarFunction const & f, double const * init_pt = nullptr);

} // namespace Algorithms
} // namespace Thea

#endif
