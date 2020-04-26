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
// First version: 2011
//
//============================================================================

#ifndef __Thea_Algorithms_ScalarConstraint_hpp__
#define __Thea_Algorithms_ScalarConstraint_hpp__

#include "../Common.hpp"
#include "ScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for an optimization constraint of the form <code>f(x) op y</code>, where <code>x</code> is a vector, <code>y</code>
 * is a scalar, <code>f</code> is a function, and <code>op</code> is a comparison operator.
 */
class ScalarConstraint
{
  public:
    THEA_DECL_SMART_POINTERS(ScalarConstraint)

    /** Get the function. */
    virtual ScalarFunction const & getFunction() const = 0;

    /** Get the comparison operator, corresponding to a value of the CompareOp enum class. */
    virtual int32 getCompareOp() const = 0;

    /** Get the right-hand side. */
    virtual float64 getRHS() const = 0;

}; // class ScalarConstraint

} // namespace Algorithms
} // namespace Thea

#endif
