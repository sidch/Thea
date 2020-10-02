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

#ifndef __Thea_Algorithms_IScalarConstraint_hpp__
#define __Thea_Algorithms_IScalarConstraint_hpp__

#include "../Common.hpp"
#include "IScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for an optimization constraint of the form <code>f(x) op y</code>, where <code>x</code> is a vector, <code>y</code>
 * is a scalar, <code>f</code> is a function, and <code>op</code> is a comparison operator.
 */
class IScalarConstraint
{
  public:
    THEA_DECL_SMART_POINTERS(IScalarConstraint)

    /** Destructor. */
    virtual ~IScalarConstraint() = 0;

    /** Get the function. */
    virtual IScalarFunction const * THEA_ICALL getFunction() const = 0;

    /** Get the comparison operator, corresponding to a value of the CompareOp enum class. */
    virtual int32 THEA_ICALL getCompareOp() const = 0;

    /** Get the right-hand side. */
    virtual float64 THEA_ICALL getRhs() const = 0;

}; // class IScalarConstraint

inline IScalarConstraint::~IScalarConstraint() {}

} // namespace Algorithms
} // namespace Thea

#endif
