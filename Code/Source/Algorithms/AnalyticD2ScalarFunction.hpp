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

#ifndef __Thea_Algorithms_AnalyticD2ScalarFunction_hpp__
#define __Thea_Algorithms_AnalyticD2ScalarFunction_hpp__

#include "../AddressableMatrix.hpp"
#include "../AbstractMatrix.hpp"
#include "AnalyticD1ScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for a function that maps each point in R^n to a real number and has analytic first (gradient) and second (Hessian)
 * derivatives.
 */
class THEA_API AnalyticD2ScalarFunction : public AnalyticD1ScalarFunction
{
  public:
    THEA_DECL_SMART_POINTERS(AnalyticD2ScalarFunction)

    /** Does the function have a sparse Hessian matrix? */
    virtual bool hasSparseHessian() const = 0;

    /** Compute the Hessian of the function at a given point. Returns true on success, false on error. */
    virtual void hessianAt(float32 const * p, AbstractMatrix<float32> & result) const = 0;

    /** Compute the Hessian of the function at a given point. Returns true on success, false on error. */
    virtual void hessianAt(float64 const * p, AbstractMatrix<float64> & result) const = 0;

}; // class AnalyticD2ScalarFunction

} // namespace Algorithms
} // namespace Thea

#endif
