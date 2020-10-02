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

#ifndef __Thea_Algorithms_IAnalyticD2ScalarFunction_hpp__
#define __Thea_Algorithms_IAnalyticD2ScalarFunction_hpp__

#include "../AddressableMatrix.hpp"
#include "../IMatrix.hpp"
#include "IAnalyticD1ScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for a function that maps each point in R^n to a real number and has analytic first (gradient) and second (Hessian)
 * derivatives.
 */
class THEA_API IAnalyticD2ScalarFunction : public virtual IAnalyticD1ScalarFunction
{
  public:
    THEA_DECL_SMART_POINTERS(IAnalyticD2ScalarFunction)

    /** Does the function have a sparse Hessian matrix? */
    virtual bool THEA_ICALL hasSparseHessian() const = 0;

    /** Compute the Hessian of the function at a given point. Returns true on success, false on error. */
    virtual void THEA_ICALL hessianAt(float32 const * p, IMatrix<float32> * result) const = 0;

    /** Compute the Hessian of the function at a given point. Returns true on success, false on error. */
    virtual void THEA_ICALL hessianAt(float64 const * p, IMatrix<float64> * result) const = 0;

}; // class IAnalyticD2ScalarFunction

} // namespace Algorithms
} // namespace Thea

#endif
