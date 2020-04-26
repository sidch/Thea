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

#ifndef __Thea_Algorithms_AnalyticD1ScalarFunction_hpp__
#define __Thea_Algorithms_AnalyticD1ScalarFunction_hpp__

#include "ScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/** Interface for a function that maps each point in R^n to a real number and has an analytic first derivative (gradient). */
class THEA_API AnalyticD1ScalarFunction : public virtual ScalarFunction
{
  public:
    THEA_DECL_SMART_POINTERS(AnalyticD1ScalarFunction)

    /** Compute the gradient of the function at a given point. Returns true on success, false on error. */
    virtual bool gradientAt(float32 const * p, float32 * result) const = 0;

    /** Compute the gradient of the function at a given point. Returns true on success, false on error. */
    virtual bool gradientAt(float64 const * p, float64 * result) const = 0;

    /**
     * If the function has an analytic second derivative, cast this object to the corresponding type. Returns null on failure.
     */
    virtual AnalyticD2ScalarFunction const * asAnalyticD2() const = 0;

}; // class AnalyticD1ScalarFunction

} // namespace Algorithms
} // namespace Thea

#endif
