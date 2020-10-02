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

#ifndef __Thea_Algorithms_IScalarFunction_hpp__
#define __Thea_Algorithms_IScalarFunction_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class IAnalyticD1ScalarFunction;

/** Interface for a function that maps each point in R^n to a real number. */
class THEA_API IScalarFunction
{
  public:
    THEA_DECL_SMART_POINTERS(IScalarFunction)

    /** Destructor. */
    virtual ~IScalarFunction() = 0;

    /** Get the number of dimensions of the function domain. */
    virtual int64 THEA_ICALL dims() const;

    /** Compute the value of the function at a given point. */
    virtual float64 THEA_ICALL valueAt(float32 const * p) const = 0;

    /** Compute the value of the function at a given point. */
    virtual float64 THEA_ICALL valueAt(float64 const * p) const = 0;

    /** If the function has an analytic first derivative, cast this object to the corresponding type. */
    virtual IAnalyticD1ScalarFunction const * THEA_ICALL asAnalyticD1() const = 0;

}; // class IScalarFunction

inline IScalarFunction::~IScalarFunction() {}

} // namespace Algorithms
} // namespace Thea

#endif
