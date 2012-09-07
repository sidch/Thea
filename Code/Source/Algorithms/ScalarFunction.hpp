//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_ScalarFunction_hpp__
#define __Thea_Algorithms_ScalarFunction_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for a function that maps each point in R^n to a real number.
 *
 * <b>IMPORTANT:</b> Always inherit virtually from this class, i.e. use:
 *
 * \code
 * class MyFunc : public virtual ScalarFunction
 * {
 *   ...
 * };
 * \endcode
 *
 * This also means that the a derived class in an inheritance hierarchy involving ScalarFunction must explicitly call the
 * ScalarFunction constructor.
 */
class THEA_API ScalarFunction
{
  public:
    THEA_DEF_POINTER_TYPES(ScalarFunction, shared_ptr, weak_ptr)

    /** Constructor. */
    ScalarFunction(int ndims_) : ndims(ndims_)
    { alwaysAssertM(ndims > 0, "ScalarFunction: Number of dimensions must be positive"); }

    /** Destructor. */
    virtual ~ScalarFunction() {}

    /** Get the number of dimensions of the function domain. */
    int numDimensions() const { return ndims; }

    /** Compute the value of the function at a given point. */
    virtual double valueAt(float const * p) const = 0;

    /** Compute the value of the function at a given point. */
    virtual double valueAt(double const * p) const = 0;

  private:
    int ndims;  ///< Number of dimensions of domain.

}; // class ScalarFunction

} // namespace Algorithms
} // namespace Thea

#endif
