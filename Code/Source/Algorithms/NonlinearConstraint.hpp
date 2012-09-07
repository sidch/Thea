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

#ifndef __Thea_Algorithms_NonlinearConstraint_hpp__
#define __Thea_Algorithms_NonlinearConstraint_hpp__

#include "../Common.hpp"
#include "ScalarFunction.hpp"

namespace Thea {
namespace Algorithms {

/** A nonlinear optimization constraint. */
class NonlinearConstraint
{
  public:
    THEA_DEF_POINTER_TYPES(NonlinearConstraint, shared_ptr, weak_ptr)

    /**
     * Constructor. Creates a nonlinear constraint of the form f(<b>x</b>) <em>op</em> rhs, where <em>op</em> is a comparison
     * operator.
     *
     * @param f_ The constraint function f.
     * @param compare_ Comparison operator <em>op</em>.
     * @param rhs_ Right-hand side constant rhs.
     */
    NonlinearConstraint(ScalarFunction::Ptr f_, CompareOp compare_, double rhs_)
    : f(f_), compare(compare_), rhs(rhs_)
    {
      alwaysAssertM(f, "NonlinearConstraint: Function cannot be null");
    }

    /** Get the number of dimensions of the linear constraint. */
    int numDimensions() const { return f->numDimensions(); }

    /** Get the function. */
    ScalarFunction::ConstPtr getFunction() const { return f; }

    /** Get the comparison operator. */
    CompareOp getCompareOp() const { return compare; }

    /** Get the right-hand side. */
    double getRHS() const { return rhs; }

  private:
    ScalarFunction::Ptr f;  ///< Constraint function.
    CompareOp compare;      ///< Comparison operator.
    double rhs;             ///< Right-hand side constant.

}; // class NonlinearConstraint

} // namespace Algorithms
} // namespace Thea

#endif
