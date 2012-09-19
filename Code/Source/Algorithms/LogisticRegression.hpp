//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_LogisticRegression_hpp__
#define __Thea_Algorithms_LogisticRegression_hpp__

#include "../Common.hpp"
#include "../Array.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class LinearLeastSquares;

/**
 * Solve logistic regression problems. The goal is to fit the logistic curve
 * <em>y</em> = 1 / (1 + exp(-<em>a</em> - <b>b.x</b>)) to data points (<b>x</b>, <em>y</em>). We follow the linearization
 * approach of:
 *
 * J. Mathews, "Bounded population growth: a curve fitting lesson", Math. and Comput. Educ. J. 26(2), pp. 169-176, Spring 1992.
 */
class THEA_API LogisticRegression
{
  public:
    THEA_DEF_POINTER_TYPES(LogisticRegression, shared_ptr, weak_ptr)

    /** Constructor. Sets the number of dimensions of the problem domain (i.e. of the vector <b>x</b>). */
    LogisticRegression(long ndim_);

    /** Destructor. */
    ~LogisticRegression();

    /** Get the number of dimensions of the problem domain. */
    long numDimensions() const;

    /**
     * Add a data point to constrain the curve. This corresponds to a pair of observed values of <b>x</b> and <em>y</em>.
     *
     * @param x The observed vector <b>x</b>.
     * @param y The observed response <em>y</em> for this \a x.
     *
     * @see clearObservations()
     */
    void addObservation(double const * x, double y);

    /**
     * Clear all observations.
     *
     * @see addObservation()
     */
    void clearObservations();

    /** Get the number of observations. */
    long numObservations() const;

    /**
     * Fit the logistic curve to the observed data.
     *
     * @param tolerance The numerical tolerance of the solution process. If a negative value is specified, a default is chosen.
     *
     * @return True if a solution was found, else false. (The same value is returned by successive calls to hasSolution().)
     */
    bool solve(double tolerance = -1);

    /** Was the curve successfully fitted by the last call to solve()? */
    bool hasSolution() const;

    /**
     * Get the parameters of the logistic curve. The returned vector has <em>a</em> as its first element and <b>b</b> as the
     * subsequent elements. Valid only if hasSolution() returns true.
     */
    TheaArray<double> const & getSolution() const;

  private:
    LinearLeastSquares * llsq;  ///< The helper class for solving the linearized problem.
    TheaArray<double> llsq_coeffs;  ///< Scratch space for passing linearized values.

}; // class LogisticRegression

} // namespace Algorithms
} // namespace Thea

#endif
