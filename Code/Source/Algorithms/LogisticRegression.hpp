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
// First version: 2012
//
//============================================================================

#ifndef __Thea_Algorithms_LogisticRegression_hpp__
#define __Thea_Algorithms_LogisticRegression_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class StdLinearSolver;

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
    THEA_DECL_SMART_POINTERS(LogisticRegression)

    /** Constructor. Sets the number of dimensions of the problem domain (i.e. of the vector <b>x</b>). */
    LogisticRegression(intx ndims_);

    /** Destructor. */
    ~LogisticRegression();

    /** Get the number of dimensions of the problem domain. */
    intx dims() const { return ndims; }

    /** Get the number of observations. */
    intx numObservations() const { return (intx)llsq_consts.size(); }

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

    /**
     * Fit the logistic curve to the observed data.
     *
     * @param tolerance The numerical tolerance of the solution process. If a negative value is specified, a default is chosen.
     *
     * @return True if a solution was found, else false. (The same value is returned by successive calls to hasSolution().)
     */
    bool solve(double tolerance = -1);

    /** Was the curve successfully fitted by the last call to solve()? */
    bool hasSolution() const { return has_solution; }

    /**
     * Get the parameters of the logistic curve. The returned vector has <em>a</em> as its first element and <b>b</b> as the
     * subsequent elements. Valid only if hasSolution() returns true.
     */
    VectorXd const & getSolution() const { return solution; }

  private:
    intx ndims;                     ///< The number of dimensions of the problem, i.e. the size of the solution vector x.
    Array<double> llsq_coeffs;  ///< Scratch space for passing coefficients. Has (dims() + 1) * llsq_consts.size() entries.
    Array<double> llsq_consts;  ///< Scratch space for passing constants.
    bool has_solution;              ///< Was the logistic regression problem successfully solved by solve()?
    VectorX<double> solution;       ///< The solution vector of the logistic regression problem.

}; // class LogisticRegression

} // namespace Algorithms
} // namespace Thea

#endif
