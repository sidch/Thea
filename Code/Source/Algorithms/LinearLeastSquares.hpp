//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_LinearLeastSquares_hpp__
#define __Thea_Algorithms_LinearLeastSquares_hpp__

#include "../Common.hpp"
#include "../Array.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Solve linear least-squares problems of the form: minimize ||A<b>x</b> - <b>b</b>||. Here ||.|| denotes the L2 norm. This
 * class provides a few basic solvers. More advanced linear least-squares solvers may be provided by plugins implementing the
 * LinearSolver interface.
 */
class THEA_API LinearLeastSquares
{
  public:
    THEA_DEF_POINTER_TYPES(LinearLeastSquares, shared_ptr, weak_ptr)

    /** The constraints to be imposed on the solution (enum class). */
    struct THEA_API Constraint
    {
      /** Supported values. */
      enum Value
      {
        UNCONSTRAINED,  ///< The solution is not constrained.
        NON_NEGATIVE    ///< The solution must have non-negative elements only.
      };

      THEA_ENUM_CLASS_BODY(Constraint)
    };

    /** Constructor. Sets the number of dimensions of the problem domain (i.e. of the vector <b>x</b>). */
    LinearLeastSquares(long ndim_);

    /** Get the number of dimensions of the problem domain. */
    long numDimensions() const { return ndim; }

    /**
     * Add an objective term that is to be satisfied in a least-squares sense. This corresponds to a row of the matrix A, and
     * the corresponding element of <b>b</b>. The equality dot(\a coeffs, <b>x</b>) = \a b is satisfied as best as possible.
     *
     * @param coeffs The coefficients of the elements of <b>x</b>.
     * @param constant The desired value of the dot product of \a coeffs and <b>x</b>.
     *
     * @see clearObjectives()
     */
    void addObjective(double const * coeffs, double constant);

    /**
     * Clear all objectives.
     *
     * @see addObjective()
     */
    void clearObjectives();

    /** Get the number of objectives (the number of rows in the matrix A). */
    long numObjectives() const { return num_objectives; }

    /**
     * Minimize the squared error of the linear system. The elements of the solution vector may be constrained (see Constraint).
     * Currently, an unconstrained solution (\a non_neg > 0) is computed via singular value decomposition (SVD) and a
     * non-negative solution is computed via NNLS (C. L. Lawson and R. J. Hanson, "Solving Least Square Problems", Prentice
     * Hall, Englewood Cliffs NJ, 1974, http://hesperia.gsfc.nasa.gov/~schmahl/nnls).
     *
     * @param constraint The constraints to be imposed on the solution.
     * @param tolerance The numerical tolerance of the solution process. For a solution via SVD, singular values smaller than
     *   the tolerance are zeroed out.
     *
     * @return True if a solution was found, else false. (The same value is returned by successive calls to hasSolution().)
     */
    bool solve(Constraint constraint = Constraint::UNCONSTRAINED, double tolerance = -1);

    /** Was the linear system successfully solved by the last call to solve()? */
    bool hasSolution() const { return has_solution; }

    /** Get the solution vector of the linear system. Valid only if hasSolution() returns true. */
    TheaArray<double> const & getSolution() const { return solution; }

    /**
     * Compute the squared error of the current solution, or return a negative value if no solution exists or the system is
     * invalid.
     */
    double squaredError() const;

  private:
    long ndim;  ///< The number of dimensions of the problem domain.
    long num_objectives;  ///< The number of objectives in the system.
    TheaArray<double> a_values;  ///< The values, in row-major form, of the matrix A.
    TheaArray<double> b;  ///< The constants vector <b>b</b>.
    bool has_solution;  ///< Was a solution computed by the last call to solve()? */
    TheaArray<double> solution;  ///< The solution vector <b>x</b>.

}; // class LinearLeastSquares

} // namespace Algorithms
} // namespace Thea

#endif
