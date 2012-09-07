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

#ifndef __Thea_Algorithms_ARPACKEigenSolver_hpp__
#define __Thea_Algorithms_ARPACKEigenSolver_hpp__

#include "ARPACKCommon.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/EigenSolver.hpp"

namespace Thea {
namespace Algorithms {

/** ARPACK-based eigensystem solver. */
class THEA_ARPACK_DLL_LOCAL ARPACKEigenSolver : public EigenSolver
{
  private:
    typedef EigenSolver BaseType;

  public:
    /** Constructor. */
    ARPACKEigenSolver(std::string const & name_);

    /**
     * {@inheritDoc}
     *
     * Valid options for the ARPACK backend, where n is the size of the operator matrix, are:
     * - <b>which</b>: Part of eigenvalue spectrum to be returned
     *   - <i>Type:</i> <code>std::string</code> in {"LM", "SM", "LR", "SR", "LI", "SI"}
     *   - <i>Default</i>: "LM"
     * - <b>ncv</b>: Number of Arnoldi vectors generated at each iteration (must be between \a num_requested_eigenpairs + 1 and
     *   n - 1)
     *   - <i>Type:</i> <code>int</code>
     *   - <i>Default</i>: min{2 * \a num_requested_eigenpairs + 1, n - 1}
     * - <b>maxit</b>: Maximum number of Arnoldi update iterations allowed
     *   - <i>Type:</i> <code>int</code>
     *   - <i>Default</i>: 100 * \a num_requested_eigenpairs
     * - <b>tol</b>: Stopping criterion, satisfied for a computed eigenvalue x if |x - x'| <= \a tol * |x|, where x' is the
     *   actual eigenvalue closest to x
     *   - <i>Type:</i> <code>double</code>
     *   - <i>Default</i>: machine precision
     */
    long solve(int num_requested_eigenpairs = -1, Options const & options = Options());

  protected:
    MatrixFormat getPreferredFormat(MatrixFormat input_format);

  private:
    /** Solve a dense system */
    long solveDense(int nev, bool shift_invert, double sigma, char * which, int ncv, double tol, int maxit, double * resid = 0,
                    bool AutoShift = true);

    /** Solve a sparse system */
    long solveSparse(int nev, bool shift_invert, double sigma, char * which, int ncv, double tol, int maxit, double * resid = 0,
                     bool AutoShift = true);

}; // class ARPACKEigenSolver

/** Factory for creating ARPACK eigensolvers. */
class THEA_ARPACK_DLL_LOCAL ARPACKEigenSolverFactory : public EigenSolverFactory
{
  public:
    /** Destructor. */
    ~ARPACKEigenSolverFactory();

    EigenSolver * createEigenSolver(std::string const & name);
    void destroyEigenSolver(EigenSolver * eigen_solver);

    /** Destroy all eigensolvers created with this factory. */
    void destroyAllEigenSolvers();

  private:
    typedef TheaSet<EigenSolver *> EigenSolverSet;  ///< Set of eigensolvers.

    EigenSolverSet eigen_solvers;  ///< All eigensolvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
