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

#ifndef __Thea_Algorithms_CSPARSELinearSolver_hpp__
#define __Thea_Algorithms_CSPARSELinearSolver_hpp__

#include "CSPARSECommon.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/LinearSolver.hpp"

namespace Thea {
namespace Algorithms {

/**
 * CSPARSE-based solver for sparse systems of linear equations.
 *
 * @see http://www.cise.ufl.edu/research/sparse/CSparse/
 */
class THEA_CSPARSE_DLL_LOCAL CSPARSELinearSolver : public LinearSolver
{
  private:
    typedef LinearSolver BaseType;

  public:
    /** Constructor. */
    CSPARSELinearSolver(std::string const & name_);

    /**
     * {@inheritDoc}
     *
     * Valid options for the CSPARSE backend are:
     * - <b>method</b>: Solution method to use
     *   - <i>Type:</i> <code>std::string</code> in {"LU", "QR", "Cholesky"}
     *   - <i>Default</i>: "LU"
     *   - <i>Note 1</i>: Cholesky factorization requires that the coefficient matrix be symmetric positive-definite.
     *   - <i>Note 2</i>: QR factorization allows rectangular coefficient matrices, in which case a least-squares solution is
     *     returned. The other methods all require square matrices.
     * - <b>symmetrize-triangular</b>: If the input matrix is upper- or lower-triangular, it is symmetrized by copying the
     *   existing triangle to its symmetric counterpart
     *   - <i>Type:</i> <code>bool</code>
     *   - <i>Default</i>: <code>false</code>
     * - <b>tol</b>: Tolerance for LU factorization
     *   - <i>Type:</i> <code>double</code>
     *   - <i>Default</i>: 1e-20
     */
    bool solve(Options const & options = Options());

  protected:
    MatrixFormat getPreferredFormat(MatrixFormat input_format);

}; // class CSPARSELinearSolver

/** Factory for creating CSPARSE linear solvers. */
class THEA_CSPARSE_DLL_LOCAL CSPARSELinearSolverFactory : public LinearSolverFactory
{
  public:
    /** Destructor. */
    ~CSPARSELinearSolverFactory();

    LinearSolver * createLinearSolver(std::string const & name);
    void destroyLinearSolver(LinearSolver * linear_solver);

    /** Destroy all linear solvers created with this factory. */
    void destroyAllLinearSolvers();

  private:
    typedef TheaSet<LinearSolver *> LinearSolverSet;  ///< Set of linear solvers.

    LinearSolverSet linear_solvers;  ///< All linear solvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
