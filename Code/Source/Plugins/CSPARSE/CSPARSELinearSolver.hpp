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
// First version: 2009
//
//============================================================================

#ifndef __Thea_Algorithms_CSPARSELinearSolver_hpp__
#define __Thea_Algorithms_CSPARSELinearSolver_hpp__

#include "CSPARSECommon.hpp"
#include "../../MatVec.hpp"
#include "../../NamedObject.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/LinearSolver.hpp"

namespace Thea {
namespace Algorithms {

/**
 * CSPARSE-based solver for sparse systems of linear equations.
 *
 * @see http://www.cise.ufl.edu/research/sparse/CSparse/
 */
class THEA_CSPARSE_DLL_LOCAL CSPARSELinearSolver : public LinearSolver, public virtual NamedObject
{
  private:
    typedef LinearSolver BaseType;

  public:
    /** Constructor. */
    CSPARSELinearSolver(std::string const & name_);

    /** Destructor. */
    ~CSPARSELinearSolver();

    /**
     * @copydoc LinearSolver::solve()
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
     *   - <i>Type:</i> <code>int8</code>
     *   - <i>Default</i>: <code>false</code>
     * - <b>tol</b>: Tolerance for LU factorization
     *   - <i>Type:</i> <code>float64</code>
     *   - <i>Default</i>: 1e-20
     */
    int8 solve(AbstractMatrix<float64> const & a, float64 const * b, AbstractOptions const * options = nullptr);

    int64 dims() const { return (int64)solution.size(); }
    int8 hasSolution() const { return has_solution; }
    float64 const * getSolution() const { return solution.data(); }
    int8 getSquaredError(float64 & err) const { return false; }

  protected:
    bool has_solution;
    VectorX<float64> solution;  ///< Solution vector (exact or approximate) x of Ax = b.

}; // class CSPARSELinearSolver

/** Factory for creating CSPARSE linear solvers. */
class THEA_CSPARSE_DLL_LOCAL CSPARSELinearSolverFactory : public LinearSolverFactory
{
  public:
    /** Destructor. */
    ~CSPARSELinearSolverFactory();

    LinearSolver * createLinearSolver(char const * name);
    void destroyLinearSolver(LinearSolver * linear_solver);

    /** Destroy all linear solvers created with this factory. */
    void destroyAllLinearSolvers();

  private:
    typedef Set<LinearSolver *> LinearSolverSet;  ///< Set of linear solvers.

    LinearSolverSet linear_solvers;  ///< All linear solvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
