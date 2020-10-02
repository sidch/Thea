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

#ifndef __Thea_Algorithms_CsparseLinearSolver_hpp__
#define __Thea_Algorithms_CsparseLinearSolver_hpp__

#include "CsparseCommon.hpp"
#include "../../MatVec.hpp"
#include "../../NamedObject.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/ILinearSolver.hpp"

namespace Thea {
namespace Algorithms {

/**
 * CSPARSE-based solver for sparse systems of linear equations.
 *
 * @see http://www.cise.ufl.edu/research/sparse/CSparse/
 */
class THEA_CSPARSE_DLL_LOCAL CsparseLinearSolver : public ILinearSolver, public virtual NamedObject
{
  private:
    typedef ILinearSolver BaseType;

  public:
    /** Constructor. */
    CsparseLinearSolver(std::string const & name_);

    /** Destructor. */
    ~CsparseLinearSolver();

    /**
     * @copydoc ILinearSolver::solve()
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
    int8 THEA_ICALL solve(IMatrix<float64> const * a, float64 const * b, IOptions const * options = nullptr);

    int64 THEA_ICALL dims() const { return (int64)solution.size(); }
    int8 THEA_ICALL hasSolution() const { return has_solution; }
    float64 const * THEA_ICALL getSolution() const { return solution.data(); }
    int8 THEA_ICALL getSquaredError(float64 * err) const { return false; }

  protected:
    bool has_solution;
    VectorX<float64> solution;  ///< Solution vector (exact or approximate) x of Ax = b.

}; // class CsparseLinearSolver

/** Factory for creating CSPARSE linear solvers. */
class THEA_CSPARSE_DLL_LOCAL CsparseLinearSolverFactory : public ILinearSolverFactory
{
  public:
    /** Destructor. */
    ~CsparseLinearSolverFactory();

    ILinearSolver * THEA_ICALL createLinearSolver(char const * name);
    int8 THEA_ICALL destroyLinearSolver(ILinearSolver * linear_solver);

    /** Destroy all linear solvers created with this factory. */
    int8 destroyAllLinearSolvers();

  private:
    typedef Set<ILinearSolver *> LinearSolverSet;  ///< Set of linear solvers.

    LinearSolverSet linear_solvers;  ///< All linear solvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
