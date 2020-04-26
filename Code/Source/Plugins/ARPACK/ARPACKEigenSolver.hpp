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

#ifndef __Thea_Algorithms_ARPACKEigenSolver_hpp__
#define __Thea_Algorithms_ARPACKEigenSolver_hpp__

#include "ARPACKCommon.hpp"
#include "../../AbstractCompressedSparseMatrix.hpp"
#include "../../AbstractDenseMatrix.hpp"
#include "../../MatVec.hpp"
#include "../../NamedObject.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/EigenSolver.hpp"

namespace Thea {
namespace Algorithms {

/** ARPACK-based eigensystem solver. */
class THEA_ARPACK_DLL_LOCAL ARPACKEigenSolver : public EigenSolver, public virtual NamedObject
{
  private:
    typedef EigenSolver BaseType;

  public:
    /** Constructor. */
    ARPACKEigenSolver(std::string const & name_);

    /** Destructor. */
    ~ARPACKEigenSolver();

    /**
     * @copydoc EigenSolver::solve()
     *
     * Valid options for the ARPACK backend, where n is the size of the operator matrix, are:
     * - <b>which</b>: Part of eigenvalue spectrum to be returned
     *   - <i>Type:</i> string in {"LM", "SM", "LR", "SR", "LI", "SI"}
     *   - <i>Default</i>: "LM"
     * - <b>ncv</b>: Number of Arnoldi vectors generated at each iteration (must be between \a num_requested_eigenpairs + 1 and
     *   n - 1)
     *   - <i>Type:</i> integer
     *   - <i>Default</i>: min{2 * \a num_requested_eigenpairs + 1, n - 1}
     * - <b>maxit</b>: Maximum number of Arnoldi update iterations allowed
     *   - <i>Type:</i> integer
     *   - <i>Default</i>: 100 * \a num_requested_eigenpairs
     * - <b>tol</b>: Stopping criterion, satisfied for a computed eigenvalue x if |x - x'| <= \a tol * |x|, where x' is the
     *   actual eigenvalue closest to x
     *   - <i>Type:</i> floating-point
     *   - <i>Default</i>: machine precision
     * - <b>sigma</b>: See shift-invert setting
     *   - <i>Type:</i> floating-point
     *   - <i>Default</i>: 0.0
     */
    int64 solve(AbstractMatrix<float64> const & m, int8 compute_eigenvectors = true, int64 num_requested_eigenpairs = -1,
                AbstractOptions const * options = nullptr);

    int64 dims() const { return ndims; }
    int64 numEigenpairs() const { return (int64)eigenvalues[0].size(); }
    int8 getEigenvalue(int64 i, float64 & re, float64 & im) const;
    int8 getEigenvector(int64 i, float64 const * & re, float64 const * & im) const;
    int8 hasRelativeErrors() const { return false; }
    int8 getRelativeError(int64 i, float64 & error) const;

  private:
    /** Solve a dense system */
    int64 solveDense(AbstractDenseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma, char * which,
                     int32 ncv, float64 tol, int32 maxit, float64 * resid = 0, int8 auto_shift = true);

    /** Solve a sparse system */
    int64 solveSparse(AbstractCompressedSparseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma,
                      char * which, int32 ncv, float64 tol, int32 maxit, float64 * resid = 0, int8 auto_shift = true);

    int64 ndims;  /**< The dimensionality of the problem, i.e. the length of each eigenvector or equivalently the size of the
                      input matrix. */
    Array<float64> eigenvalues[2];  ///< Real and imaginary parts of eigenvalues.
    Array< VectorX<float64> > eigenvectors[2];  ///< Real and imaginary parts of eigenvectors, if computed (else empty).

}; // class ARPACKEigenSolver

/** Factory for creating ARPACK eigensolvers. */
class THEA_ARPACK_DLL_LOCAL ARPACKEigenSolverFactory : public EigenSolverFactory
{
  public:
    /** Destructor. */
    ~ARPACKEigenSolverFactory();

    EigenSolver * createEigenSolver(char const * name);
    void destroyEigenSolver(EigenSolver * eigen_solver);

    /** Destroy all eigensolvers created with this factory. */
    void destroyAllEigenSolvers();

  private:
    typedef Set<EigenSolver *> EigenSolverSet;  ///< Set of eigensolvers.

    EigenSolverSet eigen_solvers;  ///< All eigensolvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
