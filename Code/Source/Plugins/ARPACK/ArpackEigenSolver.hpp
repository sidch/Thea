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

#ifndef __Thea_Algorithms_ArpackEigenSolver_hpp__
#define __Thea_Algorithms_ArpackEigenSolver_hpp__

#include "ArpackCommon.hpp"
#include "../../ICompressedSparseMatrix.hpp"
#include "../../IDenseMatrix.hpp"
#include "../../MatVec.hpp"
#include "../../NamedObject.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/IEigenSolver.hpp"

namespace Thea {
namespace Algorithms {

/** ARPACK-based eigensystem solver. */
class THEA_ARPACK_DLL_LOCAL ArpackEigenSolver : public IEigenSolver, public virtual NamedObject
{
  private:
    typedef IEigenSolver BaseType;

  public:
    /** Constructor. */
    ArpackEigenSolver(std::string const & name_);

    /** Destructor. */
    ~ArpackEigenSolver();

    /**
     * @copydoc IEigenSolver::solve()
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
    int64 THEA_ICALL solve(IMatrix<float64> const * m, int8 compute_eigenvectors = true, int64 num_requested_eigenpairs = -1,
                           IOptions const * options = nullptr);

    int64 THEA_ICALL dims() const { return ndims; }
    int64 THEA_ICALL numEigenpairs() const { return (int64)eigenvalues[0].size(); }
    int8 THEA_ICALL getEigenvalue(int64 i, float64 * re, float64 * im) const;
    int8 THEA_ICALL getEigenvector(int64 i, float64 const ** re, float64 const ** im) const;
    int8 THEA_ICALL hasRelativeErrors() const { return false; }
    int8 THEA_ICALL getRelativeError(int64 i, float64 * error) const;

  private:
    /** Solve a dense system */
    int64 solveDense(IDenseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma, char * which,
                     int32 ncv, float64 tol, int32 maxit, float64 * resid = 0, int8 auto_shift = true);

    /** Solve a sparse system */
    int64 solveSparse(ICompressedSparseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma,
                      char * which, int32 ncv, float64 tol, int32 maxit, float64 * resid = 0, int8 auto_shift = true);

    int64 ndims;  /**< The dimensionality of the problem, i.e. the length of each eigenvector or equivalently the size of the
                      input matrix. */
    Array<float64> eigenvalues[2];  ///< Real and imaginary parts of eigenvalues.
    Array< VectorX<float64> > eigenvectors[2];  ///< Real and imaginary parts of eigenvectors, if computed (else empty).

}; // class ArpackEigenSolver

/** Factory for creating ARPACK eigensolvers. */
class THEA_ARPACK_DLL_LOCAL ArpackEigenSolverFactory : public IEigenSolverFactory
{
  public:
    /** Destructor. */
    ~ArpackEigenSolverFactory();

    IEigenSolver * THEA_ICALL createEigenSolver(char const * name);
    int8 THEA_ICALL destroyEigenSolver(IEigenSolver * eigen_solver);

    /** Destroy all eigensolvers created with this factory. */
    int8 destroyAllEigenSolvers();

  private:
    typedef Set<IEigenSolver *> EigenSolverSet;  ///< Set of eigensolvers.

    EigenSolverSet eigen_solvers;  ///< All eigensolvers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
