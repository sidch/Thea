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

#include "ArpackEigenSolver.hpp"
#include "arpackpp/ardsnsym.h"

namespace Thea {
namespace Algorithms {

int64
ArpackEigenSolver::solveDense(IDenseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma, char * which,
                              int32 ncv, float64 tol, int32 maxit, float64 * resid, int8 auto_shift)
{
  try
  {
    // Create the matrix
    ARdsNonSymMatrix<float64, float64> arm(m.rows(), const_cast<float64 *>(m.data()));

    // Setup the problem
    std::shared_ptr< ARdsNonSymStdEig<float64> > eig =
        shift_invert ? std::shared_ptr< ARdsNonSymStdEig<float64> >(new ARdsNonSymStdEig<float64>(
                                                                            nev, arm, sigma, which, ncv,
                                                                            tol, maxit, resid, auto_shift))
                     : std::shared_ptr< ARdsNonSymStdEig<float64> >(new ARdsNonSymStdEig<float64>(
                                                                            nev, arm, which, ncv, tol,
                                                                            maxit, resid, auto_shift));
    eig->Trace();

    // Find eigenpairs
    size_t nconv = (size_t)eig->FindEigenvectors();

    eigenvalues [0].resize(nconv); eigenvalues [1].resize(nconv);
    eigenvectors[0].resize(nconv); eigenvectors[1].resize(nconv);

    for (size_t i = 0; i < nconv; ++i)
    {
      eigenvalues[0][i] = eig->EigenvalueReal((int)i);
      eigenvalues[1][i] = eig->EigenvalueImag((int)i);

      eigenvectors[0][i].resize(ndims); eigenvectors[1][i].resize(ndims);
      for (intx j = 0; j < ndims; ++j)
      {
        eigenvectors[0][i][j] = eig->EigenvectorReal((int)i, (int)j);
        eigenvectors[1][i][j] = eig->EigenvectorImag((int)i, (int)j);
      }
    }

    return (int64)nconv;
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s: Error solving dense eigensystem", getName())
}

} // namespace Algorithms
} // namespace Thea
