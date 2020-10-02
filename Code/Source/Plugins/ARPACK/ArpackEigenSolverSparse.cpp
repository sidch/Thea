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
#include "../../Array.hpp"

#ifdef THEA_HAVE_SUPERLU
#  include "arpackpp/arlsnsym.h"
#endif

namespace Thea {
namespace Algorithms {

namespace ArpackEigenSolverInternal {

#ifdef THEA_HAVE_SUPERLU

static int8
indicesToInt(int32 type, int64 n, void const * begin, Array<int32> & out)
{
#define THEA_ARPACK_CONVERT_AND_PUSH(numtype) \
  { \
    out.resize((size_t)n); \
    for (intx i = 0; i < n; ++i) out[(size_t)i] = ((numtype const *)begin)[i]; \
  }

  switch (type)
  {
    case NumericType::INT8   : THEA_ARPACK_CONVERT_AND_PUSH(int8);   break;
    case NumericType::INT16  : THEA_ARPACK_CONVERT_AND_PUSH(int16);  break;
    case NumericType::INT32  : THEA_ARPACK_CONVERT_AND_PUSH(int32);  break;
    case NumericType::INT64  : THEA_ARPACK_CONVERT_AND_PUSH(int64);  break;
    case NumericType::UINT8  : THEA_ARPACK_CONVERT_AND_PUSH(uint8);  break;
    case NumericType::UINT16 : THEA_ARPACK_CONVERT_AND_PUSH(uint16); break;
    case NumericType::UINT32 : THEA_ARPACK_CONVERT_AND_PUSH(uint32); break;
    case NumericType::UINT64 : THEA_ARPACK_CONVERT_AND_PUSH(uint64); break;

    default:
      THEA_ERROR << "ArpackEigenSolver: Unsupported compressed sparse matrix index type";
      return false;
  }

  return true;

#undef THEA_ARPACK_CONVERT_AND_PUSH
}

#endif // THEA_HAVE_SUPERLU

} // namespace ArpackEigenSolverInternal

int64
ArpackEigenSolver::solveSparse(ICompressedSparseMatrix<float64> const & m, int32 nev, int8 shift_invert, float64 sigma,
                               char * which, int32 ncv, float64 tol, int32 maxit, float64 * resid, int8 auto_shift)
{
#ifdef THEA_HAVE_SUPERLU

  try
  {
    // Create the matrix
    alwaysAssertM(m.isColumnMajor(), std::string(getName()) + ": Operator matrix is not in compressed column (CSC) format");
    alwaysAssertM(m.isFullyCompressed(), std::string(getName()) + ": Operator matrix is not fully compressed");

    Array<int32> irow;
    ArpackEigenSolverInternal::indicesToInt(m.getInnerIndexType(), m.numStoredElements(), m.getInnerIndices(), irow);

    Array<int32> pcol;
    ArpackEigenSolverInternal::indicesToInt(m.getOuterIndexType(), m.outerSize() + 1, m.getOuterIndices(), pcol);

    int32 nnz = (int32)m.numStoredElements();
    alwaysAssertM(nnz == pcol[pcol.size() - 1],
                  format("%s: (n + 1)th entry of pcol array should be number of non-zeros %ld, but is %ld",
                         getName(), (long)nnz, (long)pcol[pcol.size() - 1]));

    ARluNonSymMatrix<float64, float64> arm(m.rows(), nnz, const_cast<float64 *>(m.getValues()),
                                         (irow.empty() ? nullptr : &irow[0]), &pcol[0]);

    // Setup the problem
    std::shared_ptr< ARluNonSymStdEig<float64> > eig =
        shift_invert ? std::shared_ptr< ARluNonSymStdEig<float64> >(new ARluNonSymStdEig<float64>(
                                                                            nev, arm, sigma, which, ncv,
                                                                            tol, maxit, resid, auto_shift))
                     : std::shared_ptr< ARluNonSymStdEig<float64> >(new ARluNonSymStdEig<float64>(
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

#else

  THEA_ERROR << getName() << ": Sparse linear solver (SuperLU) not available";
  return -1;

#endif
}

} // namespace Algorithms
} // namespace Thea
