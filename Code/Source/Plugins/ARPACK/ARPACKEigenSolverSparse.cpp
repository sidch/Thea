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

#include "ARPACKEigenSolver.hpp"
#include "../../Array.hpp"

#ifdef THEA_HAVE_SUPERLU
#  include <arlsnsym.h>
#endif

namespace Thea {
namespace Algorithms {

namespace ARPACKEigenSolverInternal {

#ifdef THEA_HAVE_SUPERLU

static bool
indicesToInt(int type, long n, void const * begin, Array<int> & out)
{
#define THEA_ARPACK_CONVERT_AND_PUSH(numtype) \
  { \
    out.resize((size_t)n); \
    for (long i = 0; i < n; ++i) out[(size_t)i] = ((numtype const *)begin)[i]; \
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
      THEA_ERROR << "ARPACKEigenSolver: Unsupported compressed sparse matrix index type";
      return false;
  }

  return true;

#undef THEA_ARPACK_CONVERT_AND_PUSH
}

#endif // THEA_HAVE_SUPERLU

} // namespace ARPACKEigenSolverInternal

long
ARPACKEigenSolver::solveSparse(AbstractCompressedSparseMatrix<double> const & m, int nev, bool shift_invert, double sigma,
                               char * which, int ncv, double tol, int maxit, double * resid, bool auto_shift)
{
#ifdef THEA_HAVE_SUPERLU

  try
  {
    // Create the matrix
    alwaysAssertM(m.isColumnMajor(), std::string(getName()) + ": Operator matrix is not in compressed column (CSC) format");
    alwaysAssertM(m.isFullyCompressed(), std::string(getName()) + ": Operator matrix is not fully compressed");

    Array<int> irow;
    ARPACKEigenSolverInternal::indicesToInt(m.getInnerIndexType(), m.numStoredElements(), m.getInnerIndices(), irow);

    Array<int> pcol;
    ARPACKEigenSolverInternal::indicesToInt(m.getOuterIndexType(), m.outerSize(), m.getOuterIndices(), pcol);

    int nnz = (int)m.numStoredElements();
    alwaysAssertM(nnz == pcol[pcol.size() - 1],
                  std::string(getName()) + ": (n + 1)th entry of pcol array should be number of non-zeros");

    ARluNonSymMatrix<double, double> arm(m.rows(), nnz, const_cast<double *>(m.getValues()),
                                         (irow.empty() ? NULL : &irow[0]), &pcol[0]);

    // Setup the problem
    std::shared_ptr< ARluNonSymStdEig<double> > eig =
        shift_invert ? std::shared_ptr< ARluNonSymStdEig<double> >(new ARluNonSymStdEig<double>(nev, arm, sigma, which, ncv,
                                                                                                tol, maxit, resid, auto_shift))
                     : std::shared_ptr< ARluNonSymStdEig<double> >(new ARluNonSymStdEig<double>(nev, arm, which, ncv, tol,
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
      for (long j = 0; j < ndims; ++j)
      {
        eigenvectors[0][i][j] = eig->EigenvectorReal((int)i, (int)j);
        eigenvectors[1][i][j] = eig->EigenvectorImag((int)i, (int)j);
      }
    }

    return (long)nconv;
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s: Error solving dense eigensystem", getName())

#else

  THEA_ERROR << getName() << ": Sparse linear solver (SuperLU) not available";
  return -1;

#endif
}

} // namespace Algorithms
} // namespace Thea
