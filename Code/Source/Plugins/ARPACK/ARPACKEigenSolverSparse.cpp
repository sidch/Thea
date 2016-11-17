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
#include "../../MatrixUtil.hpp"
#ifdef THEA_HAVE_SUPERLU
#  include <arlsnsym.h>
#endif

namespace Thea {
namespace Algorithms {

long
ARPACKEigenSolver::solveSparse(int nev, bool shift_invert, double sigma, char * which, int ncv, double tol, int maxit,
                               double * resid, bool AutoShift)
{
#ifdef THEA_HAVE_SUPERLU

  MatrixWrapper<double>::SparseColumnMatrix const & scm = matrix.getSparseColumnMatrix();
  if (scm.isEmpty())
  {
    THEA_WARNING << getName() << ": Attempting to compute eigenvalues of an empty matrix -- no eigenpairs computed";
    return 0;
  }

  // Create the matrix
  alwaysAssertM(MatrixUtil::isSquare(scm), std::string(getName()) + ": Operator matrix is not square");
  alwaysAssertM(scm.isValid(), std::string(getName()) + ": Operator matrix has invalid internal state");

  TheaArray<int> irow(scm.getRowIndices().begin(),    scm.getRowIndices().end());
  TheaArray<int> pcol(scm.getColumnIndices().begin(), scm.getColumnIndices().end());
  int nnz = (int)scm.numSetElements();

  // Repeat the isValid() checks to check if something got screwed up during the conversion
  alwaysAssertM(irow.size() == scm.getValues().size(),
                std::string(getName()) + ": irow and nzval arrays should have same size");
  alwaysAssertM(pcol.size() == (array_size_t)scm.numRows() + 1,
                getName() + format(": pcol array should have %ld + 1 = %ld entries, instead has %ld entries",
                (long)scm.numRows(), (long)scm.numRows() + 1, (long)pcol.size()));
  alwaysAssertM(nnz == pcol[pcol.size() - 1],
                std::string(getName()) + ": (n + 1)th entry of pcol array should be number of non-zeros");

  ARluNonSymMatrix<double, double> arm(scm.numRows(), nnz, const_cast<double *>(&scm.getValues()[0]), &irow[0], &pcol[0]);

  // Setup the problem
  shared_ptr< ARluNonSymStdEig<double> > eig =
      shift_invert ? shared_ptr< ARluNonSymStdEig<double> >(new ARluNonSymStdEig<double>(nev, arm, sigma, which, ncv, tol,
                                                                                         maxit, resid, AutoShift))
                   : shared_ptr< ARluNonSymStdEig<double> >(new ARluNonSymStdEig<double>(nev, arm, which, ncv, tol, maxit,
                                                                                         resid, AutoShift));
  eig->Trace();

  // Find eigenpairs
  array_size_t nconv = (array_size_t)eig->FindEigenvectors();

  eigenvalues.resize(nconv);
  eigenvectors.resize(nconv);

  array_size_t n = (array_size_t)scm.numRows();
  for (array_size_t i = 0; i < nconv; ++i)
  {
    eigenvalues[i] = Eigenvalue(eig->EigenvalueReal((int)i), eig->EigenvalueImag((int)i));

    eigenvectors[i].resize(n);
    for (array_size_t j = 0; j < n; ++j)
      eigenvectors[i][j] = std::complex<double>(eig->EigenvectorReal((int)i, (int)j), eig->EigenvectorImag((int)i, (int)j));
  }

  return (long)nconv;

#else

  throw Error(std::string(getName()) + ": Sparse linear solver (SuperLU) not available");

#endif
}

} // namespace Algorithms
} // namespace Thea
