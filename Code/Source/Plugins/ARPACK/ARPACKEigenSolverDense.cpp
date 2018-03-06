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
#include "../../MatrixUtil.hpp"
#include <ardsnsym.h>

namespace Thea {
namespace Algorithms {

long
ARPACKEigenSolver::solveDense(int nev, bool shift_invert, double sigma, char * which, int ncv, double tol, int maxit,
                              double * resid, bool AutoShift)
{
  MatrixWrapper<double>::DenseColumnMatrix const & dcm = matrix.getDenseColumnMatrix();
  if (dcm.isEmpty())
  {
    THEA_WARNING << getName() << ": Attempting to compute eigenvalues of an empty matrix -- no eigenpairs computed";
    return 0;
  }

  // Create the matrix
  alwaysAssertM(MatrixUtil::isSquare(dcm), std::string(getName()) + ": Operator matrix is not square");
  ARdsNonSymMatrix<double, double> arm(dcm.numRows(), const_cast<double *>(&dcm.data()[0]));

  // Setup the problem
  shared_ptr< ARdsNonSymStdEig<double> > eig =
      shift_invert ? shared_ptr< ARdsNonSymStdEig<double> >(new ARdsNonSymStdEig<double>(nev, arm, sigma, which, ncv, tol,
                                                                                         maxit, resid, AutoShift))
                   : shared_ptr< ARdsNonSymStdEig<double> >(new ARdsNonSymStdEig<double>(nev, arm, which, ncv, tol, maxit,
                                                                                         resid, AutoShift));
  eig->Trace();

  // Find eigenpairs
  size_t nconv = (size_t)eig->FindEigenvectors();

  eigenvalues.resize(nconv);
  eigenvectors.resize(nconv);

  size_t n = (size_t)dcm.numRows();
  for (size_t i = 0; i < nconv; ++i)
  {
    eigenvalues[i] = Eigenvalue(eig->EigenvalueReal((int)i), eig->EigenvalueImag((int)i));

    eigenvectors[i].resize(n);
    for (size_t j = 0; j < n; ++j)
      eigenvectors[i][j] = std::complex<double>(eig->EigenvectorReal((int)i, (int)j), eig->EigenvectorImag((int)i, (int)j));
  }

  return (long)nconv;
}

} // namespace Algorithms
} // namespace Thea
