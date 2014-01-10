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
#include <algorithm>
#include <cstring>

namespace Thea {
namespace Algorithms {

ARPACKEigenSolver::ARPACKEigenSolver(std::string const & name_)
: NamedObject(name_)
{}

long
ARPACKEigenSolver::solve(int num_requested_eigenpairs, Options const & options)
{
  if (num_requested_eigenpairs < 0)
    num_requested_eigenpairs = matrix.numRows();

  std::string which_s = options.get<std::string>("which", "LM");
  char which[3];
  std::strncpy(which, which_s.c_str(), 2);
  which[2] = 0;

  int ncv     =  (int)options.get<long>  ("ncv",   0);
  double tol  =       options.get<double>("tol",   std::max(1e-7, 5e-9 * matrix.numRows()));
  int maxit   =  (int)options.get<long>  ("maxit", 100);

  double sigma = 0;
  bool shift_invert = options.hasOption("sigma");
  if (shift_invert)
    sigma = options.get<double>("sigma", 0);

  THEA_DEBUG << "nev = " << num_requested_eigenpairs << ", which = " << which << ", ncv = " << ncv << ", tol = " << tol
             << ", maxit = " << maxit;

  eigenvalues.clear();
  eigenvectors.clear();

  switch (matrix.getFormat())
  {
    case MatrixFormat::DENSE_ROW_MAJOR:
    case MatrixFormat::DENSE_COLUMN_MAJOR:
      return solveDense(num_requested_eigenpairs, shift_invert, sigma, which, ncv, tol, maxit);

    case MatrixFormat::SPARSE_ROW_MAJOR:
    case MatrixFormat::SPARSE_COLUMN_MAJOR:
      return solveSparse(num_requested_eigenpairs, shift_invert, sigma, which, ncv, tol, maxit);

    default:
      throw Error(std::string(getName()) + ": Unknown format of cached matrix");
  }
}

MatrixFormat
ARPACKEigenSolver::getPreferredFormat(MatrixFormat input_format)
{
  switch (input_format)
  {
    case MatrixFormat::SPARSE_ROW_MAJOR:
    case MatrixFormat::SPARSE_COLUMN_MAJOR:
      return MatrixFormat::SPARSE_COLUMN_MAJOR;

    default:
      return MatrixFormat::DENSE_COLUMN_MAJOR;
  }
}

ARPACKEigenSolverFactory::~ARPACKEigenSolverFactory()
{
  destroyAllEigenSolvers();
}

EigenSolver *
ARPACKEigenSolverFactory::createEigenSolver(std::string const & name)
{
  ARPACKEigenSolver * es = new ARPACKEigenSolver(name);
  eigen_solvers.insert(es);
  return es;
}

void
ARPACKEigenSolverFactory::destroyEigenSolver(EigenSolver * eigen_solver)
{
  eigen_solvers.erase(eigen_solver);
  delete eigen_solver;
}

void
ARPACKEigenSolverFactory::destroyAllEigenSolvers()
{
  for (EigenSolverSet::iterator ei = eigen_solvers.begin(); ei != eigen_solvers.end(); ++ei)
    delete *ei;

  eigen_solvers.clear();
}

} // namespace Algorithms
} // namespace Thea
