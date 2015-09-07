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

#include "CSPARSELinearSolver.hpp"
#include "../../Algorithms/FastCopy.hpp"
#include "../../Array.hpp"
#include "csparse.h"
#include <cmath>

namespace Thea {
namespace Algorithms {

CSPARSELinearSolver::CSPARSELinearSolver(std::string const & name_)
: NamedObject(name_)
{}

namespace CSPARSELinearSolverInternal {

// Return value: 1: symmetric, 2: upper-triangular, 3: lower-triangular, 0: other
int
isSymmetric(cs * A)
{
  if (A->m != A->n)
    return 0;

  if (A->n <= 0 || A->nzmax <= 0 || A->p[A->n] <= 0)  // empty matrix is symmetric by default
    return 1;

  // Checking for upper or lower triangular is fast
  bool is_upper = true;
  bool is_lower = true;
  for (int j = 0 ; j < A->n; ++j)
  {
    for (int p = A->p[j]; p < A->p[j + 1]; ++p)
    {
      if (A->i[p] > j) is_upper = false;
      if (A->i[p] < j) is_lower = false;
    }
  }

  if (is_upper) return 2;
  if (is_lower) return 3;

  // Checking for symmetric is relatively slower since we compute the transpose first
  static double const EPSILON = 1.0e-200;
  cs * At = cs_transpose(A, 1);
  bool is_sym = true;

  // First check that all the column pointers are the same
  for (int j = 0 ; j < A->n; ++j)
    if (A->p[j] != At->p[j])
    {
      is_sym = false;
      break;
    }

  // Now check that all the row indices are the same
  if (is_sym)
  {
    int nnz = A->p[A->n];
    for (int p = 0; p < nnz; ++p)
      if (A->i[p] != At->i[p])
      {
        is_sym = false;
        break;
      }
  }

  // Finally check that all the values are the same (floating point comparisons so most expensive)
  if (is_sym)
  {
    int nnz = A->p[A->n];
    for (int p = 0; p < nnz; ++p)
      if (std::fabs(A->x[p] - At->x[p]) > EPSILON)
      {
        is_sym = false;
        break;
      }
  }

  // Cleanup
  cs_spfree(At);

  return is_sym ? 1 : 0;
}

// Returns true for off-diagonal entries.
int
dropDiag(int i, int j, double aij, void * other)
{
  return (i != j);
}

cs *
makeSymmetric(cs * A)
{
  cs * At, * C;
  At = cs_transpose(A, 1);
  cs_fkeep(At, &dropDiag, NULL);  // drop diagonal entries from At
  C = cs_add(A, At, 1, 1);        // C = A + At
  cs_spfree(At) ;
  return C;
}

} // namespace CSPARSELinearSolverInternal

bool
CSPARSELinearSolver::solve(Options const & options)
{
  using namespace CSPARSELinearSolverInternal;

  has_solution = false;
  solution.clear();

  switch (coeffs.getFormat())
  {
    case MatrixFormat::DENSE_ROW_MAJOR:
    case MatrixFormat::DENSE_COLUMN_MAJOR:
      throw Error(std::string(getName()) + ": The CSPARSE linear solver does not support dense coefficient matrices");

    case MatrixFormat::SPARSE_ROW_MAJOR:  // should never be encountered since we always cache as sparse column-major
      throw Error(std::string(getName()) + ": The CSPARSE linear solver does not support sparse row-major matrices");

    case MatrixFormat::SPARSE_COLUMN_MAJOR:
    {
      MatrixWrapper<double>::SparseColumnMatrix const & scm = coeffs.getSparseColumnMatrix();
      alwaysAssertM(scm.numRows() == (int)constants.size(),
                    std::string(getName()) + ": Size of coefficient matrix does not match number of constants");

      if (scm.isEmpty())
      {
        THEA_WARNING << getName() << ": Attempting to solve empty system";
        has_solution = true;
        return has_solution;
      }

      std::string method = toUpper(options.get<std::string>("method", "LU"));

      // Only QR factorization allows rectangular matrices, returning least-squares solutions
      alwaysAssertM(method == "QR" || MatrixUtil::isSquare(scm),
                    std::string(getName()) + ": Coefficient matrix for method '" + method + "' must be square");

      // Create the coefficient matrix
      // (The row indices are ints by default, but just in case the parent class changes the default format let's cache them in
      // an explicit int array.)
      TheaArray<int> pcol(scm.getColumnIndices().begin(), scm.getColumnIndices().end());
      TheaArray<int> irow(scm.getRowIndices().begin(), scm.getRowIndices().end());
      cs A;
      A.nzmax = (int)scm.numSetElements();
      A.m = scm.numRows();
      A.n = scm.numColumns();
      A.p = &pcol[0];
      A.i = &irow[0];
      A.x = const_cast<double *>(&scm.getValues()[0]);  // assume the values are always double-precision
      A.nz = -1;

      // Check if the matrix is symmetric, and if necessary symmetrize it if is triangular
      cs * Ap = &A;
      int sym = isSymmetric(&A);
      if (sym == 2 || sym == 3)
      {
        if (options.get<bool>("symmetrize-triangular", false))
          Ap = makeSymmetric(&A);
        else
          sym = 0;
      }

      // Initialize the solutions vector with the constants vector. This will be overwritten by the solver.
      solution.resize(std::max(constants.size(), (array_size_t)scm.numColumns()));
      Algorithms::fastCopy(constants.begin(), constants.end(), solution.begin());
      if (solution.size() > constants.size())
        std::fill(solution.begin() + constants.size(), solution.end(), 0);

      int ok = 0;
      try
      {
        if (method == "CHOLESKY")
        {
          THEA_DEBUG << getName() << ": Trying to solve linear system by Cholesky factorization";

          if (!sym)
            throw Error(std::string(getName())
                      + ": Cholesky factorization requires a symmetric positive-definite coefficient matrix");

          ok = cs_cholsol(Ap, &solution[0], 1);

          THEA_DEBUG << getName() << ": Solution by Cholesky factorization failed:"
                                     " please check that the coefficient matrix is symmetric positive-definite";
        }
        else if (method == "QR")
        {
          THEA_DEBUG << getName() << ": Trying to solve linear system by QR factorization";

          ok = cs_qrsol(Ap, &solution[0], 3);
        }
        else if (method == "LU")
        {
          THEA_DEBUG << getName() << ": Trying to solve linear system by LU factorization";

          double tol = options.get<double>("tol", 1e-10);
          ok = cs_lusol(Ap, &solution[0], 1, tol);
        }
        else
          throw Error(std::string(getName()) + ": Unknown solution method '" + method + '\'');
      }
      catch (...)
      {
        // Cleanup
        if (Ap != &A)
          cs_spfree(Ap);

        throw;
      }

      // Cleanup
      if (Ap != &A)
        cs_spfree(Ap);

      has_solution = (ok != 0);
      break;
    }

    default:
      throw Error(std::string(getName()) + ": Unknown format of cached matrix");
  }

  return has_solution;
}

MatrixFormat
CSPARSELinearSolver::getPreferredFormat(MatrixFormat input_format)
{
  switch (input_format)
  {
    case MatrixFormat::SPARSE_ROW_MAJOR:
    case MatrixFormat::SPARSE_COLUMN_MAJOR:
      return MatrixFormat::SPARSE_COLUMN_MAJOR;

    default:
      return input_format;
  }
}

CSPARSELinearSolverFactory::~CSPARSELinearSolverFactory()
{
  destroyAllLinearSolvers();
}

LinearSolver *
CSPARSELinearSolverFactory::createLinearSolver(std::string const & name)
{
  CSPARSELinearSolver * ls = new CSPARSELinearSolver(name);
  linear_solvers.insert(ls);
  return ls;
}

void
CSPARSELinearSolverFactory::destroyLinearSolver(LinearSolver * linear_solver)
{
  linear_solvers.erase(linear_solver);
  delete linear_solver;
}

void
CSPARSELinearSolverFactory::destroyAllLinearSolvers()
{
  for (LinearSolverSet::iterator li = linear_solvers.begin(); li != linear_solvers.end(); ++li)
    delete *li;

  linear_solvers.clear();
}

} // namespace Algorithms
} // namespace Thea
