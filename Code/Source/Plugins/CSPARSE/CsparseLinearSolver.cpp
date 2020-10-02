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

#include "CsparseLinearSolver.hpp"
#include "../../Algorithms/FastCopy.hpp"
#include "../../ICompressedSparseMatrix.hpp"
#include "../../Array.hpp"
#include "csparse.h"
#include <cmath>

namespace Thea {
namespace Algorithms {

namespace CsparseLinearSolverInternal {

static bool
indicesToInt(int32 type, int64 n, void const * begin, Array<int32> & out)
{
#define THEA_CSPARSE_CONVERT_AND_PUSH(numtype) \
  { \
    out.resize((size_t)n); \
    for (intx i = 0; i < n; ++i) out[(size_t)i] = ((numtype const *)begin)[i]; \
  }

  switch (type)
  {
    case NumericType::INT8   : THEA_CSPARSE_CONVERT_AND_PUSH(int8);   break;
    case NumericType::INT16  : THEA_CSPARSE_CONVERT_AND_PUSH(int16);  break;
    case NumericType::INT32  : THEA_CSPARSE_CONVERT_AND_PUSH(int32);  break;
    case NumericType::INT64  : THEA_CSPARSE_CONVERT_AND_PUSH(int64);  break;
    case NumericType::UINT8  : THEA_CSPARSE_CONVERT_AND_PUSH(uint8);  break;
    case NumericType::UINT16 : THEA_CSPARSE_CONVERT_AND_PUSH(uint16); break;
    case NumericType::UINT32 : THEA_CSPARSE_CONVERT_AND_PUSH(uint32); break;
    case NumericType::UINT64 : THEA_CSPARSE_CONVERT_AND_PUSH(uint64); break;

    default:
      THEA_ERROR << "CsparseLinearSolver: Unsupported compressed sparse matrix index type";
      return false;
  }

  return true;

#undef THEA_CSPARSE_CONVERT_AND_PUSH
}

// Return value: 1: symmetric, 2: upper-triangular, 3: lower-triangular, 0: other
int32
isSymmetric(cs * A)
{
  if (A->m != A->n)
    return 0;

  if (A->n <= 0 || A->nzmax <= 0 || A->p[A->n] <= 0)  // empty matrix is symmetric by default
    return 1;

  // Checking for upper or lower triangular is fast
  bool is_upper = true;
  bool is_lower = true;
  for (int32 j = 0 ; j < A->n; ++j)
  {
    for (int32 p = A->p[j]; p < A->p[j + 1]; ++p)
    {
      if (A->i[p] > j) is_upper = false;
      if (A->i[p] < j) is_lower = false;
    }
  }

  if (is_upper) return 2;
  if (is_lower) return 3;

  // Checking for symmetric is relatively slower since we compute the transpose first
  static float64 const EPSILON = 1.0e-200;
  cs * At = cs_transpose(A, 1);
  bool is_sym = true;

  // First check that all the column pointers are the same
  for (int32 j = 0 ; j < A->n; ++j)
    if (A->p[j] != At->p[j])
    {
      is_sym = false;
      break;
    }

  // Now check that all the row indices are the same
  if (is_sym)
  {
    int32 nnz = A->p[A->n];
    for (int32 p = 0; p < nnz; ++p)
      if (A->i[p] != At->i[p])
      {
        is_sym = false;
        break;
      }
  }

  // Finally check that all the values are the same (floating point comparisons so most expensive)
  if (is_sym)
  {
    int32 nnz = A->p[A->n];
    for (int32 p = 0; p < nnz; ++p)
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
int32
dropDiag(int32 i, int32 j, float64 aij, void * other)
{
  return (i != j);
}

cs *
makeSymmetric(cs * A)
{
  cs * At, * C;
  At = cs_transpose(A, 1);
  cs_fkeep(At, &dropDiag, nullptr);  // drop diagonal entries from At
  C = cs_add(A, At, 1, 1);        // C = A + At
  cs_spfree(At) ;
  return C;
}

} // namespace CsparseLinearSolverInternal

CsparseLinearSolver::CsparseLinearSolver(std::string const & name_)
: NamedObject(name_), has_solution(false)
{}

CsparseLinearSolver::~CsparseLinearSolver()
{}

int8
CsparseLinearSolver::solve(IMatrix<float64> const * a, float64 const * b, IOptions const * options)
{
  using namespace CsparseLinearSolverInternal;

  has_solution = false;
  solution.resize(0);

  try
  {
    alwaysAssertM(a, std::string(getName()) + ": Coefficient matrix is null");
    alwaysAssertM(b, std::string(getName()) + ": Constant matrix is null");

    if (a->asSparse() && a->asSparse()->asCompressed())
    {
      if (a->rows() <= 0 || a->cols() <= 0)
      {
        THEA_WARNING << getName() << ": Attempting to solve empty system";
        has_solution = true;
        return has_solution;
      }

      std::string method = toUpper(options ? options->getString("method", "LU") : "LU");

      // Only QR factorization allows rectangular matrices, returning least-squares solutions
      alwaysAssertM(method == "QR" || Math::isSquare(*a),
                    std::string(getName()) + ": Coefficient matrix for method '" + method + "' must be square");

      ICompressedSparseMatrix<float64> const & sm = *a->asSparse()->asCompressed();

      // Convert the coefficient matrix to CSPARSE format
      alwaysAssertM(sm.isColumnMajor(),
                    std::string(getName()) + ": Coefficient matrix is not in compressed column (CSC) format");
      alwaysAssertM(sm.isFullyCompressed(), std::string(getName()) + ": Operator matrix is not fully compressed");

      // TODO: avoid the copy when the input indices are actually ints
      Array<int32> irow;
      indicesToInt(sm.getInnerIndexType(), sm.numStoredElements(), sm.getInnerIndices(), irow);

      Array<int32> pcol;
      indicesToInt(sm.getOuterIndexType(), sm.outerSize(), sm.getOuterIndices(), pcol);

      int32 nnz = (int32)sm.numStoredElements();
      alwaysAssertM(nnz == pcol[pcol.size() - 1],
                    std::string(getName()) + ": (n + 1)th entry of pcol array should be number of non-zeros");

      // Create the coefficient matrix
      cs C;
      C.nzmax = nnz;
      C.m = sm.rows();
      C.n = sm.cols();
      C.p = &pcol[0];
      C.i = irow.empty() ? nullptr : &irow[0];
      C.x = const_cast<float64 *>(sm.getValues());
      C.nz = -1;

      // Check if the matrix is symmetric, and if necessary symmetrize it if is triangular
      cs * Cp = &C;
      int32 sym = isSymmetric(&C);
      if (sym == 2 || sym == 3)
      {
        if (options && options->getInteger("symmetrize-triangular", 0) != 0)
          Cp = makeSymmetric(&C);
        else
          sym = 0;
      }

      // Initialize the solutions vector with the constants vector. This will be overwritten by the solver.
      solution.resize(std::max(C.m, C.n));
      Algorithms::fastCopy(b, b + C.m, solution.data());
      if (solution.size() > C.m)  // underdetermined system
        solution.tail(solution.size() - C.m).setZero();

      int32 ok = 0;
      try
      {
        if (method == "CHOLESKY")
        {
          THEA_DEBUG << getName() << ": Trying to solve linear system by Cholesky factorization";

          if (!sym)
            throw Error("Cholesky factorization requires a symmetric positive-definite coefficient matrix");

          ok = cs_cholsol(Cp, solution.data(), 1);

          THEA_DEBUG << getName() << ": Solution by Cholesky factorization failed:"
                                     " please check that the coefficient matrix is symmetric positive-definite";
        }
        else if (method == "QR")
        {
          THEA_DEBUG << getName() << "Trying to solve linear system by QR factorization";

          ok = cs_qrsol(Cp, solution.data(), 3);
        }
        else if (method == "LU")
        {
          THEA_DEBUG << getName() << ": Trying to solve linear system by LU factorization";

          static float64 const DEFAULT_TOLERANCE = 1e-10;
          float64 tol = options ? options->getFloat("tol", DEFAULT_TOLERANCE) : DEFAULT_TOLERANCE;
          ok = cs_lusol(Cp, solution.data(), 1, tol);
        }
        else
          throw Error("Unknown solution method '" + method + '\'');
      }
      catch (...)
      {
        // Cleanup
        if (Cp != &C)
          cs_spfree(Cp);

        throw;
      }

      // Cleanup
      if (Cp != &C)
        cs_spfree(Cp);

      has_solution = (ok != 0);
    }
    else
       throw Error("Unsupported coefficient matrix format: only compressed column (CSC) allowed)");
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s: Error solving linear system", getName())

  return has_solution;
}

CsparseLinearSolverFactory::~CsparseLinearSolverFactory()
{
  destroyAllLinearSolvers();
}

ILinearSolver *
CsparseLinearSolverFactory::createLinearSolver(char const * name)
{
  CsparseLinearSolver * ls = new CsparseLinearSolver(name);
  linear_solvers.insert(ls);
  return ls;
}

int8
CsparseLinearSolverFactory::destroyLinearSolver(ILinearSolver * linear_solver)
{
  linear_solvers.erase(linear_solver);
  delete linear_solver;

  return true;
}

int8
CsparseLinearSolverFactory::destroyAllLinearSolvers()
{
  for (LinearSolverSet::iterator li = linear_solvers.begin(); li != linear_solvers.end(); ++li)
    delete *li;

  linear_solvers.clear();

  return true;
}

} // namespace Algorithms
} // namespace Thea
