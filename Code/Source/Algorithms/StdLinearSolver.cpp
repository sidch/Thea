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

#include "StdLinearSolver.hpp"
#include "../AbstractDenseMatrix.hpp"
#include "../AbstractCompressedSparseMatrix.hpp"
#include "NNLS/nnls.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/QR>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SVD>
#include <type_traits>

namespace Thea {
namespace Algorithms {

namespace StdLinearSolverInternal {

// Implementation of StdLinearSolver functions.
class THEA_DLL_LOCAL StdLinearSolverImpl
{
  public:
    // Constructor.
    StdLinearSolverImpl(StdLinearSolver::Method method_, StdLinearSolver::Constraint constraint_)
    : method(method_), constraint(constraint_), tolerance(-1), max_iters(-1), ndims(0), has_solution(false)
    {}

    // Solve the linear system Ax = b for a dense double-precision matrix A.
    template <typename MatrixT>
    bool solve(Eigen::MatrixBase<MatrixT> const & a, double const * b, AbstractOptions const * options = NULL,
               std::enable_if< std::is_same<typename MatrixT::value_type, double>::value > * dummy = NULL)
    {
      if (a.rows() < a.cols())
        THEA_WARNING << "StdLinearSolver: Fewer objectives than dimensions -- the solution will not be unique";

      has_solution = false;
      try
      {
        if (a.rows() <= 0 || a.cols() <= 0)
          throw Error("Empty coefficient matrix");

        long num_objectives = a.rows();
        ndims = a.cols();

        switch (constraint)
        {
          case StdLinearSolver::Constraint::NON_NEGATIVE:
          {
            if (method != StdLinearSolver::Method::DEFAULT && method != StdLinearSolver::Method::NNLS)
              throw Error("Unsupported method for non-negative least squares problems");

            MatrixX<double, MatrixLayout::COLUMN_MAJOR> nnls_a = a;  // NNLS requires a COLUMN-MAJOR (Fortran-style) matrix.
                                                                     // Plus values will be overwritten, so need a copy anyway.
            VectorXd nnls_b = VectorXdConstMap(b, num_objectives);  // will also be overwritten, so make a copy
            solution.resize(ndims);

            double         rnorm;
            Array<double>  w((size_t)ndims);
            Array<double>  zz((size_t)num_objectives);
            Array<int>     index((size_t)ndims);
            int            mode;

            int mda = (int)num_objectives, im = (int)num_objectives, in = (int)ndims;
            nnls_c(nnls_a.data(), &mda, &im, &in, nnls_b.data(), solution.data(), &rnorm, &w[0], &zz[0], &index[0], &mode);

            if (mode == 1)
              has_solution = true;
            else
            {
              switch (mode)
              {
                // Should never be 2 since we've checked for this above, but recheck all the same
                case 2:  THEA_DEBUG << "StdLinearSolver: NNLS error (bad problem dimensions)"; break;
                case 3:  THEA_DEBUG << "StdLinearSolver: NNLS error (iteration count exceeded)"; break;
                default: THEA_DEBUG << "StdLinearSolver: Unknown NNLS error";
              }
            }

            break;
          }

          case StdLinearSolver::Constraint::UNCONSTRAINED:
          {
            switch (method)
            {
              case StdLinearSolver::Method::HOUSEHOLDER_QR:
              {
                Eigen::HouseholderQR<typename MatrixT::PlainObject> solver(a);
                solution = solver.solve(VectorXdConstMap(b, num_objectives));
                has_solution = true;
                break;
              }

              case StdLinearSolver::Method::DEFAULT:  // slower than plain Householder, but more accurate
              case StdLinearSolver::Method::COL_PIV_HOUSEHOLDER_QR:
              {
                Eigen::ColPivHouseholderQR<typename MatrixT::PlainObject> solver(a);
                if (tolerance >= 0) solver.setThreshold(tolerance);
                solution = solver.solve(VectorXdConstMap(b, num_objectives));
                has_solution = true;
                break;
              }

              case StdLinearSolver::Method::FULL_PIV_HOUSEHOLDER_QR:
              {
                Eigen::FullPivHouseholderQR<typename MatrixT::PlainObject> solver(a);
                if (tolerance >= 0) solver.setThreshold(tolerance);
                solution = solver.solve(VectorXdConstMap(b, num_objectives));
                has_solution = true;
                break;
              }

              case StdLinearSolver::Method::COMPLETE_ORTHOGONAL_DECOMPOSITION:
              {
                Eigen::CompleteOrthogonalDecomposition<typename MatrixT::PlainObject> solver(a);
                if (tolerance >= 0) solver.setThreshold(tolerance);
                solution = solver.solve(VectorXdConstMap(b, num_objectives));
                has_solution = true;
                break;
              }

              case StdLinearSolver::Method::BDCSVD:
              {
                Eigen::BDCSVD<typename MatrixT::PlainObject> solver(a);
                if (tolerance >= 0) solver.setThreshold(tolerance);
                solution = solver.solve(VectorXdConstMap(b, num_objectives));
                has_solution = true;
                break;
              }

              default:
                throw Error("Unsupported method for unconstrained dense least squares problems");
            }

            break;
          }

          default:
            throw Error("Unsupported constraint");
        }
      }
      THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s",
                                 "StdLinearSolver: Error solving dense linear least-squares system")

      return has_solution;
    }

    // Solve the linear system Ax = b for a sparse double-precision matrix A.
    template <typename MatrixT>
    bool solve(Eigen::SparseMatrixBase<MatrixT> const & a, double const * b, AbstractOptions const * options = NULL,
               std::enable_if< std::is_same<typename MatrixT::value_type, double>::value > * dummy = NULL)
    {
      if (a.rows() < a.cols())
        THEA_WARNING << "StdLinearSolver: Fewer objectives than dimensions -- the solution will not be unique";

      has_solution = false;
      try
      {
        if (a.rows() <= 0 || a.cols() <= 0)
          throw Error("Empty coefficient matrix");

        ndims = a.cols();

        switch (constraint)
        {
          case StdLinearSolver::Constraint::UNCONSTRAINED:
          {
            // Return true if a matching solver was found
            if (solveSparseFactorize(a, b)) break;
            if (solveIterative(a, b)) break;

            throw Error("Unsupported method for unconstrained sparse least squares problems");
          }

          default:
            throw Error("Unsupported constraint");
        }
      }
      THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s",
                                 "StdLinearSolver: Error solving sparse linear least-squares system")

      return has_solution;
    }

  private:
    // Use one of the iterative solvers to solve the dense or sparse problem. Returns true if a suitable solver was found, NOT
    // if the problem was successfully solved.
    template <typename MatrixT> bool solveIterative(MatrixT const & a, double const * b)
    {
      switch (method)
      {
        case StdLinearSolver::Method::CONJUGATE_GRADIENT:
        {
          Eigen::ConjugateGradient<typename MatrixT::PlainObject> solver;
          if (tolerance >= 0) solver.setTolerance(tolerance);
          if (max_iters > 0) solver.setMaxIterations(max_iters);
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        case StdLinearSolver::Method::DEFAULT:
        case StdLinearSolver::Method::LEAST_SQUARES_CONJUGATE_GRADIENT:
        {
          Eigen::LeastSquaresConjugateGradient<typename MatrixT::PlainObject> solver;
          if (tolerance >= 0) solver.setTolerance(tolerance);
          if (max_iters > 0) solver.setMaxIterations(max_iters);
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        case StdLinearSolver::Method::BICGSTAB:
        {
          Eigen::BiCGSTAB<typename MatrixT::PlainObject> solver;
          if (tolerance >= 0) solver.setTolerance(tolerance);
          if (max_iters > 0) solver.setMaxIterations(max_iters);
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        default: return false;
      }

      return true;
    }

    // Use a solver based on sparse factorization, for column-major matrices. Returns true if a suitable solver was found, NOT
    // if the problem was successfully solved.
    template <typename MatrixT>
    bool solveSparseFactorize(MatrixT const & a, double const * b,
                              typename std::enable_if< !(MatrixT::Flags & Eigen::RowMajorBit) >::type * dummy = NULL)
    {
      switch (method)
      {
        case StdLinearSolver::Method::SIMPLICIALT_LLT:
        {
          Eigen::SimplicialLLT<typename MatrixT::PlainObject> solver;
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        case StdLinearSolver::Method::SIMPLICIALT_LDLT:
        {
          Eigen::SimplicialLDLT<typename MatrixT::PlainObject> solver;
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        case StdLinearSolver::Method::SPARSE_LU:
        {
          Eigen::SparseLU<typename MatrixT::PlainObject> solver;
          if (tolerance >= 0) solver.setPivotThreshold(tolerance);
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        case StdLinearSolver::Method::DEFAULT:
        case StdLinearSolver::Method::SPARSE_QR:
        {
          Eigen::SparseQR<typename MatrixT::PlainObject, Eigen::AMDOrdering<typename MatrixT::StorageIndex> > solver;
          if (tolerance >= 0) solver.setPivotThreshold(tolerance);
          solver.compute(a);
          if (solver.info() == Eigen::Success)
          {
            solution = solver.solve(VectorXdMap(const_cast<double *>(b), a.rows()));
            has_solution = (solver.info() == Eigen::Success);
          }
          break;
        }

        default: return false;
      }

      return true;
    }

    // Use a solver based on sparse factorization, for row-major matrices. Eigen does not provide any such solvers, so this
    // method is empty. Returns false to indicate no solver was found.
    template <typename MatrixT>
    bool solveSparseFactorize(MatrixT const & a, double const * b,
                              typename std::enable_if< (MatrixT::Flags & Eigen::RowMajorBit) >::type * dummy = NULL)
    {
      return false;
    }

  public:  // no need for accessor fns for this internal class, and friend declarations are complicated by namespaces
    StdLinearSolver::Method method;          // Solution method.
    StdLinearSolver::Constraint constraint;  // Solution constraint.
    double tolerance;                        // Solution tolerance/threshold.
    long max_iters;                          // Maximum number of solver iterations, if solver is iterative.
    long ndims;                              // Solution dimensions.
    bool has_solution;                       // Was a solution computed by the last call to solve()? */
    VectorX<double> solution;                // The solution vector <b>x</b>.
};

// Shorthand
template <MatrixLayout::Value L, typename StorageIndex> using SM = Eigen::Map< SparseMatrix<double, L, StorageIndex> >;

// Given a storage index type, dispatch a call to StdLinearSolverImpl::solve() with the correct pointer conversions
template <MatrixLayout::Value L>
bool
implSolve(StdLinearSolverImpl * impl, int storage_type, long nr, long nc, long nnz, void const * in, void const * out,
          double const * val, void const * nzc, double const * b, AbstractOptions const * opt)
{
  void   * in2  = const_cast<void *>(in);
  void   * out2 = const_cast<void *>(out);
  double * val2 = const_cast<double *>(val);
  void   * nzc2 = const_cast<void *>(nzc);

  switch (storage_type)
  {
// A casting bug in the current Eigen prevents us using StorageIndex shorter than an int.
// Specifically, OrderingMethods/Amd.h, L104: dense = (std::min)(n-2, dense);
// Here n and dense are of type StorageIndex, but n - 2 is of type int for StorageIndex shorter than an int, so std::min fails
// because the arguments are not of the same type.
//
//     case NumericType::INT8:
//       return impl->solve(SM<L, int8  >(nr, nc, nnz, (int8   *)in2, (int8   *)out2, val2, (int8   *)nzc2), b, opt);
//     case NumericType::INT16:
//       return impl->solve(SM<L, int16 >(nr, nc, nnz, (int16  *)in2, (int16  *)out2, val2, (int16  *)nzc2), b, opt);
    case NumericType::INT32:
      return impl->solve(SM<L, int32 >(nr, nc, nnz, (int32  *)in2, (int32  *)out2, val2, (int32  *)nzc2), b, opt);
    case NumericType::INT64:
      return impl->solve(SM<L, int64 >(nr, nc, nnz, (int64  *)in2, (int64  *)out2, val2, (int64  *)nzc2), b, opt);
    default: THEA_ERROR << "StdLinearSolver: Unsupported index type";
  }

  return false;
}

} // namespace StdLinearSolverInternal

StdLinearSolver::StdLinearSolver(Method method_, Constraint constraint_)
: NamedObject("StdLinearSolver"), impl(new StdLinearSolverInternal::StdLinearSolverImpl(method_, constraint_))
{
}

StdLinearSolver::~StdLinearSolver()
{
  delete impl;
}

StdLinearSolver::Method
StdLinearSolver::getMethod() const
{
  return impl->method;
}

StdLinearSolver::Constraint
StdLinearSolver::getConstraint() const
{
  return impl->constraint;
}

double
StdLinearSolver::getTolerance() const
{
  return impl->tolerance;
}

long
StdLinearSolver::maxIterations() const
{
  return impl->max_iters;
}

void
StdLinearSolver::setMethod(StdLinearSolver::Method method_)
{
  impl->method = method_;
}

void
StdLinearSolver::setConstraint(StdLinearSolver::Constraint constraint_)
{
  impl->constraint = constraint_;
}

void
StdLinearSolver::setTolerance(double tol)
{
  impl->tolerance = tol;
}

void
StdLinearSolver::setMaxIterations(long max_iters_)
{
  impl->max_iters = max_iters_;
}

bool
StdLinearSolver::solve(Eigen::Ref<MatrixXd> const & a, double const * b, AbstractOptions const * options)
{
  return impl->solve(a, b, options);
}

bool
StdLinearSolver::solve(Eigen::Ref< SparseMatrix<double> > const & a, double const * b, AbstractOptions const * options)
{
  return impl->solve(a, b, options);
}

bool
StdLinearSolver::solve(AbstractMatrix<double> const & a, double const * b, AbstractOptions const * options)
{
  if (a.asAddressable() && a.asAddressable()->asDense())
  {
    AbstractDenseMatrix<double> const & dm = *a.asAddressable()->asDense();
    if (dm.isRowMajor())
    {
      Eigen::Map< MatrixX<double, MatrixLayout::ROW_MAJOR> const > wrapped(dm.data(), dm.rows(), dm.cols());
      return impl->solve(wrapped, b, options);
    }
    else  // col-major
    {
      Eigen::Map< MatrixX<double, MatrixLayout::COLUMN_MAJOR> const > wrapped(dm.data(), dm.rows(), dm.cols());
      return impl->solve(wrapped, b, options);
    }
  }
  else if (a.asSparse() && a.asSparse()->asCompressed())
  {
    AbstractCompressedSparseMatrix<double> const & sm = *a.asSparse()->asCompressed();
    int storage_type = sm.getInnerIndexType();
    if (storage_type != sm.getOuterIndexType() || storage_type != sm.getNonZeroCountType())
    {
      // TODO: Convert to integer arrays of consistent type to work around this problem

      THEA_ERROR << "StdLinearSolver: Different indices have different storage types -- cannot convert to SparseMatrix";
      return false;
    }

    if (sm.isRowMajor())
      return StdLinearSolverInternal::implSolve<MatrixLayout::ROW_MAJOR>(impl, storage_type, sm.rows(), sm.cols(),
                                                                         sm.numStoredElements(), sm.getOuterIndices(),
                                                                         sm.getInnerIndices(), sm.getValues(),
                                                                         sm.getNonZeroCounts(), b, options);
    else
      return StdLinearSolverInternal::implSolve<MatrixLayout::COLUMN_MAJOR>(impl, storage_type, sm.rows(), sm.cols(),
                                                                            sm.numStoredElements(), sm.getOuterIndices(),
                                                                            sm.getInnerIndices(), sm.getValues(),
                                                                            sm.getNonZeroCounts(), b, options);
  }
  else
  {
    THEA_ERROR << "StdLinearSolver: Unsupported matrix type";
    return false;
  }
}

long
StdLinearSolver::dims() const
{
  return impl->ndims;
}

bool
StdLinearSolver::hasSolution() const
{
  return impl->has_solution;
}

double const *
StdLinearSolver::getSolution() const
{
  return impl->solution.data();
}

bool
StdLinearSolver::getSquaredError(double & err) const
{
  return false;
}

} // namespace Algorithms
} // namespace Thea
