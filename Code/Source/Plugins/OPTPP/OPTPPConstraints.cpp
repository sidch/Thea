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

#include "OPTPPConstraints.hpp"
#include "../../Array.hpp"
#include <LinearEquation.h>
#include <LinearInequality.h>
#include <NonLinearEquation.h>
#include <NonLinearInequality.h>

namespace Thea {
namespace Algorithms {

template <typename T, typename I2, typename I1>
void
toNEWMATMatrix(MatrixWrapper<T, I2, I1> const & src, NEWMAT::Matrix & dst)
{
  typedef MatrixWrapper<T, I2, I1> MW;

  dst.ReSize(src.rows(), src.cols());

  switch (src.getFormat())
  {
    case MatrixFormat::DENSE_ROW_MAJOR:
    {
      typename MW::DenseRowMatrix const & drm = src.getDenseRowMatrix();
      for (int r = 0; r < drm.rows(); ++r)
        for (int c = 0; c < drm.cols(); ++c)
          dst(r, c) = static_cast<NEWMAT::Real>(drm(r, c));

      break;
    }

    case MatrixFormat::DENSE_COLUMN_MAJOR:
    {
      typename MW::DenseColumnMatrix const & dcm = src.getDenseColumnMatrix();
      for (int c = 0; c < dcm.cols(); ++c)
        for (int r = 0; r < dcm.rows(); ++r)
          dst(r, c) = static_cast<NEWMAT::Real>(dcm(r, c));

      break;
    }

    case MatrixFormat::SPARSE_ROW_MAJOR:
    {
      typename MW::SparseRowMatrix const & srm = src.getSparseRowMatrix();
      Array<I1> const & irow  =  srm.getRowIndices();
      Array<I2> const & icol  =  srm.getColumnIndices();
      Array<T>  const & val   =  srm.getValues();

      dst = 0;  // set all entries to zero
      for (size_t r = 0; r < (size_t)srm.rows(); ++r)
        for (I1 i = irow[r]; i < irow[r + 1]; ++i)
        {
          int col = (int)icol[(size_t)i];
          dst(r, col) = static_cast<NEWMAT::Real>(val[(size_t)i]);
        }

      break;
    }

    case MatrixFormat::SPARSE_COLUMN_MAJOR:
    {
      typename MW::SparseColumnMatrix const & scm = src.getSparseColumnMatrix();
      Array<I1> const & icol  =  scm.getColumnIndices();
      Array<I2> const & irow  =  scm.getRowIndices();
      Array<T>  const & val   =  scm.getValues();

      dst = 0;  // set all entries to zero
      for (size_t c = 0; c < (size_t)scm.cols(); ++c)
        for (I1 i = icol[c]; i < icol[c + 1]; ++i)
        {
          int row = (int)irow[(size_t)i];
          dst(row, c) = static_cast<NEWMAT::Real>(val[(size_t)i]);
        }

      break;
    }

    default:
      throw Error("OPTPPNumericalOptimizer: Can't convert to NEWMAT matrix -- unknown source format");
  }
}

template <typename T>
void
toNEWMATColumnVector(Array<T> const & src, NEWMAT::ColumnVector & dst)
{
  dst.ReSize((int)src.size());
  for (size_t i = 0; i < src.size(); ++i)
    dst((int)i) = static_cast<NEWMAT::Real>(src[i]);
}

OPTPP::LinearConstraint *
toOPTPPConstraint(LinearConstraint const & constraint)
{
  alwaysAssertM(Math::isSquare(constraint.getCoefficients())
             && constraint.getCoefficients().rows() == (int)constraint.getConstants().size(),
                "OPTPPNumericalOptimizer: Dimension mismatch in constraint");

  NEWMAT::Matrix lhs;
  toNEWMATMatrix(constraint.getCoefficients(), lhs);

  NEWMAT::ColumnVector rhs;
  toNEWMATColumnVector(constraint.getConstants(), rhs);

  switch (constraint.getCompareOp())
  {
    case CompareOp::EQUAL:
      return new OPTPP::LinearEquation(lhs, rhs);

    case CompareOp::GEQUAL:
      return new OPTPP::LinearInequality(lhs, rhs);  // standard form

    case CompareOp::LEQUAL:
    {
      // Negate LHS
      for (int r = 0; r < lhs.Nrows(); ++r)
        for (int c = 0; c < lhs.Ncols(); ++c)
          lhs(r, c) = -lhs(r, c);

      // Negate RHS
      for (int i = 0; i < rhs.Ncols(); ++i)
        rhs(i) = -rhs(i);

      return new OPTPP::LinearInequality(lhs, rhs);
    }

    default:
      throw Error("OPTPPNumericalOptimizer: Unsupported comparison operator in constraint");
  }
}

OPTPP::NonLinearConstraint *
toOPTPPConstraint(NonlinearConstraint const & constraint)
{
  return nullptr;
}

} // namespace Algorithms
} // namespace Thea
