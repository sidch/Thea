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

  dst.ReSize(src.numRows(), src.numColumns());

  switch (src.getFormat())
  {
    case MatrixFormat::DENSE_ROW_MAJOR:
    {
      typename MW::DenseRowMatrix const & drm = src.getDenseRowMatrix();
      for (int r = 0; r < drm.numRows(); ++r)
        for (int c = 0; c < drm.numColumns(); ++c)
          dst(r, c) = static_cast<NEWMAT::Real>(drm(r, c));

      break;
    }

    case MatrixFormat::DENSE_COLUMN_MAJOR:
    {
      typename MW::DenseColumnMatrix const & dcm = src.getDenseColumnMatrix();
      for (int c = 0; c < dcm.numColumns(); ++c)
        for (int r = 0; r < dcm.numRows(); ++r)
          dst(r, c) = static_cast<NEWMAT::Real>(dcm(r, c));

      break;
    }

    case MatrixFormat::SPARSE_ROW_MAJOR:
    {
      typename MW::SparseRowMatrix const & srm = src.getSparseRowMatrix();
      TheaArray<I1> const & irow  =  srm.getRowIndices();
      TheaArray<I2> const & icol  =  srm.getColumnIndices();
      TheaArray<T>  const & val   =  srm.getValues();

      dst = 0;  // set all entries to zero
      for (array_size_t r = 0; r < (array_size_t)srm.numRows(); ++r)
        for (I1 i = irow[r]; i < irow[r + 1]; ++i)
        {
          int col = (int)icol[(array_size_t)i];
          dst(r, col) = static_cast<NEWMAT::Real>(val[(array_size_t)i]);
        }

      break;
    }

    case MatrixFormat::SPARSE_COLUMN_MAJOR:
    {
      typename MW::SparseColumnMatrix const & scm = src.getSparseColumnMatrix();
      TheaArray<I1> const & icol  =  scm.getColumnIndices();
      TheaArray<I2> const & irow  =  scm.getRowIndices();
      TheaArray<T>  const & val   =  scm.getValues();

      dst = 0;  // set all entries to zero
      for (array_size_t c = 0; c < (array_size_t)scm.numColumns(); ++c)
        for (I1 i = icol[c]; i < icol[c + 1]; ++i)
        {
          int row = (int)irow[(array_size_t)i];
          dst(row, c) = static_cast<NEWMAT::Real>(val[(array_size_t)i]);
        }

      break;
    }

    default:
      throw Error("OPTPPNumericalOptimizer: Can't convert to NEWMAT matrix -- unknown source format");
  }
}

template <typename T>
void
toNEWMATColumnVector(TheaArray<T> const & src, NEWMAT::ColumnVector & dst)
{
  dst.ReSize((int)src.size());
  for (array_size_t i = 0; i < src.size(); ++i)
    dst((int)i) = static_cast<NEWMAT::Real>(src[i]);
}

OPTPP::LinearConstraint *
toOPTPPConstraint(LinearConstraint const & constraint)
{
  alwaysAssertM(constraint.getCoefficients().isSquare()
             && constraint.getCoefficients().numRows() == (int)constraint.getConstants().size(),
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
  return NULL;
}

} // namespace Algorithms
} // namespace Thea
