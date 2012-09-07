//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_MatrixWrapper_hpp__
#define __Thea_MatrixWrapper_hpp__

#include "Common.hpp"
#include "CompressedSparseMatrix.hpp"
#include "Matrix.hpp"
#include "MatrixFormat.hpp"

namespace Thea {

/** A wrapper for storing a single dense or sparse matrix, in row or column-major form. */
template <typename T, typename Index2DT = int, typename Index1DT = long>
class MatrixWrapper : public virtual BasicMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(MatrixWrapper, shared_ptr, weak_ptr)

    typedef Index2DT  Index2D;  ///< The type of 2D element indices for sparse matrices.
    typedef Index1DT  Index1D;  ///< The type of 1D (flat) element indices for sparse matrices.

    typedef Matrix<T, MatrixLayout::ROW_MAJOR>           DenseRowMatrix;      ///< A matrix in dense row-major form.
    typedef Matrix<T, MatrixLayout::COLUMN_MAJOR>        DenseColumnMatrix;   ///< A matrix in dense column-major form.
    typedef CompressedRowMatrix<T, Index2D, Index1D>     SparseRowMatrix;     ///< A matrix in sparse row-major form.
    typedef CompressedColumnMatrix<T, Index2D, Index1D>  SparseColumnMatrix;  ///< A matrix in sparse column-major form.

    /** Default constructor. */
    MatrixWrapper() : format(MatrixFormat::UNKNOWN) {}

    /** Constructor. Stores a copy of the given matrix in the specified format. */
    template <typename MatrixT> MatrixWrapper(MatrixT const & matrix, MatrixFormat dst_format)
    {
      setMatrix(matrix, dst_format);
    }

    /** Copy constructor. Does a deep copy of the wrapped matrix. */
    MatrixWrapper(MatrixWrapper const & src) { operator=(src); }

    /** Templated copy constructor. Does a deep copy of the wrapped matrix. */
    template <typename U, typename I2, typename I1> MatrixWrapper(MatrixWrapper<U, I2, I1> const & src)
    {
      operator=(src);
    }

    /** Assignment operator. Does a deep copy of the wrapped matrix. */
    MatrixWrapper & operator=(MatrixWrapper const & src)
    {
      MatrixFormat src_format = src.getFormat();
      switch (src_format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:      setMatrix(src.getDenseRowMatrix(), src_format);     break;
        case MatrixFormat::DENSE_COLUMN_MAJOR:   setMatrix(src.getDenseColumnMatrix(), src_format);  break;
        case MatrixFormat::SPARSE_ROW_MAJOR:     setMatrix(src.getSparseRowMatrix(), src_format);    break;
        case MatrixFormat::SPARSE_COLUMN_MAJOR:  setMatrix(src.getSparseColumnMatrix(), src_format); break;
        default: {}
      }

      return *this;
    }

    /** Templated assignment operator. Does a deep copy of the wrapped matrix. */
    template <typename U, typename I2, typename I1> MatrixWrapper & operator=(MatrixWrapper<U, I2, I1> const & src)
    {
      MatrixFormat src_format = src.getFormat();
      switch (src_format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:      setMatrix(src.getDenseRowMatrix(), src_format);     break;
        case MatrixFormat::DENSE_COLUMN_MAJOR:   setMatrix(src.getDenseColumnMatrix(), src_format);  break;
        case MatrixFormat::SPARSE_ROW_MAJOR:     setMatrix(src.getSparseRowMatrix(), src_format);    break;
        case MatrixFormat::SPARSE_COLUMN_MAJOR:  setMatrix(src.getSparseColumnMatrix(), src_format); break;
        default: {}
      }

      return *this;
    }

    /** Store a copy of the given matrix in the specified format. */
    template <typename MatrixT> void setMatrix(MatrixT const & matrix, MatrixFormat dst_format)
    {
      dense_row_matrix.reset();
      dense_column_matrix.reset();
      sparse_column_matrix.reset();
      sparse_row_matrix.reset();

      // Choose a dense row-major storage format by default
      if (dst_format == MatrixFormat::UNKNOWN)
        format = MatrixFormat::DENSE_ROW_MAJOR;
      else
        format = dst_format;

      switch (dst_format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:
          dense_row_matrix = typename DenseRowMatrix::Ptr(new DenseRowMatrix(matrix)); break;

        case MatrixFormat::DENSE_COLUMN_MAJOR:
          dense_column_matrix = typename DenseColumnMatrix::Ptr(new DenseColumnMatrix(matrix)); break;

        case MatrixFormat::SPARSE_ROW_MAJOR:
          sparse_row_matrix = typename SparseRowMatrix::Ptr(new SparseRowMatrix(matrix)); break;

        case MatrixFormat::SPARSE_COLUMN_MAJOR:
          sparse_column_matrix = typename SparseColumnMatrix::Ptr(new SparseColumnMatrix(matrix)); break;

        default: {}
      }
    }

    /** Get the format of the wrapped matrix. */
    MatrixFormat getFormat() const { return format; }

    long numRows() const
    {
      switch (format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:      return dense_row_matrix->numRows();
        case MatrixFormat::DENSE_COLUMN_MAJOR:   return dense_column_matrix->numRows();
        case MatrixFormat::SPARSE_ROW_MAJOR:     return sparse_row_matrix->numRows();
        case MatrixFormat::SPARSE_COLUMN_MAJOR:  return sparse_column_matrix->numRows();
        default: throw Error("MatrixWrapper: Can't get number of rows of unknown matrix type");
      }
    }

    long numColumns() const
    {
      switch (format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:     return dense_row_matrix->numColumns();
        case MatrixFormat::DENSE_COLUMN_MAJOR:  return dense_column_matrix->numColumns();
        case MatrixFormat::SPARSE_ROW_MAJOR:    return sparse_row_matrix->numColumns();
        case MatrixFormat::SPARSE_COLUMN_MAJOR: return sparse_column_matrix->numColumns();
        default: throw Error("MatrixWrapper: Can't get number of columns of unknown matrix type");
      }
    }

    void makeZero()
    {
      switch (format)
      {
        case MatrixFormat::DENSE_ROW_MAJOR:     dense_row_matrix->makeZero(); return;
        case MatrixFormat::DENSE_COLUMN_MAJOR:  dense_column_matrix->makeZero(); return;
        case MatrixFormat::SPARSE_ROW_MAJOR:    sparse_row_matrix->makeZero(); return;
        case MatrixFormat::SPARSE_COLUMN_MAJOR: sparse_column_matrix->makeZero(); return;
        default: throw Error("MatrixWrapper: Can't zero unknown matrix type");
      }
    }

    /** Get the wrapped dense row-major matrix. Call only if getFormat() returns MatrixFormat::DENSE_ROW_MAJOR. */
    DenseRowMatrix const & getDenseRowMatrix() const
    {
      debugAssertM(format == MatrixFormat::DENSE_ROW_MAJOR, "MatrixWrapper: Wrapped matrix is not dense row-major");
      return *dense_row_matrix;
    }

    /** Get the wrapped dense row-major matrix. Call only if getFormat() returns MatrixFormat::DENSE_ROW_MAJOR. */
    DenseRowMatrix & getDenseRowMatrix()
    {
      debugAssertM(format == MatrixFormat::DENSE_ROW_MAJOR, "MatrixWrapper: Wrapped matrix is not dense row-major");
      return *dense_row_matrix;
    }

    /** Get the wrapped dense column-major matrix. Call only if getFormat() returns MatrixFormat::DENSE_COLUMN_MAJOR. */
    DenseColumnMatrix const & getDenseColumnMatrix() const
    {
      debugAssertM(format == MatrixFormat::DENSE_COLUMN_MAJOR, "MatrixWrapper: Wrapped matrix is not dense column-major");
      return *dense_column_matrix;
    }

    /** Get the wrapped dense column-major matrix. Call only if getFormat() returns MatrixFormat::DENSE_COLUMN_MAJOR. */
    DenseColumnMatrix & getDenseColumnMatrix()
    {
      debugAssertM(format == MatrixFormat::DENSE_COLUMN_MAJOR, "MatrixWrapper: Wrapped matrix is not dense column-major");
      return *dense_column_matrix;
    }

    /** Get the wrapped sparse row-major matrix. Call only if getFormat() returns MatrixFormat::SPARSE_ROW_MAJOR. */
    SparseRowMatrix const & getSparseRowMatrix() const
    {
      debugAssertM(format == MatrixFormat::SPARSE_ROW_MAJOR, "MatrixWrapper: Wrapped matrix is not sparse row-major");
      return *sparse_row_matrix;
    }

    /** Get the wrapped sparse row-major matrix. Call only if getFormat() returns MatrixFormat::SPARSE_ROW_MAJOR. */
    SparseRowMatrix & getSparseRowMatrix()
    {
      debugAssertM(format == MatrixFormat::SPARSE_ROW_MAJOR, "MatrixWrapper: Wrapped matrix is not sparse row-major");
      return *sparse_row_matrix;
    }

    /** Get the wrapped sparse column-major matrix. Call only if getFormat() returns MatrixFormat::SPARSE_COLUMN_MAJOR. */
    SparseColumnMatrix const & getSparseColumnMatrix() const
    {
      debugAssertM(format == MatrixFormat::SPARSE_COLUMN_MAJOR, "MatrixWrapper: Wrapped matrix is not sparse column-major");
      return *sparse_column_matrix;
    }

    /** Get the wrapped sparse column-major matrix. Call only if getFormat() returns MatrixFormat::SPARSE_COLUMN_MAJOR. */
    SparseColumnMatrix & getSparseColumnMatrix()
    {
      debugAssertM(format == MatrixFormat::SPARSE_COLUMN_MAJOR, "MatrixWrapper: Wrapped matrix is not sparse column-major");
      return *sparse_column_matrix;
    }

  private:
    MatrixFormat format;
    typename DenseColumnMatrix::Ptr   dense_column_matrix;
    typename DenseRowMatrix::Ptr      dense_row_matrix;
    typename SparseColumnMatrix::Ptr  sparse_column_matrix;
    typename SparseRowMatrix::Ptr     sparse_row_matrix;

}; // class MatrixWrapper

} // namespace Thea

#endif
