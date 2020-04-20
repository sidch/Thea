//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2019, Siddhartha Chaudhuri
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

#include "AbstractDenseMatrix.hpp"
#include "MatVec.hpp"
#include <type_traits>

namespace Thea {

/**
 * Wrapper for a dense 2D matrix, implementing a pure virtual interface for passing across shared library boundaries. MatrixT
 * must be an instance of the Eigen::Matrix template.
 */
template <typename MatrixT>
class /* THEA_API */ MatrixWrapper : public AbstractDenseMatrix<typename MatrixT::value_type>
{
  private:
    static_assert(std::is_base_of<Eigen::MatrixBase<MatrixT>, MatrixT>::value,
                  "MatrixWrapper: Wrapped matrix must be Eigen::Matrix");

    typedef AbstractDenseMatrix<typename MatrixT::value_type> BaseT;  ///< The base type

  public:
    THEA_DECL_SMART_POINTERS(MatrixWrapper)

    typedef MatrixT WrappedMatrix;  ///< The wrapped matrix type.
    using typename BaseT::Value;
    using typename BaseT::value_type;

    /**
     * Constructor. The passed matrix must persist as long as this class is being actively used, since it is accessed via a
     * pointer.
     */
    MatrixWrapper(MatrixT * mat)
    : m(mat)
    {
      alwaysAssertM(mat, "MatrixWrapper: Must be initialized with an existing matrix");
    }

    /** Destructor. */
    ~MatrixWrapper() {}

    /** Get the wrapped matrix. */
    MatrixT const & getMatrix() const { return *m; }

    /** Get the wrapped matrix. */
    MatrixT & getMatrix() { return *m; }

    // Functions from AbstractMatrix
    int64 rows() const { return (int64)m->rows(); }
    int64 cols() const { return (int64)m->cols(); }
    void setZero() { m->setZero(); }

    int8 isResizable() const
    {
      return MatrixT::RowsAtCompileTime == Eigen::Dynamic || MatrixT::ColsAtCompileTime == Eigen::Dynamic;
    }

    int8 resize(int64 nrows, int64 ncols)
    {
      try  // no exceptions should cross shared library boundaries
      {
        m->resize(nrows, ncols);
        return (m->rows() == nrows && m->cols() == ncols);  // check if it failed without throwing an exception
      }
      THEA_STANDARD_CATCH_BLOCKS(return 0;, ERROR, "%s", "MatrixWrapper: Could not resize matrix")
    }

    // Functions from AbstractAddressableMatrix
    Value const & at(int64 row, int64 col) const { return (*m)(row, col); }
    Value & mutableAt(int64 row, int64 col) { return (*m)(row, col); }

    // Functions from AbstractDenseMatrix
    int8 isRowMajor() const { return (MatrixT::Flags & Eigen::RowMajorBit); }
    int8 isColumnMajor() const { return isRowMajor(); }
    Value const * data() const { return m->data(); }
    Value * data() { return m->data(); }
    void fill(Value const & value) { m->fill(value); }

    void getRow(int64 row, Value * values) const
    {
      for (intx c = 0, ncols = m->cols(); c < ncols; ++c)
        values[c] = (*m)(row, c);
    }

    void setRow(int64 row, Value const * values)
    {
      for (intx c = 0, ncols = m->cols(); c < ncols; ++c)
        (*m)(row, c) = values[c];
    }

    void getColumn(int64 col, Value * values) const
    {
      for (intx r = 0, nrows = m->rows(); r < nrows; ++r)
        values[r] = (*m)(r, col);
    }

    void setColumn(int64 col, Value const * values)
    {
      for (intx r = 0, nrows = m->rows(); r < nrows; ++r)
        (*m)(r, col) = values[r];
    }

    // Type-casting functions
    AbstractAddressableMatrix<Value> const * asAddressable() const { return this; }
    AbstractAddressableMatrix<Value> * asAddressable() { return this; }
    AbstractSparseMatrix<Value> const * asSparse() const { return nullptr; }
    AbstractSparseMatrix<Value> * asSparse() { return nullptr; }
    AbstractDenseMatrix<Value> const * asDense() const { return this; }
    AbstractDenseMatrix<Value> * asDense() { return this; }

  private:
    MatrixT * m;  ///< The wrapped matrix.

}; // class MatrixWrapper

} // namespace Thea

#endif
