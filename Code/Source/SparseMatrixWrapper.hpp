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

#ifndef __Thea_SparseMatrixWrapper_hpp__
#define __Thea_SparseMatrixWrapper_hpp__

#include "AbstractCompressedSparseMatrix.hpp"
#include "SparseMatVec.hpp"
#include <type_traits>

namespace Thea {

/**
 * Wrapper for a sparse 2D matrix, implementing a pure virtual interface for passing across shared library boundaries. MatrixT
 * must be an instance of the Eigen::SparseMatrix template.
 */
template <typename MatrixT>
class /* THEA_API */ SparseMatrixWrapper : public AbstractCompressedSparseMatrix<typename MatrixT::value_type>
{
  private:
    static_assert(std::is_base_of<Eigen::SparseMatrixBase<MatrixT>, MatrixT>::value,
                  "SparseMatrixWrapper: Wrapped matrix must be Eigen::SparseMatrix");

    typedef AbstractCompressedSparseMatrix<typename MatrixT::value_type> BaseT;  ///< The base type.

  public:
    THEA_DEF_POINTER_TYPES(SparseMatrixWrapper, std::shared_ptr, std::weak_ptr)

    typedef MatrixT WrappedMatrix;  ///< The wrapped matrix type.
    using typename BaseT::Value;
    using typename BaseT::value_type;

    /**
     * Constructor. The passed matrix must persist as int64 as this class is being actively used, since it is accessed via a
     * pointer.
     */
    SparseMatrixWrapper(MatrixT * mat)
    : m(mat)
    {
      alwaysAssertM(mat, "SparseMatrixWrapper: Must be initialized with an existing matrix");
    }

    /** Destructor. */
    ~SparseMatrixWrapper() {}

    /** Get the wrapped matrix. */
    MatrixT const & getMatrix() const { return *m; }

    /** Get the wrapped matrix. */
    MatrixT & getMatrix() { return *m; }

    // Functions from AbstractMatrix
    int64 rows() const { return (int64)m->rows(); }
    int64 cols() const { return (int64)m->cols(); }
    void setZero() { m->setZero(); }
    int8 isResizable() const { return true; }
    int8 resize(int64 nrows, int64 ncols)
    {
      try  // no exceptions should cross shared library boundaries
      {
        m->resize(nrows, ncols);
        return (m->rows() == nrows && m->cols() == ncols);  // check if it failed without throwing an exception
      }
      THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "SparseMatrixWrapper: Could not resize matrix")
    }

    // Functions from AbstractSparseMatrix
    int64 numStoredElements() const { return m->nonZeros(); }

    // Functions from AbstractCompressedSparseMatrix
    int8 isRowMajor() const { return MatrixT::Flags & Eigen::RowMajorBit; }
    int8 isColumnMajor() const { return !isRowMajor(); }
    int64 innerSize() const { return (int64)m->innerSize(); }
    int64 outerSize() const { return (int64)m->outerSize(); }
    int8 isFullyCompressed() const { return m->isCompressed(); }
    int32 getInnerIndexType() const   { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    int32 getOuterIndexType() const   { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    int32 getNonZeroCountType() const { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    void const * getInnerIndices() const { return m->innerIndexPtr(); }
    void * getInnerIndices() { return m->innerIndexPtr(); }
    void const * getOuterIndices() const { return m->outerIndexPtr(); }
    void * getOuterIndices() { return m->outerIndexPtr(); }
    void const * getNonZeroCounts() const { return m->innerNonZeroPtr(); }
    void * getNonZeroCounts() { return m->innerNonZeroPtr(); }
    Value const * getValues() const { return m->valuePtr(); }
    Value * getValues() { return m->valuePtr(); }

    // Type-casting functions
    AbstractAddressableMatrix<Value> const * asAddressable() const { return NULL; }
    AbstractAddressableMatrix<Value> * asAddressable() { return NULL; }
    AbstractSparseMatrix<Value> const * asSparse() const { return this; }
    AbstractSparseMatrix<Value> * asSparse() { return this; }
    AbstractCompressedSparseMatrix<Value> const * asCompressed() const { return this; }
    AbstractCompressedSparseMatrix<Value> * asCompressed() { return this; }

  private:
    MatrixT * m;  ///< The wrapped matrix.

}; // class SparseMatrixWrapper

} // namespace Thea

#endif
