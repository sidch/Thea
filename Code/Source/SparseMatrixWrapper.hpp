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
// First version: 2019
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
    THEA_DECL_SMART_POINTERS(SparseMatrixWrapper)

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
    AbstractAddressableMatrix<Value> const * asAddressable() const { return nullptr; }
    AbstractAddressableMatrix<Value> * asAddressable() { return nullptr; }
    AbstractSparseMatrix<Value> const * asSparse() const { return this; }
    AbstractSparseMatrix<Value> * asSparse() { return this; }
    AbstractCompressedSparseMatrix<Value> const * asCompressed() const { return this; }
    AbstractCompressedSparseMatrix<Value> * asCompressed() { return this; }

  private:
    MatrixT * m;  ///< The wrapped matrix.

}; // class SparseMatrixWrapper

} // namespace Thea

#endif
