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

#include "ICompressedSparseMatrix.hpp"
#include "SparseMatVec.hpp"
#include <type_traits>

namespace Thea {

/**
 * Wrapper for a sparse 2D matrix, implementing a pure virtual interface for passing across shared library boundaries. MatrixT
 * must be an instance of the Eigen::SparseMatrix template.
 */
template <typename MatrixT>
class /* THEA_API */ SparseMatrixWrapper : public ICompressedSparseMatrix<typename MatrixT::value_type>
{
  private:
    static_assert(std::is_base_of<Eigen::SparseMatrixBase<MatrixT>, MatrixT>::value,
                  "SparseMatrixWrapper: Wrapped matrix must be Eigen::SparseMatrix");

    typedef ICompressedSparseMatrix<typename MatrixT::value_type> BaseT;  ///< The base type.

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

    // Functions from IMatrix
    int64 THEA_ICALL rows() const { return (int64)m->rows(); }
    int64 THEA_ICALL cols() const { return (int64)m->cols(); }
    void THEA_ICALL setZero() { m->setZero(); }
    int8 THEA_ICALL isResizable() const { return true; }
    int8 THEA_ICALL resize(int64 nrows, int64 ncols)
    {
      try  // no exceptions should cross shared library boundaries
      {
        m->resize(nrows, ncols);
        return (m->rows() == nrows && m->cols() == ncols);  // check if it failed without throwing an exception
      }
      THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "SparseMatrixWrapper: Could not resize matrix")
    }

    // Functions from ISparseMatrix
    int64 THEA_ICALL numStoredElements() const { return m->nonZeros(); }

    // Functions from IRowOrColumnMajorMatrix
    int8 THEA_ICALL isRowMajor() const { return (MatrixT::Flags & Eigen::RowMajorBit); }
    int8 THEA_ICALL isColumnMajor() const { return !isRowMajor(); }

    // Functions from ICompressedSparseMatrix
    int64 THEA_ICALL innerSize() const { return (int64)m->innerSize(); }
    int64 THEA_ICALL outerSize() const { return (int64)m->outerSize(); }
    int8 THEA_ICALL isFullyCompressed() const { return m->isCompressed(); }
    int32 THEA_ICALL getInnerIndexType() const   { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    int32 THEA_ICALL getOuterIndexType() const   { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    int32 THEA_ICALL getNonZeroCountType() const { return NumericType::From<typename MatrixT::StorageIndex>::value; }
    void const * THEA_ICALL getInnerIndices() const { return m->innerIndexPtr(); }
    void * THEA_ICALL getInnerIndices() { return m->innerIndexPtr(); }
    void const * THEA_ICALL getOuterIndices() const { return m->outerIndexPtr(); }
    void * THEA_ICALL getOuterIndices() { return m->outerIndexPtr(); }
    void const * THEA_ICALL getNonZeroCounts() const { return m->innerNonZeroPtr(); }
    void * THEA_ICALL getNonZeroCounts() { return m->innerNonZeroPtr(); }
    Value const * THEA_ICALL getValues() const { return m->valuePtr(); }
    Value * THEA_ICALL getValues() { return m->valuePtr(); }

    // Type-casting functions
    IAddressableMatrix<Value> const * THEA_ICALL asAddressable() const { return nullptr; }
    IAddressableMatrix<Value> * THEA_ICALL asAddressable() { return nullptr; }
    ISparseMatrix<Value> const * THEA_ICALL asSparse() const { return this; }
    ISparseMatrix<Value> * THEA_ICALL asSparse() { return this; }
    ICompressedSparseMatrix<Value> const * THEA_ICALL asCompressed() const { return this; }
    ICompressedSparseMatrix<Value> * THEA_ICALL asCompressed() { return this; }

  private:
    MatrixT * m;  ///< The wrapped matrix.

}; // class SparseMatrixWrapper

namespace Math {

/**
 * Convenience function for creating a MatrixWrapper object from a matrix reference, without needing to specify template
 * parameters.
 */
template <typename MatrixT>
SparseMatrixWrapper<MatrixT>
wrapMatrix(Eigen::SparseMatrixBase<MatrixT> & m)
{
  return SparseMatrixWrapper<MatrixT>(static_cast<MatrixT *>(&m));
}

/**
 * Convenience function for creating a MatrixWrapper object from a matrix pointer, without needing to specify template
 * parameters.
 */
template <typename MatrixT>
SparseMatrixWrapper<MatrixT>
wrapMatrix(Eigen::SparseMatrixBase<MatrixT> * m)
{
  return SparseMatrixWrapper<MatrixT>(static_cast<MatrixT *>(m));
}

} // namespace Math

} // namespace Thea

#endif
