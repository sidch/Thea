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

#ifndef __Thea_MatrixWrapper_hpp__
#define __Thea_MatrixWrapper_hpp__

#include "IDenseMatrix.hpp"
#include "MatVec.hpp"
#include <type_traits>

namespace Thea {

/**
 * Wrapper for a dense 2D matrix, implementing a pure virtual interface for passing across shared library boundaries. MatrixT
 * must be an instance of the Eigen::Matrix template.
 */
template <typename MatrixT>
class /* THEA_API */ MatrixWrapper : public virtual IDenseMatrix<typename MatrixT::value_type>
{
  private:
    static_assert(std::is_base_of<Eigen::MatrixBase<MatrixT>, MatrixT>::value,
                  "MatrixWrapper: Wrapped matrix must be an Eigen matrix or matrix expression");

    typedef IDenseMatrix<typename MatrixT::value_type> BaseT;  ///< The base type

  public:
    THEA_DECL_SMART_POINTERS(MatrixWrapper)

    typedef MatrixT WrappedMatrix;  ///< The wrapped matrix type.
    using typename BaseT::Value;
    using typename BaseT::value_type;

    /**
     * Constructor. The passed matrix must persist as long as this class is being actively used, since it is accessed via a
     * pointer. If a null pointer is supplied, the functions of this class will have undefined behavior until a valid matrix is
     * supplied via setMatrix or the copy constructor.
     */
    MatrixWrapper(MatrixT * mat = nullptr) : m(mat) {}

    /** Destructor. */
    ~MatrixWrapper() {}

    /** Get the wrapped matrix. */
    MatrixT const & getMatrix() const { return *m; }

    /** Get the wrapped matrix. */
    MatrixT & getMatrix() { return *m; }

    /** Set the wrapped matrix. */
    void setMatrix(MatrixT * mat) { m = mat; }

    // Functions from IMatrix
    int64 THEA_ICALL rows() const { return (int64)m->rows(); }
    int64 THEA_ICALL cols() const { return (int64)m->cols(); }
    void THEA_ICALL setZero() { m->setZero(); }

    int8 THEA_ICALL isResizable() const
    {
      return std::is_base_of<Eigen::PlainObjectBase<MatrixT>, MatrixT>::value  // plain objects (Matrix/Array) own their storage
          && (MatrixT::RowsAtCompileTime == Eigen::Dynamic || MatrixT::ColsAtCompileTime == Eigen::Dynamic);
    }

    int8 THEA_ICALL resize(int64 nrows, int64 ncols)
    {
      if (!isResizable()) { return (m->rows() == nrows && m->cols() == ncols); }
      if (MatrixT::RowsAtCompileTime != Eigen::Dynamic && m->rows() != nrows) { return false; }
      if (MatrixT::ColsAtCompileTime != Eigen::Dynamic && m->cols() != ncols) { return false; }

      try  // no exceptions should cross shared library boundaries
      {
        m->resize(nrows, ncols);
        return (m->rows() == nrows && m->cols() == ncols);  // check if it failed without throwing an exception
      }
      THEA_CATCH(return 0;, ERROR, "%s", "MatrixWrapper: Could not resize matrix")
    }

    // Functions from IAddressableMatrix
    Value const & THEA_ICALL at(int64 row, int64 col) const { return (*m)(row, col); }
    Value & THEA_ICALL mutableAt(int64 row, int64 col) { return (*m)(row, col); }

    // Functions from IRowOrColumnMajorMatrix
    int8 THEA_ICALL isRowMajor() const { return (MatrixT::Flags & Eigen::RowMajorBit); }
    int8 THEA_ICALL isColumnMajor() const { return !isRowMajor(); }

    // Functions from IDenseMatrix
    Value const * THEA_ICALL data() const { return m->data(); }
    Value * THEA_ICALL data() { return m->data(); }
    void THEA_ICALL fill(Value value) { m->fill(value); }

    void THEA_ICALL getRow(int64 row, Value * values) const
    {
      for (intx c = 0, ncols = m->cols(); c < ncols; ++c)
        values[c] = (*m)(row, c);
    }

    void THEA_ICALL setRow(int64 row, Value const * values)
    {
      for (intx c = 0, ncols = m->cols(); c < ncols; ++c)
        (*m)(row, c) = values[c];
    }

    void THEA_ICALL getColumn(int64 col, Value * values) const
    {
      for (intx r = 0, nrows = m->rows(); r < nrows; ++r)
        values[r] = (*m)(r, col);
    }

    void THEA_ICALL setColumn(int64 col, Value const * values)
    {
      for (intx r = 0, nrows = m->rows(); r < nrows; ++r)
        (*m)(r, col) = values[r];
    }

    // Type-casting functions
    IAddressableMatrix<Value> const * THEA_ICALL asAddressable() const { return this; }
    IAddressableMatrix<Value> * THEA_ICALL asAddressable() { return this; }
    ISparseMatrix<Value> const * THEA_ICALL asSparse() const { return nullptr; }
    ISparseMatrix<Value> * THEA_ICALL asSparse() { return nullptr; }
    IDenseMatrix<Value> const * THEA_ICALL asDense() const { return this; }
    IDenseMatrix<Value> * THEA_ICALL asDense() { return this; }

  private:
    MatrixT * m;  ///< The wrapped matrix.

}; // class MatrixWrapper

namespace Math {

/**
 * Convenience function for creating a MatrixWrapper object from a matrix reference, without needing to specify template
 * parameters.
 */
template <typename MatrixT>
MatrixWrapper<MatrixT>
wrapMatrix(Eigen::MatrixBase<MatrixT> & m)
{
  return MatrixWrapper<MatrixT>(static_cast<MatrixT *>(&m));
}

/**
 * Convenience function for creating a MatrixWrapper object from a matrix pointer, without needing to specify template
 * parameters.
 */
template <typename MatrixT>
MatrixWrapper<MatrixT>
wrapMatrix(Eigen::MatrixBase<MatrixT> * m)
{
  return MatrixWrapper<MatrixT>(static_cast<MatrixT *>(m));
}

} // namespace Math

} // namespace Thea

#endif
