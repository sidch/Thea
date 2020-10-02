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

#ifndef __Thea_IMatrix_hpp__
#define __Thea_IMatrix_hpp__

#include "Common.hpp"

namespace Thea {

// Forward declarations
template <typename T> class IAddressableMatrix;
template <typename T> class ISparseMatrix;

/** Interface for a 2D matrix, to allow passing matrices and vectors across shared library boundaries. */
template <typename T>
class /* THEA_API */ IMatrix
{
  public:
    THEA_DECL_SMART_POINTERS(IMatrix)

    typedef T Value;       ///< Type of values stored in the matrix.
    typedef T value_type;  ///< Type of values stored in the matrix (STL convention).

    /** Destructor. */
    virtual ~IMatrix() {}

    /** Get the number of rows. */
    virtual int64 THEA_ICALL rows() const = 0;

    /** Get the number of columns. */
    virtual int64 THEA_ICALL cols() const = 0;

    /** Set all elements to zero. */
    virtual void THEA_ICALL setZero() = 0;

    /** Check if the matrix can be resized. */
    virtual int8 THEA_ICALL isResizable() const = 0;

    /**
     * Resize the matrix to new dimensions, if isResizable() returns true. Existing entries will in general <b>not</b> be
     * preserved in the resized matrix.
     *
     * @return True if the matrix was successfully resized, else false.
     */
    virtual int8 THEA_ICALL resize(int64 nrows, int64 ncols) = 0;

    /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual IAddressableMatrix<T> const * THEA_ICALL asAddressable() const = 0;

    /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual IAddressableMatrix<T> * THEA_ICALL asAddressable() = 0;

     /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual ISparseMatrix<T> const * THEA_ICALL asSparse() const = 0;

    /**
     * If the matrix is sparse, get a pointer to a derived interface giving sparse-specific access. Else, return null. Note that
     * a sparse matrix can still be addressable, e.g. if it is stored as a map of <code>(row, col) --> value</code> pairs.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual ISparseMatrix<T> * THEA_ICALL asSparse() = 0;

}; // class IMatrix

} // namespace Thea

#endif
