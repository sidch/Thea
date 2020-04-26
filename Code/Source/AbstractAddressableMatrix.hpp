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

#ifndef __Thea_AbstractAddressableMatrix_hpp__
#define __Thea_AbstractAddressableMatrix_hpp__

#include "AbstractMatrix.hpp"

namespace Thea {

// Forward declarations
template <typename T> class AbstractDenseMatrix;

/**
 * Abstract base interface for a 2D matrix with elements that may be directly accessed via a (row, column) pair. Useful for
 * passing matrices across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ AbstractAddressableMatrix : public virtual AbstractMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(AbstractAddressableMatrix)

    /**
     * Get a read-only element. Most derived/underlying classes define operator() to access an element quicker, without the
     * virtual function overhead. Use this function only in generic algorithms that need polymorphic access to matrices without
     * using templates, or when accessing matrices across shared library boundaries.
     */
    virtual T const & at(int64 row, int64 col) const = 0;

    /**
     * Get an element that can be directly modified. Most derived/underlying classes define operator() to access an element
     * quicker, without the virtual function overhead. Use this function only in generic algorithms that need polymorphic access
     * to matrices without using templates, or when accessing matrices across shared library boundaries.
     */
    virtual T & mutableAt(int64 row, int64 col) = 0;

    /** Get a row of the matrix. \a values must be preallocated with cols() elements. */
    virtual void getRow(int64 row, T * values) const = 0;

    /** Set a row of the matrix. \a values must contain cols() elements. */
    virtual void setRow(int64 row, T const * values) = 0;

    /** Get a column of the matrix. \a values must be preallocated with rows() elements. */
    virtual void getColumn(int64 col, T * values) const = 0;

    /** Set a column of the matrix. \a values must contain rows() elements. */
    virtual void setColumn(int64 col, T const * values) = 0;

    /**
     * If the matrix is stored as a dense array, get a pointer to a derived interface supporting dense-specific access. Else,
     * return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractDenseMatrix<T> const * asDense() const = 0;

    /**
     * If the matrix is stored as a dense array, get a pointer to a derived interface supporting dense-specific access. Else,
     * return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractDenseMatrix<T> * asDense() = 0;

}; // class AbstractAddressableMatrix

} // namespace Thea

#endif
