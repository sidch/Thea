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

#ifndef __Thea_IAddressableMatrix_hpp__
#define __Thea_IAddressableMatrix_hpp__

#include "IMatrix.hpp"

namespace Thea {

// Forward declarations
template <typename T> class IDenseMatrix;

/**
 * Interface for a 2D matrix with elements that may be directly accessed via a (row, column) pair. Useful for passing matrices
 * across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ IAddressableMatrix : public virtual IMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(IAddressableMatrix)

    /**
     * Get a read-only element. Most derived/underlying classes define operator() to access an element quicker, without the
     * virtual function overhead. Use this function only in generic algorithms that need polymorphic access to matrices without
     * using templates, or when accessing matrices across shared library boundaries.
     */
    virtual T const & THEA_ICALL at(int64 row, int64 col) const = 0;

    /**
     * Get an element that can be directly modified. Most derived/underlying classes define operator() to access an element
     * quicker, without the virtual function overhead. Use this function only in generic algorithms that need polymorphic access
     * to matrices without using templates, or when accessing matrices across shared library boundaries.
     */
    virtual T & THEA_ICALL mutableAt(int64 row, int64 col) = 0;

    /** Get a row of the matrix. \a values must be preallocated with cols() elements. */
    virtual void THEA_ICALL getRow(int64 row, T * values) const = 0;

    /** Set a row of the matrix. \a values must contain cols() elements. */
    virtual void THEA_ICALL setRow(int64 row, T const * values) = 0;

    /** Get a column of the matrix. \a values must be preallocated with rows() elements. */
    virtual void THEA_ICALL getColumn(int64 col, T * values) const = 0;

    /** Set a column of the matrix. \a values must contain rows() elements. */
    virtual void THEA_ICALL setColumn(int64 col, T const * values) = 0;

    /**
     * If the matrix is stored as a dense array, get a pointer to a derived interface supporting dense-specific access. Else,
     * return null.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual IDenseMatrix<T> const * THEA_ICALL asDense() const = 0;

    /**
     * If the matrix is stored as a dense array, get a pointer to a derived interface supporting dense-specific access. Else,
     * return null.
     *
     * @note <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid it
     *   and only use <code>static_cast</code> is dangerous.
     */
    virtual IDenseMatrix<T> * THEA_ICALL asDense() = 0;

}; // class IAddressableMatrix

} // namespace Thea

#endif
