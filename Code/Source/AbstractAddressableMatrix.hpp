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
    virtual T & at(int64 row, int64 col) = 0;

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
