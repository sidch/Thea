//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_AbstractMatrix_hpp__
#define __Thea_AbstractMatrix_hpp__

#include "Common.hpp"

namespace Thea {

// Forward declarations
template <typename T> class AbstractAddressableMatrix;
template <typename T> class AbstractSparseMatrix;

/** Abstract base interface for a 2D matrix, to allow passing matrices and vectors across shared library boundaries. */
template <typename T>
class /* THEA_API */ AbstractMatrix
{
  public:
    THEA_DEF_POINTER_TYPES(AbstractMatrix, std::shared_ptr, std::weak_ptr)

    typedef T Value;       ///< Type of values stored in the matrix.
    typedef T value_type;  ///< Type of values stored in the matrix (STL convention).

    /** Destructor. */
    virtual ~AbstractMatrix() {}

    /** Get the number of rows. */
    virtual int64 rows() const = 0;

    /** Get the number of columns. */
    virtual int64 cols() const = 0;

    /** Set all elements to zero. */
    virtual void setZero() = 0;

    /** Check if the matrix can be resized. */
    virtual int8 isResizable() const = 0;

    /**
     * Resize the matrix to new dimensions.
     *
     * @return True if the matrix was successfully resized, else false.
     */
    virtual int8 resize(int64 nrows, int64 ncols) = 0;

    /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractAddressableMatrix<T> const * asAddressable() const = 0;

    /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractAddressableMatrix<T> * asAddressable() = 0;

     /**
     * If the matrix elements are addressable by (row, col) pairs, get a pointer to a derived interface supporting such access.
     * Else, return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractSparseMatrix<T> const * asSparse() const = 0;

    /**
     * If the matrix is sparse, get a pointer to a derived interface giving sparse-specific access. Else, return null. Note that
     * a sparse matrix can still be addressable, e.g. if it is stored as a map of <code>(row, col) --> value</code> pairs.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractSparseMatrix<T> * asSparse() = 0;

}; // class AbstractMatrix

} // namespace Thea

#endif
