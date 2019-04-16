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

#ifndef __Thea_AbstractSparseMatrix_hpp__
#define __Thea_AbstractSparseMatrix_hpp__

#include "AbstractMatrix.hpp"

namespace Thea {

// Forward declarations
template <typename T> class AbstractCompressedSparseMatrix;

/** Abstract base interface for a 2D sparse matrix. Useful for passing matrices across shared library boundaries. */
template <typename T>
class /* THEA_API */ AbstractSparseMatrix : public virtual AbstractMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(AbstractSparseMatrix, std::shared_ptr, std::weak_ptr)

    /**
     * Get the number of entries actually stored in the matrix. These are often called "non-zeros", though they may actually
     * have the numeric value 0.
     */
    virtual long numStoredElements() const = 0;

    /**
     * If the matrix is stored in compressed column or row format, get a pointer to a derived interface supporting access
     * specific to that format. Else, return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractCompressedSparseMatrix<T> const * asCompressed() const = 0;

    /**
     * If the matrix is stored in compressed column or row format, get a pointer to a derived interface supporting access
     * specific to that format. Else, return null.
     *
     * @note: <code>dynamic_cast</code> does not work reliably across shared library boundaries, and relying on users to avoid
     *   it and only use <code>static_cast</code> is dangerous.
     */
    virtual AbstractCompressedSparseMatrix<T> * asCompressed() = 0;

}; // class AbstractSparseMatrix

} // namespace Thea

#endif
