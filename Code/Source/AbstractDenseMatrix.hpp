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

#ifndef __Thea_AbstractDenseMatrix_hpp__
#define __Thea_AbstractDenseMatrix_hpp__

#include "AbstractAddressableMatrix.hpp"

namespace Thea {

/**
 * Abstract base interface for a 2D dense matrix, assumed to be packed with no gaps or padding in a contiguous memory block.
 * Useful for passing matrices across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ AbstractDenseMatrix : public virtual AbstractAddressableMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(AbstractDenseMatrix, std::shared_ptr, std::weak_ptr)

    /** Is the matrix stored in row-major format? */
    virtual bool isRowMajor() const = 0;

    /** Is the matrix stored in column-major format? */
    virtual bool isColumnMajor() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T const * data() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T * data() = 0;

    /** Set all elements of the matrix to a given value. */
    virtual void fill(T const & value) = 0;

}; // class AbstractDenseMatrix

} // namespace Thea

#endif
