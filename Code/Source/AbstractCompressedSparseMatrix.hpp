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

#ifndef __Thea_AbstractCompressedSparseMatrix_hpp__
#define __Thea_AbstractCompressedSparseMatrix_hpp__

#include "AbstractSparseMatrix.hpp"

namespace Thea {

/**
 * Abstract base interface for a 2D sparse matrix in compressed column or compressed row format. Useful for passing matrices
 * across shared library boundaries.
 *
 * @see Eigen::SparseMatrix
 */
template <typename T>
class /* THEA_API */ AbstractCompressedSparseMatrix : public virtual AbstractSparseMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(AbstractCompressedSparseMatrix, std::shared_ptr, std::weak_ptr)

    /** Is the matrix stored in compressed row format? */
    virtual bool isRowMajor() const = 0;

    /** Is the matrix stored in compressed column format? */
    virtual bool isColumnMajor() const = 0;

    /**
     * Make sure the matrix actually is fully compressed. If not, it will be stored in Eigen's custom format which allows gaps
     * between non-zeros. In the latter case, the array returned by getNonZeroCounts() contains the number of actual
     * non-zeros in each inner segment.
     */
    virtual bool isFullyCompressed() const = 0;

    /** Get the size of the matrix along the inner dimension (rows if column-major, or columns if row-major). */
    virtual long innerSize() const = 0;

    /** Get the size of the matrix long the outer dimension (columns if column-major, or row if row-major). */
    virtual long outerSize() const = 0;

    /** The integer type used to store inner indices, as per the values in NumericType. */
    virtual int getInnerIndexType() const = 0;

    /** The integer type used to store outer indices, as per the values in NumericType. */
    virtual int getOuterIndexType() const = 0;

    /** The integer type used to store per-segment non-zero counts, as per the values in NumericType. */
    virtual int getNonZeroCountType() const = 0;

    /**
     * Get the array of inner indices. The array has numStoredElements() <i>meaningful</i> entries, but these may not be
     * contiguous if isFullyCompressed() is false. Convert to the correct integer type using getInnerIndexType().
     */
    virtual void const * getInnerIndices() const = 0;

    /**
     * Get the array of inner indices. The array has numStoredElements() <i>meaningful</i> entries, but these may not be
     * contiguous if isFullyCompressed() is false. Convert to the correct integer type using getInnerIndexType().
     */
    virtual void * getInnerIndices() = 0;

    /**
     * Get the array of outer indices. The array has outerSize() + 1 entries. Convert to the correct integer type using
     * getOuterIndexType().
     */
    virtual void const * getOuterIndices() const = 0;

    /**
     * Get the array of outer indices. The array has outerSize() + 1 entries. Convert to the correct integer type using
     * getOuterIndexType().
     */
    virtual void * getOuterIndices() = 0;

    /**
     * Get the array of non-zero counts. Each entry is the number of actual non-zeros in the corresponding inner segment. The
     * return value is <b>undefined</b> if isFullyCompressed() is true. Convert to the correct integer type using
     * getNonZeroCountType().
     */
    virtual void const * getNonZeroCounts() const = 0;

    /**
     * Get the array of non-zero counts. Each entry is the number of actual non-zeros in the corresponding inner segment. The
     * return value is <b>undefined</b> if isFullyCompressed() is true. Convert to the correct integer type using
     * getNonZeroCountType().
     */
    virtual void * getNonZeroCounts() = 0;

    /**
     * Get the array storing the non-zero values. The array has numStoredElements() <i>meaningful</i> entries, but these may not
     * be contiguous if isFullyCompressed() is false.
     */
    virtual T const * getValues() const = 0;

    /**
     * Get the array storing the non-zero values. The array has numStoredElements() <i>meaningful</i> entries, but these may not
     * be contiguous if isFullyCompressed() is false.
     */
    virtual T * getValues() = 0;

}; // class AbstractCompressedSparseMatrix

} // namespace Thea

#endif
