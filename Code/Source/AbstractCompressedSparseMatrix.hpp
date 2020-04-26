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
    THEA_DECL_SMART_POINTERS(AbstractCompressedSparseMatrix)

    /** Is the matrix stored in compressed row format? */
    virtual int8 isRowMajor() const = 0;

    /** Is the matrix stored in compressed column format? */
    virtual int8 isColumnMajor() const = 0;

    /**
     * Make sure the matrix actually is fully compressed. If not, it will be stored in Eigen's custom format which allows gaps
     * between non-zeros. In the latter case, the array returned by getNonZeroCounts() contains the number of actual
     * non-zeros in each inner segment.
     */
    virtual int8 isFullyCompressed() const = 0;

    /** Get the size of the matrix along the inner dimension (rows if column-major, or columns if row-major). */
    virtual int64 innerSize() const = 0;

    /** Get the size of the matrix int64 the outer dimension (columns if column-major, or row if row-major). */
    virtual int64 outerSize() const = 0;

    /** The integer type used to store inner indices, as per the values in NumericType. */
    virtual int32 getInnerIndexType() const = 0;

    /** The integer type used to store outer indices, as per the values in NumericType. */
    virtual int32 getOuterIndexType() const = 0;

    /** The integer type used to store per-segment non-zero counts, as per the values in NumericType. */
    virtual int32 getNonZeroCountType() const = 0;

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
