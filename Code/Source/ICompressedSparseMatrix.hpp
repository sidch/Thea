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

#ifndef __Thea_ICompressedSparseMatrix_hpp__
#define __Thea_ICompressedSparseMatrix_hpp__

#include "ISparseMatrix.hpp"
#include "IRowOrColumnMajorMatrix.hpp"

namespace Thea {

/**
 * Interface for a 2D sparse matrix in compressed column or compressed row format. Useful for passing matrices across shared
 * library boundaries.
 *
 * @see Eigen::SparseMatrix
 */
template <typename T>
class /* THEA_API */ ICompressedSparseMatrix : public virtual ISparseMatrix<T>, public virtual IRowOrColumnMajorMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(ICompressedSparseMatrix)

    /**
     * Make sure the matrix actually is fully compressed. If not, it will be stored in Eigen's custom format which allows gaps
     * between non-zeros. In the latter case, the array returned by getNonZeroCounts() contains the number of actual
     * non-zeros in each inner segment.
     */
    virtual int8 THEA_ICALL isFullyCompressed() const = 0;

    /** Get the size of the matrix along the inner dimension (rows if column-major, or columns if row-major). */
    virtual int64 THEA_ICALL innerSize() const = 0;

    /** Get the size of the matrix int64 the outer dimension (columns if column-major, or row if row-major). */
    virtual int64 THEA_ICALL outerSize() const = 0;

    /** The integer type used to store inner indices, as per the values in NumericType. */
    virtual int32 THEA_ICALL getInnerIndexType() const = 0;

    /** The integer type used to store outer indices, as per the values in NumericType. */
    virtual int32 THEA_ICALL getOuterIndexType() const = 0;

    /** The integer type used to store per-segment non-zero counts, as per the values in NumericType. */
    virtual int32 THEA_ICALL getNonZeroCountType() const = 0;

    /**
     * Get the array of inner indices. The array has numStoredElements() <i>meaningful</i> entries, but these may not be
     * contiguous if isFullyCompressed() is false. Convert to the correct integer type using getInnerIndexType().
     */
    virtual void const * THEA_ICALL getInnerIndices() const = 0;

    /**
     * Get the array of inner indices. The array has numStoredElements() <i>meaningful</i> entries, but these may not be
     * contiguous if isFullyCompressed() is false. Convert to the correct integer type using getInnerIndexType().
     */
    virtual void * THEA_ICALL getInnerIndices() = 0;

    /**
     * Get the array of outer indices. The array has outerSize() + 1 entries. Convert to the correct integer type using
     * getOuterIndexType().
     */
    virtual void const * THEA_ICALL getOuterIndices() const = 0;

    /**
     * Get the array of outer indices. The array has outerSize() + 1 entries. Convert to the correct integer type using
     * getOuterIndexType().
     */
    virtual void * THEA_ICALL getOuterIndices() = 0;

    /**
     * Get the array of non-zero counts. Each entry is the number of actual non-zeros in the corresponding inner segment. The
     * return value is <b>undefined</b> if isFullyCompressed() is true. Convert to the correct integer type using
     * getNonZeroCountType().
     */
    virtual void const * THEA_ICALL getNonZeroCounts() const = 0;

    /**
     * Get the array of non-zero counts. Each entry is the number of actual non-zeros in the corresponding inner segment. The
     * return value is <b>undefined</b> if isFullyCompressed() is true. Convert to the correct integer type using
     * getNonZeroCountType().
     */
    virtual void * THEA_ICALL getNonZeroCounts() = 0;

    /**
     * Get the array storing the non-zero values. The array has numStoredElements() <i>meaningful</i> entries, but these may not
     * be contiguous if isFullyCompressed() is false.
     */
    virtual T const * THEA_ICALL getValues() const = 0;

    /**
     * Get the array storing the non-zero values. The array has numStoredElements() <i>meaningful</i> entries, but these may not
     * be contiguous if isFullyCompressed() is false.
     */
    virtual T * THEA_ICALL getValues() = 0;

}; // class ICompressedSparseMatrix

} // namespace Thea

#endif
