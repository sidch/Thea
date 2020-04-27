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
// First version: 2020
//
//============================================================================

#ifndef __Thea_RowOrColumnMajorMatrix_hpp__
#define __Thea_RowOrColumnMajorMatrix_hpp__

#include "AbstractMatrix.hpp"

namespace Thea {

/** Abstract base interface for a matrix that has either a row-major or a column-major layout. */
template <typename T>
class /* THEA_API */ RowOrColumnMajorMatrix : public virtual AbstractMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(RowOrColumnMajorMatrix)

    /** Is the matrix stored in row-major format? */
    virtual int8 isRowMajor() const = 0;

    /** Is the matrix stored in column-major format? */
    virtual int8 isColumnMajor() const = 0;

}; // class RowOrColumnMajorMatrix

} // namespace Thea

#endif
