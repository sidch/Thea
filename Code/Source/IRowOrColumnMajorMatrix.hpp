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

#ifndef __Thea_IRowOrColumnMajorMatrix_hpp__
#define __Thea_IRowOrColumnMajorMatrix_hpp__

#include "IMatrix.hpp"

namespace Thea {

/** Interface for a matrix that has either a row-major or a column-major layout. */
template <typename T>
class /* THEA_API */ IRowOrColumnMajorMatrix : public virtual IMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(IRowOrColumnMajorMatrix)

    /** Is the matrix stored in row-major format? */
    virtual int8 THEA_ICALL isRowMajor() const = 0;

    /** Is the matrix stored in column-major format? */
    virtual int8 THEA_ICALL isColumnMajor() const = 0;

}; // class IRowOrColumnMajorMatrix

} // namespace Thea

#endif
