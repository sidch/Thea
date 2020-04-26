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
// First version: 2011
//
//============================================================================

#ifndef __Thea_ImageMatrix_hpp__
#define __Thea_ImageMatrix_hpp__

#include "AbstractAddressableMatrix.hpp"
#include "Image.hpp"
#include <algorithm>
#include <climits>
#include <cstring>

namespace Thea {

/**
 * A wrapper class that allows a single-channel (luminance) image to be treated as an addressable matrix. Please ensure the type
 * T matches the type of pixels in the image! The constructor will do a series of checks for this, but checks can be defeated
 * when custom types are involved.
 */
template <typename T>
class /* THEA_API */ ImageMatrix : public AbstractAddressableMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(ImageMatrix)

    using typename AbstractAddressableMatrix<T>::Value;       // == T
    using typename AbstractAddressableMatrix<T>::value_type;  // == T

    /** Constructor. The image pointer must remain valid until the matrix is destroyed. */
    ImageMatrix(Image * image_) : image(image_)
    {
      alwaysAssertM(image, "ImageMatrix: Cannot initialize from a null image");
      alwaysAssertM(image->getDepth() == 1, "ImageMatrix: The image must be 2D");
      alwaysAssertM(image->hasByteAlignedPixels(), "ImageMatrix: The image must have byte-aligned pixels");
      alwaysAssertM(image->getBitsPerPixel() / CHAR_BIT == sizeof(T),
                    "ImageMatrix: The number of bytes per pixel does not match the size of the matrix data type");
    }

    /** Destructor. */
    ~ImageMatrix() {}

    int64 rows() const { return image->getHeight(); }
    int64 cols() const { return image->getWidth(); }
    int8 isResizable() const { return true; }
    int8 resize(int64 num_rows, int64 num_cols) { image->resize(image->getType(), num_cols, num_rows); return true; }

    void setZero()
    {
      // Assume all image channels are zero when they are bitwise zero. Scan width is in bytes.
      std::memset(image->getData(), 0, image->getScanWidth() * image->getHeight());
    }

    Value const & at(int64 row, int64 col) const { return ((Value const *)image->getScanLine(row))[col]; }
    Value & mutableAt(int64 row, int64 col) { return ((Value *)image->getScanLine(row))[col]; }

    // TODO: These could perhaps be made faster by avoiding at() or mutableAt()
    void getRow(int64 row, Value * values) const
    {
      for (intx c = 0, ncols = (intx)cols(); c < ncols; ++c)
        values[c] = at(row, c);
    }

    void setRow(int64 row, Value const * values)
    {
      for (intx c = 0, ncols = (intx)cols(); c < ncols; ++c)
        mutableAt(row, c) = values[c];
    }

    void getColumn(int64 col, Value * values) const
    {
      for (intx r = 0, nrows = (intx)rows(); r < nrows; ++r)
        values[r] = at(r, col);
    }

    void setColumn(int64 col, Value const * values)
    {
      for (intx r = 0, nrows = (intx)rows(); r < nrows; ++r)
        mutableAt(r, col) = values[r];
    }

    // Type-casting functions
    AbstractAddressableMatrix<Value> const * asAddressable() const { return this; }
    AbstractAddressableMatrix<Value> * asAddressable() { return this; }
    AbstractSparseMatrix<Value> const * asSparse() const { return nullptr; }
    AbstractSparseMatrix<Value> * asSparse() { return nullptr; }
    AbstractDenseMatrix<Value> const * asDense() const { return nullptr; }
    AbstractDenseMatrix<Value> * asDense() { return nullptr; }

  private:
    Image * image;  ///< The wrapped image

}; // class ImageMatrix

} // namespace Thea

#endif
