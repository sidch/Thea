//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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
    THEA_DEF_POINTER_TYPES(ImageMatrix, std::shared_ptr, std::weak_ptr)

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

    long rows() const { return image->getHeight(); }
    long cols() const { return image->getWidth(); }
    bool isResizable() const { return true; }
    bool resize(long num_rows, long num_cols) { image->resize(image->getType(), (int)num_cols, (int)num_rows); return true; }
    void setZero()
    {
      // Assume all image channels are zero when they are bitwise zero. Scan width is in bytes.
      std::memset(image->getData(), 0, image->getScanWidth() * image->getHeight());
    }

    Value const & at(long row, long col) const { return ((Value const *)image->getScanLine(row))[col]; }
    Value & at(long row, long col) { return ((Value *)image->getScanLine(row))[col]; }

    // TODO: These could perhaps be made faster by avoiding at()
    void getRow(long row, Value * values) const
    {
      for (long c = 0, ncols = cols(); c < ncols; ++c)
        values[c] = at(row, c);
    }

    void setRow(long row, Value const * values)
    {
      for (long c = 0, ncols = cols(); c < ncols; ++c)
        at(row, c) = values[c];
    }

    void getColumn(long col, Value * values) const
    {
      for (long r = 0, nrows = rows(); r < nrows; ++r)
        values[r] = at(r, col);
    }

    void setColumn(long col, Value const * values)
    {
      for (long r = 0, nrows = rows(); r < nrows; ++r)
        at(r, col) = values[r];
    }

    // Type-casting functions
    AbstractAddressableMatrix<Value> const * asAddressable() const { return this; }
    AbstractAddressableMatrix<Value> * asAddressable() { return this; }
    AbstractSparseMatrix<Value> const * asSparse() const { return NULL; }
    AbstractSparseMatrix<Value> * asSparse() { return NULL; }
    AbstractDenseMatrix<Value> const * asDense() const { return NULL; }
    AbstractDenseMatrix<Value> * asDense() { return NULL; }

  private:
    Image * image;  ///< The wrapped image

}; // class ImageMatrix

} // namespace Thea

#endif
