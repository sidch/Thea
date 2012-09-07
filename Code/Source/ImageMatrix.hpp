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

#include "AddressableMatrix.hpp"
#include "Image.hpp"
#include "ResizableMatrix.hpp"
#include <climits>

namespace Thea {

/**
 * A wrapper class that allows a single-channel (luminance) image to be treated as an addressable matrix. Please ensure the type
 * T matches the type of pixels in the image! The constructor will do a series of checks for this, but checks can be defeated
 * when custom types are involved.
 */
template <typename T>
class /* THEA_API */ ImageMatrix : public AddressableMatrix<T>, public ResizableMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(ImageMatrix, shared_ptr, weak_ptr)

    /** Constructor. The image pointer must remain valid until the matrix is destroyed. */
    ImageMatrix(Image * image_) : image(image_)
    {
      alwaysAssertM(image, "ImageMatrix: Cannot initialize from a null image");
      alwaysAssertM(image->hasByteAlignedPixels(), "ImageMatrix: The image must have byte-aligned pixels");
      alwaysAssertM(image->getBitsPerPixel() / CHAR_BIT == sizeof(T),
                    "ImageMatrix: The number of bytes per pixel does not match the size of the matrix data type");
    }

    long numRows() const { return image->getHeight(); }
    long numColumns() const { return image->getWidth(); }

    void resize(long num_rows, long num_cols) { image->resize(image->getType(), (int)num_cols, (int)num_rows); }

    // This is needed to avoid a Visual Studio warning about inheriting via dominance
    void makeZero() { AddressableMatrix<T>::makeZero(); }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get(). */
    T const & operator()(long row, long col) const { return ((T const *)image->getScanLine(row))[col]; }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get() or set(). */
    T & operator()(long row, long col) { return ((T *)image->getScanLine(row))[col]; }

    T const & get(long row, long col) const { return (*this)(row, col); }
    T & getMutable(long row, long col) { return (*this)(row, col); }
    void set(long row, long col, T const & value) { (*this)(row, col) = value; }

    void fill(T const & value)
    {
      long nr = numRows(), nc = numColumns();

      for (long r = 0; r < nr; ++r)
      {
        T * scanline = static_cast<T *>(image->getScanLine(r));
        std::fill(scanline, scanline + nc, value);
      }
    }

  private:
    Image * image;  ///< The wrapped image

}; // class ImageMatrix

} // namespace Thea

#endif
