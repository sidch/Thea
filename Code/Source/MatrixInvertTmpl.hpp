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

#include "Common.hpp"
#include "AddressableMatrix.hpp"
#include <boost/type_traits/is_integral.hpp>
#include <algorithm>
#include <cmath>

namespace Thea {
namespace Internal {

/**
 * <b>[Internal]</b> Invert a square matrix in place. This function is called by various matrix classes and should not be used
 * directly. Code borrowed from G3D::Matrix::inverseInPlaceGaussJordan(), originally from Numerical Recipes in C (gaussj). The
 * integer arrays \a pivot, \a row_index, and \a col_index are used for bookkeeping on the pivoting and must be
 * preallocated to n elements each, where n is the size (number of rows or number of columns) in the matrix to be inverted.
 *
 * Depending on which preprocessor flags are defined, this function inverts addressable matrices (in which case it uses virtual
 * accessors), or matrices with non-virtual element access via operator(). The latter is expected to be faster.
 *
 * @param mat Input matrix.
 * @param col_index A book-keeping array. Should be preallocated to n elements.
 * @param row_index A book-keeping array. Should be preallocated to n elements.
 * @param pivot A book-keeping array. Should be preallocated to n elements.
 */
template <typename THEA_INVERT_MATRIX_TEMPLATE_TYPE> /* THEA_DLL_LOCAL */
void
THEA_INVERT_MATRIX_FN( THEA_INVERT_MATRIX_TEMPLATE_TYPE & mat, long * col_index, long * row_index, long * pivot )
{
  typedef typename THEA_INVERT_MATRIX_TEMPLATE_TYPE::Value Value;

  alwaysAssertM(mat.isSquare(), "Can't invert non-square matrices");
  alwaysAssertM(!boost::is_integral<Value>::value, "Can't invert integer matrices");

  long n = mat.numRows();

  // Initialize the pivot array to default values.
  static long const NO_PIVOT = -1;
  for (long i = 0; i < n; ++i)
    pivot[i] = NO_PIVOT;

  long col = 0, row = 0;

  // This is the main loop over the columns to be reduced
  //
  // Loop over the columns
  for (long c = 0; c < n; ++c)
  {
    // Find the largest element and use that as a pivot
    Value largest_magnitude = 0;

    // This is the outer loop of the search for a pivot element
    for (long r = 0; r < n; ++r)
    {
      // Unless we've already found the pivot, keep going
      if (pivot[r] != 0)
      {
        // Find the largest pivot
        for (long k = 0; k < n; ++k)
        {
          if (pivot[k] == NO_PIVOT)
          {
            Value mag = std::abs(THEA_MATRIX_GET(mat, r, k));
            if (mag >= largest_magnitude)
            {
              largest_magnitude = mag;
              row = r; col = k;
            }
          }
        }
      }
    }

    pivot[col]++;

    // Interchange columns so that the pivot element is on the diagonal (we'll have to undo this at the end)
    if (row != col)
      for (long k = 0; k < n; ++k)
        std::swap(THEA_MATRIX_GET_MUTABLE(mat, row, k), THEA_MATRIX_GET_MUTABLE(mat, col, k));

    // The pivot is now at [row, col]
    row_index[c] = row;
    col_index[c] = col;

    Value piv = THEA_MATRIX_GET(mat, col, col);
    if (piv == 0)  // FIXME: A more robust comparison?
      throw Error("invertMatrix: Matrix is singular");

    // Divide everything by the pivot (avoid computing the division multiple times).
    Value pivot_inverse = 1 / piv;
    THEA_MATRIX_SET(mat, col, col, 1);

    for (long k = 0; k < n; ++k)
      THEA_MATRIX_GET_MUTABLE(mat, col, k) *= pivot_inverse;

    // Reduce all rows
    for (long r = 0; r < n; ++r)
    {
      // Skip over the pivot row
      if (r != col)
      {
        Value old_value = THEA_MATRIX_GET(mat, r, col);
        THEA_MATRIX_SET(mat, r, col, 0);

        for (long k = 0; k < n; ++k)
          THEA_MATRIX_GET_MUTABLE(mat, r, k) -= THEA_MATRIX_GET(mat, col, k) * old_value;
      }
    }
  }

  // Put the columns back in the correct locations
  for (int i = (int)(n - 1); i >= 0; --i)  // int to handle case when i is -1
    if (row_index[i] != col_index[i])
      for (long k = 0; k < n; ++k)
        std::swap(THEA_MATRIX_GET_MUTABLE(mat, k, row_index[i]), THEA_MATRIX_GET_MUTABLE(mat, k, col_index[i]));
}

} // namespace Internal
} // namespace Thea
