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

#ifndef __Thea_Matrix4_hpp__
#define __Thea_Matrix4_hpp__

#include "Math.hpp"
#include "Matrix3.hpp"
#include "MatrixMN.hpp"
#include "VectorN.hpp"
#include <limits>

namespace Thea {

/**
 * Square 4 x 4 matrices on a field T. The matrices are stored internally in row-major form, so row-major access is recommended.
 *
 * @note Several of the functions in this class are taken from the G3D Matrix4 class.
 */
template <typename T>
class /* THEA_API */ MatrixMN<4, 4, T> : public Internal::SquareMatrixN<4, T>
{
  private:
    typedef Internal::SquareMatrixN<4, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /** Initialize all 16 components of the matrix. */
    MatrixMN(T const & m00, T const & m01, T const & m02, T const & m03,
             T const & m10, T const & m11, T const & m12, T const & m13,
             T const & m20, T const & m21, T const & m22, T const & m23,
             T const & m30, T const & m31, T const & m32, T const & m33)
    {
      (*this)(0, 0) = m00; (*this)(0, 1) = m01; (*this)(0, 2) = m02; (*this)(0, 3) = m03;
      (*this)(1, 0) = m10; (*this)(1, 1) = m11; (*this)(1, 2) = m12; (*this)(1, 3) = m13;
      (*this)(2, 0) = m20; (*this)(2, 1) = m21; (*this)(2, 2) = m22; (*this)(2, 3) = m23;
      (*this)(3, 0) = m30; (*this)(3, 1) = m31; (*this)(3, 2) = m32; (*this)(3, 3) = m33;
    }

    /**
     * Initialize from the upper 3x3 submatrix and the first three entries of the last column ("translation" terms). The last
     * row is set to [0, 0, 0, 1].
     */
    MatrixMN(MatrixMN<3, 3, T> const & u3x3, VectorN<3, T> const & trans)
    {
      (*this)(0, 0) = u3x3(0, 0); (*this)(0, 1) = u3x3(0, 1); (*this)(0, 2) = u3x3(0, 2); (*this)(0, 3) = trans[0];
      (*this)(1, 0) = u3x3(1, 0); (*this)(1, 1) = u3x3(1, 1); (*this)(1, 2) = u3x3(1, 2); (*this)(1, 3) = trans[1];
      (*this)(2, 0) = u3x3(2, 0); (*this)(2, 1) = u3x3(2, 1); (*this)(2, 2) = u3x3(2, 2); (*this)(2, 3) = trans[2];
      (*this)(3, 0) = 0; (*this)(3, 1) = 0; (*this)(3, 2) = 0; (*this)(3, 3) = 1;
    }

    /** Get the matrix of cofactors. */
    MatrixMN cofactor() const
    {
      MatrixMN result;

      // We'll use i to incrementally compute -1^(r + c)
      long i = 1;
      for (long r = 0; r < 4; ++r)
      {
        for (long c = 0; c < 4; ++c)
        {
          // Compute the determinant of the 3x3 submatrix
          T det = subDeterminant(r, c);
          result(r, c) = i * det;
          i = -i;
        }

        i = -i;
      }

      return result;
    }

    /** Get the adjoint of the matrix. */
    MatrixMN adjoint() const
    {
      return cofactor().transpose();
    }

    /** Get the determinant of the matrix. */
    T determinant() const
    {
      // Determinant is the dot product of the first row and the first row of cofactors (the first column of the adjoint matrix)
      VectorN<4, T> cofactor0;
      cofactor0[0] =  subDeterminant(0, 0);
      cofactor0[1] = -subDeterminant(0, 1);
      cofactor0[2] =  subDeterminant(0, 2);
      cofactor0[3] = -subDeterminant(0, 3);

      return cofactor0.dot(this->getRow(0));
    }

    /**
     * Invert the matrix. All computations are done using type T, so do <b>not</b> call this function on integer matrices
     * (built-in (POD) integer types will generate assertion failures).
     */
    void invert()
    {
      long col_index[4], row_index[4], pivot[4];
      Internal::invertMatrix(*this, col_index, row_index, pivot);
    }

    /**
     * Get the inverse of the matrix. All computations are done using type T, so do <b>not</b> call this function on integer
     * matrices (built-in (POD) integer types will generate assertion failures). Creates a new matrix for the result, so you
     * might prefer in-place inversion using invert().
     */
    MatrixMN inverse() const { MatrixMN result = *this; result.invert(); return result; }

    /** Get the upper 3x3 submatrix. */
    MatrixMN<3, 3, T> upper3x3() const
    {
      return MatrixMN<3, 3, T>((*this)(0, 0), (*this)(0, 1), (*this)(0, 2),
                               (*this)(1, 0), (*this)(1, 1), (*this)(1, 2),
                               (*this)(2, 0), (*this)(2, 1), (*this)(2, 2));
    }

    /**
     * Create a matrix from its columns.
     *
     * @param cv0 First column of the matrix.
     * @param cv1 Second column of the matrix.
     * @param cv2 Third column of the matrix.
     * @param cv3 Fourth column of the matrix.
     */
    static MatrixMN fromColumns(VectorN<4, T> const & cv0, VectorN<4, T> const & cv1, VectorN<4, T> const & cv2,
                                VectorN<4, T> const & cv3)
    {
      return MatrixMN(cv0[0], cv1[0], cv2[0], cv3[0],
                      cv0[1], cv1[1], cv2[1], cv3[1],
                      cv0[2], cv1[2], cv2[2], cv3[2],
                      cv0[3], cv1[3], cv2[3], cv3[3]);
    }

    /**
     * Create a matrix from its rows.
     *
     * @param rv0 First row of the matrix.
     * @param rv1 Second row of the matrix.
     * @param rv2 Third row of the matrix.
     * @param rv3 Fourth row of the matrix.
     */
    static MatrixMN fromRows(VectorN<4, T> const & rv0, VectorN<4, T> const & rv1, VectorN<4, T> const & rv2,
                             VectorN<4, T> const & rv3)
    {
      return MatrixMN(rv0[0], rv0[1], rv0[2], rv0[3],
                      rv1[0], rv1[1], rv1[2], rv1[3],
                      rv2[0], rv2[1], rv2[2], rv2[3],
                      rv3[0], rv3[1], rv3[2], rv3[3]);
    }

    /** 3D scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(T const & sx, T const & sy, T const & sz)
    {
      return MatrixMN(sx,  0,  0,  0,
                       0, sy,  0,  0,
                       0,  0, sz,  0,
                       0,  0,  0,  1);
    }

    /** 3D scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(VectorN<3, T> const & v)
    {
      return homScaling(v[0], v[1], v[2]);
    }

    /** 3D uniform scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(T const & s)
    {
      return homScaling(s, s, s);
    }

    /** 3D translation matrix in homogeneous coordinates. */
    static MatrixMN homTranslation(T const & tx, T const & ty, T const & tz)
    {
      return MatrixMN(1, 0, 0, tx,
                      0, 1, 0, ty,
                      0, 0, 1, tz,
                      0, 0, 0, 1);
    }

    /** 3D translation matrix in homogeneous coordinates. */
    static MatrixMN homTranslation(VectorN<3, T> const & v)
    {
      return homTranslation(v[0], v[1], v[2]);
    }

    /**
     * Constructs an orthogonal projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
     * the near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases
     * downwards instead, e.g. for screen pixel space.
     */
    static MatrixMN orthogonalProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval,
                                         T const & farval, bool y_increases_upwards = true)
    {
      // Adapted from Mesa.
      // Note that Microsoft (http://msdn.microsoft.com/library/default.asp?url=/library/en-us/opengl/glfunc03_8qnj.asp) and
      // Linux (http://www.xfree86.org/current/glOrtho.3.html) have different matrices shown in their documentation.

      T x =  2 / (right - left);
      T y =  2 / (top - bottom);
      T z = -2 / (farval - nearval);
      T tx = -(right + left) / (right - left);
      T ty = -(top + bottom) / (top - bottom);
      T tz = -(farval + nearval) / (farval - nearval);

      if (!y_increases_upwards)
      {
        y  *= -1;
        ty *= -1;
      }

      return MatrixMN(x, 0, 0, tx,
                      0, y, 0, ty,
                      0, 0, z, tz,
                      0, 0, 0, 1);
    }

    /**
     * Constructs a perspective projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
     * the near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases
     * downwards instead, e.g. for screen pixel space.
     */
    static MatrixMN perspectiveProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval,
                                          T const & farval, bool y_increases_upwards = true)
    {
      T x = (2 * nearval) / (right - left);
      T y = (2 * nearval) / (top - bottom);
      T a = (right + left) / (right - left);
      T b = (top + bottom) / (top - bottom);

      T c, d;
      if (Math::isInfinite(farval))
      {
        // Infinite view frustum
        c = -1;
        d = -2 * nearval;
      }
      else
      {
        c = -(farval + nearval) / (farval - nearval);
        d = -(2 * farval * nearval) / (farval - nearval);
      }

      if (!y_increases_upwards)
      {
        y *= -1;
        b *= -1;
      }

      return MatrixMN(x,  0,  a,  0,
                      0,  y,  b,  0,
                      0,  0,  c,  d,
                      0,  0, -1,  0);
    }

  private:
    /** Get the determinant of the submatrix formed by excluding a row and a column. */
    T subDeterminant(long exclude_row, long exclude_col) const
    {
      // Compute non-excluded row and column indices
      long row[3], col[3];
      for (long i = 0; i < 3; ++i)
      {
        row[i] = i;
        col[i] = i;

        if (i >= exclude_row) ++row[i];
        if (i >= exclude_col) ++col[i];
      }

      // Compute the first row of cofactors
      T cofactor00 = (*this)(row[1], col[1]) * (*this)(row[2], col[2])
                   - (*this)(row[1], col[2]) * (*this)(row[2], col[1]);

      T cofactor10 = (*this)(row[1], col[2]) * (*this)(row[2], col[0])
                   - (*this)(row[1], col[0]) * (*this)(row[2], col[2]);

      T cofactor20 = (*this)(row[1], col[0]) * (*this)(row[2], col[1])
                   - (*this)(row[1], col[1]) * (*this)(row[2], col[0]);

      // Product of the first row and the cofactors along the first row
      return (*this)(row[0], col[0]) * cofactor00
           + (*this)(row[0], col[1]) * cofactor10
           + (*this)(row[0], col[2]) * cofactor20;
    }

}; // class MatrixMN<4, 4, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API MatrixMN<4, 4, Real>;
#endif

/** The default 4x4 real matrix class. */
typedef MatrixMN<4, 4, Real> Matrix4;

/**
 * Multiply a 3-vector by a homogenous 4x4 matrix (assuming the vector has w = 1) and scale back to non-homogeneous coordinates.
 */
template <typename T>
VectorN<3, T>
operator*(MatrixMN<4, 4, T> const & m, VectorN<3, T> const & v)
{
  VectorN<4, T> result = m * VectorN<4, T>(v, 1);

  if (Math::fuzzyEq(result.w(), static_cast<T>(0)))
    return VectorN<3, T>(std::numeric_limits<T>::infinity(),
                         std::numeric_limits<T>::infinity(),
                         std::numeric_limits<T>::infinity());
  else
    return result.xyz() / result.w();
}

} // namespace Thea

#endif
