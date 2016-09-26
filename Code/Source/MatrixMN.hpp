//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_MatrixMN_hpp__
#define __Thea_MatrixMN_hpp__

#include "Common.hpp"
#include "Algorithms/FastCopy.hpp"
#include "AddressableMatrix.hpp"
#include "MatrixInvert.hpp"
#include "VectorN.hpp"
#include <algorithm>
#include <cmath>

namespace Thea {

// Forward declaration
template <long M, long N, typename T> class MatrixMN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for fixed-size M x N matrices, where M (rows) and N (columns) are any <b>positive</b> (non-zero)
 * integers. The matrices are stored internally in row-major form, so row-major access is recommended.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long M, long N, typename T>
class /* THEA_DLL_LOCAL */ MatrixMNBase : public AddressableMatrix<T>
{
    typedef AddressableMatrix<T> BaseT;

  public:
    typedef MatrixMN<M, N, T> MatrixT;  ///< M x N matrix type.

    THEA_DEF_POINTER_TYPES(MatrixT, shared_ptr, weak_ptr)

    /** Default constructor (does not initialize anything). */
    MatrixMNBase() {}

    /** Initialize all components to a single value. */
    explicit MatrixMNBase(T const & fill_value) { fill(fill_value); }

    long numRows() const { return M; }
    long numColumns() const { return N; }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get(). */
    T const & operator()(long row, long col) const { return m[row][col]; }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get() or set(). */
    T & operator()(long row, long col) { return m[row][col]; }

    /** Element access on the matrix unrolled by rows. */
    T const & operator()(long index) const { return m[index / N][index % N]; }

    /** Element access on the matrix unrolled by rows. */
    T & operator()(long index) { return m[index / N][index % N]; }

    /** Element access on the matrix unrolled by rows. */
    T const & operator[](long index) const { return m[index / N][index % N]; }

    /** Element access on the matrix unrolled by rows. */
    T & operator[](long index) { return m[index / N][index % N]; }

    T const & get(long row, long col) const { return m[row][col]; }
    T & getMutable(long row, long col) { return m[row][col]; }
    void set(long row, long col, T const & value) { m[row][col] = value; }

    /** Get the elements of the matrix in row-major ordering. The output array must have at least M * N elements. */
    template <typename S> void getElementsRowMajor(S * elems) const
    {
      long counter = 0;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          elems[counter++] = static_cast<S>(m[i][j]);
    }

    /** Get the elements of the matrix in column-major ordering. The output array must have at least M * N elements. */
    template <typename S> void getElementsColumnMajor(S * elems) const
    {
      long counter = 0;
      for (long j = 0; j < N; ++j)
        for (long i = 0; i < M; ++i)
          elems[counter++] = static_cast<S>(m[i][j]);
    }

    void fill(T const & fill_value)
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          m[i][j] = fill_value;
    }

    /** Equality test. */
    bool operator==(MatrixT const & rhs) const
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          if (m[i][j] != rhs.m[i][j])
            return false;

      return true;
    }

    /** Inequality test. */
    bool operator!=(MatrixT const & rhs) const
    {
      return !operator==(rhs);
    }

    /** Transpose. */
    MatrixMN<N, M, T> transpose() const
    {
      MatrixMN<N, M, T> result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result(j, i) = m[i][j];

      return result;
    }

    void getRow(long row, T * values) const
    {
      for (long j = 0; j < N; ++j)
        values[j] = m[row][j];
    }

    void setRow(long row, T const * values)
    {
      for (long j = 0; j < N; ++j)
        m[row][j] = values[j];
    }

    void getColumn(long col, T * values) const
    {
      for (long i = 0; i < N; ++i)
        values[i] = m[i][col];
    }

    void setColumn(long col, T const * values)
    {
      for (long i = 0; i < N; ++i)
        m[i][col] = values[i];
    }

    /** Get a row of the matrix. */
    VectorN<N, T> getRow(long r) const
    {
      VectorN<N, T> row;
      for (long j = 0; j < N; ++j)
        row[j] = m[r][j];

      return row;
    }

    /** Set a row of the matrix. */
    void setRow(long r, VectorN<N, T> const & values)
    {
      for (long j = 0; j < N; ++j)
        m[r][j] = values[j];
    }

    /** Get a column of the matrix. */
    VectorN<M, T> getColumn(long c) const
    {
      VectorN<M, T> col;
      for (long i = 0; i < M; ++i)
        col[i] = m[i][c];

      return col;
    }

    /** Set a column of the matrix. */
    void setColumn(long c, VectorN<M, T> const & values)
    {
      for (long i = 0; i < M; ++i)
        m[i][c] = values[i];
    }

    /** Negation. */
    MatrixT operator-() const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = -m[i][j];

      return result;
    }

    /** Addition (per-component). */
    MatrixT operator+(MatrixT const & rhs) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = m[i][j] + rhs.m[i][j];

      return result;
    }

    /** Subtraction (per-component). */
    MatrixT operator-(MatrixT const & rhs) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = m[i][j] - rhs.m[i][j];

      return result;
    }

    /** Multiplication. */
    template <long K> MatrixMN<M, K, T> operator*(MatrixMN<N, K, T> const & rhs) const
    {
      MatrixMN<M, K, T> result;
      for (long i = 0; i < M; ++i)
        for (long k = 0; k < K; ++k)
        {
          result(i, k) = m[i][0] * rhs(0, k);

          for (long j = 1; j < N; ++j)
            result(i, k) += m[i][j] * rhs(j, k);
        }

      return result;
    }

    /** Post-multiply by a vector. */
    VectorN<M, T> operator*(VectorN<N, T> const & v) const
    {
      VectorN<M, T> result;
      for (long i = 0; i < M; ++i)
      {
        result[i] = m[i][0] * v[0];

        for (long j = 1; j < N; ++j)
          result[i] += m[i][j] * v[j];
      }

      return result;
    }

    /** Multiply by a scalar. */
    MatrixT operator*(T const & t) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = m[i][j] * t;

      return result;
    }

    /** Divide by a scalar. */
    MatrixT operator/(T const & t) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = m[i][j] / t;

      return result;
    }

    /** Add-and-assign. */
    MatrixT & operator+=(MatrixT const & rhs)
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          m[i][j] += rhs.m[i][j];

      return *static_cast<MatrixT *>(this);
    }

    /** Subtract-and-assign. */
    MatrixT & operator-=(MatrixT const & rhs)
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          m[i][j] -= rhs.m[i][j];

      return *static_cast<MatrixT *>(this);
    }

    /** Multiply in-place by a scalar. */
    MatrixT & operator*=(T const & t)
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          m[i][j] *= t;

      return *static_cast<MatrixT *>(this);
    }

    /** Divide in-place by a scalar. */
    MatrixT & operator/=(T const & t)
    {
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          m[i][j] /= t;

      return *static_cast<MatrixT *>(this);
    }

    T const & min() const { return BaseT::min(); }
    T const & max() const { return BaseT::max(); }

    /** Return a matrix containing the component-wise minima of this matrix and another. */
    MatrixT min(MatrixT const & other) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = std::min(m[i][j], other.m[i][j]);

      return result;
    }

    /** Return a matrix containing the component-wise maxima of this matrix and another. */
    MatrixT max(MatrixT const & other) const
    {
      MatrixT result;
      for (long i = 0; i < M; ++i)
        for (long j = 0; j < N; ++j)
          result.m[i][j] = std::max(m[i][j], other.m[i][j]);

      return result;
    }

    /** Get a matrix containing only zeroes. */
    static MatrixT const & zero() { static MatrixT const z(static_cast<T>(0)); return z; }

  private:
    T m[M][N];  ///< Wrapped 2D array.

}; // class MatrixMNBase

/**
 * <b>[Internal]</b> Base class for fixed-size square N x N matrices, where N is any <b>positive</b> (non-zero) integer. The
 * matrices are stored internally in row-major form, so row-major access is recommended.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ SquareMatrixN : public MatrixMNBase<N, N, T>
{
  private:
    typedef MatrixMNBase<N, N, T> BaseT;

  public:
    typedef typename BaseT::MatrixT MatrixT;

    /** Default constructor (does not initialize anything). */
    SquareMatrixN() {}

    /** Initialize all components to a single value. */
    explicit SquareMatrixN(T const & fill_value) : BaseT(fill_value) {}

    /** Make this an identity matrix. */
    void makeIdentity()
    {
      for (long i = 0; i < N; ++i)
        for (long j = 0; j < N; ++j)
          (*this)(i, j) = (i == j ? 1 : 0);
    }

    /** Get the identity matrix. */
    static MatrixT const & identity() { static MatrixT const id = getIdentityByValue(); return id; }

    /**
     * Construct a diagonal matrix.
     *
     * @param diag A vector containing the elements on the main diagonal.
     */
    static MatrixT fromDiagonal(VectorN<N, T> const & diag)
    {
      MatrixT d = BaseT::zero();
      for (long i = 0; i < N; ++i)
        d(i, i) = diag[i];

      return d;
    }

    /** Replace the matrix with its transpose. */
    void makeTranspose()
    {
      for (long i = 0; i < N; ++i)
        for (long j = i + 1; j < N; ++j)
          std::swap((*this)(i, j), (*this)(j, i));
    }

  private:
    /** Return an identity matrix by value (assertion failure if matrix is not square). */
    static MatrixT getIdentityByValue()
    {
      MatrixT result;
      result.makeIdentity();
      return result;
    }

}; // class SquareMatrixN<N, T>

} // namespace Internal

/**
 * Fixed-size M x N matrices, where M (rows) and N (columns) are any <b>positive</b> (non-zero) integers and T is a field. The
 * matrices are stored internally in row-major form, so row-major access is recommended.
 */
template <long M, long N, typename T>
class /* THEA_API */ MatrixMN : public Internal::MatrixMNBase<M, N, T>
{
  private:
    typedef Internal::MatrixMNBase<M, N, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

}; // class MatrixMN

/**
 * Fixed-size square N x N matrices, where N is any <b>positive</b> (non-zero) integer and T is a field. The matrices are stored
 * internally in row-major form, so row-major access is recommended.
 */
template <long N, typename T>
class /* THEA_API */ MatrixMN<N, N, T> : public Internal::SquareMatrixN<N, T>
{
  private:
    typedef Internal::SquareMatrixN<N, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /**
     * Invert the matrix. All computations are done using type T, so do <b>not</b> call this function on integer matrices
     * (built-in (POD) integer types will generate assertion failures).
     */
    void invert()
    {
      long col_index[N], row_index[N], pivot[N];
      Internal::invertMatrix(*this, &col_index[0], &row_index[0], &pivot[0]);
    }

    /**
     * Get the inverse of the matrix. All computations are done using type T, so do <b>not</b> call this function on integer
     * matrices (built-in (POD) integer types will generate assertion failures). Creates a new matrix for the result, so you
     * might prefer in-place inversion using invert().
     */
    MatrixMN inverse() const { MatrixMN result = *this; result.invert(); return result; }

    /**
     * Post-multiply by a vector of one lower dimension, using homogenous coordinates. The last (homogenous) coordinate of the
     * vector is assumed to be 1, i.e. this is assumed to be a position vector.
     */
    VectorN<N - 1, T> operator*(VectorN<N - 1, T> const & v) const
    {
      VectorN<N - 1, T> result;
      for (long i = 0; i < N - 1; ++i)
      {
        result[i] = (*this)(i, N - 1);
        for (long j = 0; j < N - 1; ++j)
          result[i] += (*this)(i, j) * v[j];
      }

      double w = (*this)(N - 1, N - 1);
      for (long j = 0; j < N - 1; ++j)
        w += (*this)(N - 1, j) * v[j];

      return result / w;
    }

}; // class MatrixMN<M, M, T>

/**
 * A column vector as an M x 1 matrix, where M is any <b>positive</b> (non-zero) integer and T is a field. Generally used to
 * reinterpret VectorN objects as matrices and vice versa. This involves a copy and cannot be used to modify the original
 * objects.
 */
template <long M, typename T>
class /* THEA_API */ MatrixMN<M, 1, T> : public Internal::MatrixMNBase<M, 1, T>
{
  private:
    typedef Internal::MatrixMNBase<M, 1, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /** Initialize from a vector. */
    template <typename U> explicit MatrixMN(Internal::VectorNBase<M, U> const & v)
    {
      Algorithms::fastCopy(&v[0], &v[0] + M, &(*this)(0, 0));
    }

    /** Convert to a vector. */
    VectorN<M, T> toVector() const { return VectorN<M, T>(*this); }

}; // class MatrixMN<M, 1, T>

/**
 * A row vector as a 1 x N matrix, where N is any <b>positive</b> (non-zero) integer and T is a field. Generally used to
 * reinterpret VectorN objects as matrices and vice versa. This involves a copy and cannot be used to modify the original
 * objects.
 */
template <long N, typename T>
class /* THEA_API */ MatrixMN<1, N, T> : public Internal::MatrixMNBase<1, N, T>
{
  private:
    typedef Internal::MatrixMNBase<1, N, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /** Initialize from a vector. */
    template <typename U> explicit MatrixMN(Internal::VectorNBase<N, U> const & v)
    {
      Algorithms::fastCopy(&v[0], &v[0] + N, &(*this)(0, 0));
    }

    /** Convert to a vector. */
    VectorN<N, T> toVector() const { return VectorN<N, T>(*this); }

}; // class MatrixMN<1, N, T>

namespace Internal {

// Definitions of VectorN member functions to convert to N x 1 and 1 x N matrices, and evaluate the outer product.
template <long N, typename T>
MatrixMN<N, 1, T>
VectorNBase<N, T>::toColumnMatrix() const
{
  return MatrixMN<N, 1, T>(*this);
}

template <long N, typename T>
MatrixMN<1, N, T>
VectorNBase<N, T>::toRowMatrix() const
{
  return MatrixMN<1, N, T>(*this);
}

template <long N, typename T>
MatrixMN<N, N, T>
VectorNBase<N, T>::outerProduct(VectorN<N, T> const & v) const
{
  MatrixMN<N, N, T> op;
  for (long i = 0; i < N; ++i)
    for (long j = 0; j < N; ++j)
      op(i, j) = (*this)[i] * v[j];

  return op;
}

} // namespace Internal

/** Pre-multiply by a vector. */
template <long M, long N, typename T>
Thea::VectorN<N, T>
operator*(Thea::VectorN<M, T> const & v, Thea::MatrixMN<M, N, T> const & m)
{
  Thea::VectorN<N, T> result;
  for (long j = 0; j < N; ++j)
  {
    result[j] = m(0, j) * v[0];

    for (long i = 1; i < M; ++i)
      result[j] += m(i, j) * v[i];
  }

  return result;
}

/** Pre-multiply by a scalar. */
template <long M, long N, typename T>
Thea::MatrixMN<M, N, T>
operator*(T const & t, Thea::MatrixMN<M, N, T> const & m)
{
  return m * t;
}

} // namespace Thea

#include "Matrix2.hpp"
#include "Matrix3.hpp"
#include "Matrix4.hpp"

#endif
