//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, * except for portions
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

#ifndef __Thea_Matrix_hpp__
#define __Thea_Matrix_hpp__

#include "Common.hpp"
#include "Algorithms/FastCopy.hpp"
#include "AddressableMatrix.hpp"
#include "Array.hpp"
#include "Math.hpp"
#include "MatrixInvert.hpp"
#include "ResizableMatrix.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <algorithm>
#include <memory>

namespace Thea {

// Forward declarations
template <typename T, MatrixLayout::Value L, typename Index2DT, typename Index1DT> class CompressedSparseMatrix;
template <typename T, typename Index2DT, typename Index1DT> class CompressedRowMatrix;
template <typename T, typename Index2DT, typename Index1DT> class CompressedColumnMatrix;

/** A standard dense 2D matrix, either row-major or column-major. */
template < typename T, MatrixLayout::Value L = MatrixLayout::ROW_MAJOR, typename AllocT = std::allocator<T> >
class /* THEA_API */ Matrix : public AddressableMatrix<T>, public ResizableMatrix<T>
{
  private:
    typedef AddressableMatrix<T> AddressableBaseT;

    AllocT allocator;

  public:
    THEA_DEF_POINTER_TYPES(Matrix, shared_ptr, weak_ptr)

    static MatrixLayout::Value const Layout = L;  ///< The layout of the matrix (row-major or column-major).

    /** Constructs an empty matrix. */
    Matrix() : num_rows(0), num_cols(0), values(NULL), owns_memory(true) {}

    /** Constructs a matrix of a given size, filled with a given value. */
    Matrix(long num_rows_, long num_cols_) : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      resize(num_rows_, num_cols_);
    }

    /** Constructs a matrix of a given size, filled with a given value. */
    Matrix(long num_rows_, long num_cols_, T const & fill_value) : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      alwaysAssertM(num_rows_ >= 0 && num_cols_ >= 0, "Matrix: Dimensions must be non-negative");
      resize(num_rows_, num_cols_, fill_value);
    }

    /**
     * Constructs this matrix as a wrapper for an existing block of storage. The matrix thus created is <b>not resizable</b>,
     * and the memory block will <b>not be freed</b> when the Matrix object is destroyed.
     */
    Matrix(T * existing_data, long num_rows_, long num_cols_)
    : num_rows(num_rows_), num_cols(num_cols_), values(existing_data), owns_memory(false)
    {
      alwaysAssertM(num_rows_ >= 0 && num_cols_ >= 0, "Matrix: Dimensions must be non-negative");
    }

    /** Construct a square diagonal matrix, given the values on the diagonal. */
    static Matrix fromDiagonal(TheaArray<T> const & diagonal)
    {
      Matrix m((long)diagonal.size(), (long)diagonal.size(), static_cast<T>(0));
      for (size_t i = 0; i < diagonal.size(); ++i)
        m((long)i, (long)i) = diagonal[i];

      return m;
    }

    /** Copy constructor. */
    Matrix(Matrix const & m) : num_rows(0), num_cols(0), values(NULL), owns_memory(true) { *this = m; }

    /** Initialize from an addressable matrix. */
    template <typename MatrixT> explicit
    Matrix(MatrixT const & src, typename boost::enable_if< boost::is_base_of< AddressableMatrix<typename MatrixT::Value>,
                                                                              MatrixT > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(src.numRows(), src.numColumns());
      else
      {
        alwaysAssertM(num_rows == src.numRows() && num_cols == src.numColumns(),
                      "Matrix: A wrapper matrix cannot be assigned a value of different dimensions");
      }

      for (long r = 0; r < num_rows; ++r)
        for (long c = 0; c < num_cols; ++c)
          (*this)(r, c) = static_cast<T>(src(r, c));
    }

    /**
     * Initialize from a compressed row matrix. The class T must be <code>static_cast</code>able from an integer, specifically
     * from zero.
     *
     * <b>Implementation note:</b> Since this is a templated function, the forward declaration is enough. The body of
     * CompressedRowMatrix will necessarily be available at the point of instantiation.
     */
    template <typename MatrixT> explicit
    Matrix(MatrixT const & crm, typename boost::enable_if< boost::is_base_of< CompressedRowMatrix< typename MatrixT::Value,
                                                                                                   typename MatrixT::Index2D,
                                                                                                   typename MatrixT::Index1D >,
                                                                              MatrixT > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(crm.numRows(), crm.numColumns());
      else
      {
        alwaysAssertM(num_rows == crm.numRows() && num_cols == crm.numColumns(),
                      "Matrix: A wrapper matrix cannot be assigned a value of different dimensions");
      }

      fill(static_cast<T>(0));

      size_t curr_pos = 0;
      for (long row = 0; row < num_rows; ++row)
      {
        long num_elems = (long)(crm.getRowIndices()[row + 1] - crm.getRowIndices()[row]);
        for (long e = 0; e < num_elems; ++e, ++curr_pos)
        {
          long col = (long)crm.getColumnIndices()[curr_pos];
          (*this)(row, col) = static_cast<T>(crm.getValues()[curr_pos]);
        }
      }
    }

    /**
     * Initialize from a compressed column matrix. The class T must be <code>static_cast</code>able from an integer,
     * specifically from zero.
     *
     * <b>Implementation note:</b> Since this is a templated function, the forward declaration is enough. The body of
     * CompressedColumnMatrix will necessarily be available at the point of instantiation.
     */
    template <typename MatrixT> explicit
    Matrix(MatrixT const & ccm,
           typename boost::enable_if< boost::is_base_of< CompressedColumnMatrix< typename MatrixT::Value,
                                                                                 typename MatrixT::Index2D,
                                                                                 typename MatrixT::Index1D >,
                                                         MatrixT > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(ccm.numRows(), ccm.numColumns());
      else
      {
        alwaysAssertM(num_rows == ccm.numRows() && num_cols == ccm.numColumns(),
                      "Matrix: A wrapper matrix cannot be assigned a value of different dimensions");
      }

      fill(static_cast<T>(0));

      size_t curr_pos = 0;
      for (long col = 0; col < num_cols; ++col)
      {
        long num_elems = (long)(ccm.getColumnIndices()[col + 1] - ccm.getColumnIndices()[col]);
        for (long e = 0; e < num_elems; ++e, ++curr_pos)
        {
          long row = (long)ccm.getRowIndices()[curr_pos];
          (*this)(row, col) = static_cast<T>(ccm.getValues()[curr_pos]);
        }
      }
    }

    /** Destructor. */
    ~Matrix()
    {
      if (owns_memory) allocator.deallocate(values, (size_t)this->numElements());
    }

    /** Assignment operator. */
    Matrix & operator=(Matrix const & src)
    {
      if (owns_memory)
        resize(src.numRows(), src.numColumns());
      else
      {
        alwaysAssertM(num_rows == src.numRows() && num_cols == src.numColumns(),
                      "Matrix: A wrapper matrix cannot be assigned a value of different dimensions");
      }

      Algorithms::fastCopy(src.values, src.values + this->numElements(), values);
      return *this;
    }

    /** Get the layout of the matrix (row or column major) */
    static MatrixLayout getLayout() { return L; }

    long numRows() const { return num_rows; }
    long numColumns() const { return num_cols; }

    void resize(long num_rows_, long num_cols_)
    {
      alwaysAssertM(num_rows_ >= 0 && num_cols_ >= 0, "Matrix: Dimensions must be non-negative");
      alwaysAssertM(owns_memory, "Matrix: This matrix does not own its storage block and can't be resized");

      if (num_rows != num_rows_ || num_cols != num_cols_)
      {
        long old_num_elems = this->numElements();
        long new_num_elems = num_rows_ * num_cols_;
        if (old_num_elems != new_num_elems)
        {
          allocator.deallocate(values, (size_t)old_num_elems);

          if (new_num_elems > 0) values = allocator.allocate((size_t)new_num_elems);
          else                   values = NULL;
        }

        num_rows = num_rows_;
        num_cols = num_cols_;
      }
    }

    /** Resize the matrix and set all entries to a given value. */
    void resize(long num_rows_, long num_cols_, T const & fill_value)
    {
      resize(num_rows_, num_cols_);
      fill(fill_value);
    }

    /**
     * Append a single (uninitialized) row to a row-major matrix. If the matrix is not row-major, or if the matrix does not own
     * its memory block, an assertion failure occurs.
     *
     * @warning This is a slow operation especially for large matrices! Avoid calling it repeatedly. Instead, preallocate memory
     *   for the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to
     *   resize().
     *
     * @see appendRows()
     */
    void appendRow()
    {
      appendRows(1);
    }

    /**
     * Append one or more (uninitialized) rows to a row-major matrix. If the matrix is not row-major, or if the matrix does not
     * own its memory block, an assertion failure occurs. \a num_rows_to_append must be non-negative.
     *
     * @warning This is very slow if you repeatedly add small numbers of rows to a large matrix. Instead, preallocate memory for
     *   the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to resize().
     *
     * @see appendRow()
     */
    void appendRows(long num_rows_to_append)
    {
      alwaysAssertM(L == MatrixLayout::ROW_MAJOR, "Matrix: Cannot append row(s) to a matrix that is not row-major");
      alwaysAssertM(owns_memory, "Matrix: Cannot append row(s) to a matrix that does not own its memory block");
      alwaysAssertM(num_rows_to_append >= 0, "Matrix: Cannot append a negative number of rows to a matrix");

      if (num_rows_to_append <= 0)
        return;

      long new_num_rows = num_rows + num_rows_to_append;
      long new_num_elems = new_num_rows * num_cols;
      T * new_values = allocator.allocate((size_t)new_num_elems);

      if (!this->isEmpty())
      {
        size_t old_num_elems = (size_t)this->numElements();
        Algorithms::fastCopy(values, values + old_num_elems, new_values);
        allocator.deallocate(values, old_num_elems);
      }

      values = new_values;
      num_rows = new_num_rows;
    }

    /**
     * Append a single (uninitialized) column to a column-major matrix. If the matrix is not column-major, or if the matrix does
     * not own its memory block, an assertion failure occurs.
     *
     * @warning This is a slow operation especially for large matrices! Avoid calling it repeatedly. Instead, preallocate memory
     *   for the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to
     *   resize().
     *
     * @see appendColumns()
     */
    void appendColumn()
    {
      appendColumns(1);
    }

    /**
     * Append one or more (uninitialized) columns to a column-major matrix. If the matrix is not column-major, or if the matrix
     * does not own its memory block, an assertion failure occurs. \a num_cols_to_append must be non-negative.
     *
     * @warning This is very slow if you repeatedly add small numbers of columns to a large matrix. Instead, preallocate memory
     *   for the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to
     *   resize().
     *
     * @see appendColumn()
     */
    void appendColumns(long num_cols_to_append)
    {
      alwaysAssertM(L == MatrixLayout::COLUMN_MAJOR, "Matrix: Cannot append column(s) to a matrix that is not column-major");
      alwaysAssertM(owns_memory, "Matrix: Cannot append column(s) to a matrix that does not own its memory block");
      alwaysAssertM(num_cols >= 0, "Matrix: Cannot append a row to a matrix with no columns");
      alwaysAssertM(num_cols_to_append >= 0, "Matrix: Cannot append a negative number of columns to a matrix");

      if (num_cols_to_append <= 0)
        return;

      long new_num_cols = num_cols + num_cols_to_append;
      long new_num_elems = num_rows * new_num_cols;
      T * new_values = allocator.allocate((size_t)new_num_elems);

      if (!this->isEmpty())
      {
        size_t old_num_elems = (size_t)this->numElements();
        Algorithms::fastCopy(values, values + old_num_elems, new_values);
        allocator.deallocate(values, old_num_elems);
      }

      values = new_values;
      num_cols = new_num_cols;
    }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get(). */
    T const & operator()(long row, long col) const
    {
      // Compiler should optimize to remove the conditional
      if (L == MatrixLayout::ROW_MAJOR)
        return values[row * num_cols + col];
      else
        return values[col * num_rows + row];
    }

    /** Element access. Use this whenever possible to avoid the virtual function overhead of get() or set(). */
    T & operator()(long row, long col)
    {
      // Compiler should optimize to remove the conditional
      if (L == MatrixLayout::ROW_MAJOR)
        return values[row * num_cols + col];
      else
        return values[col * num_rows + row];
    }

    /**
     * Element access on the unrolled matrix (by columns if column-major, or rows if row-major). This is the same as
     * data()[index].
     */
    T const & operator()(long index) const { return values[index]; }

    /**
     * Element access on the unrolled matrix (by columns if column-major, or rows if row-major). This is the same as
     * data()[index].
     */
    T & operator()(long index) { return values[index]; }

    /**
     * Element access on the unrolled matrix (by columns if column-major, or rows if row-major). This is the same as
     * data()[index].
     */
    T const & operator[](long index) const { return values[index]; }

    /**
     * Element access on the unrolled matrix (by columns if column-major, or rows if row-major). This is the same as
     * data()[index].
     */
    T & operator[](long index) { return values[index]; }

    T const & get(long row, long col) const { return (*this)(row, col); }
    T & getMutable(long row, long col) { return (*this)(row, col); }
    void set(long row, long col, T const & value) { (*this)(row, col) = value; }

    void getRow(long row, T * values) const
    {
      for (long c = 0; c < num_cols; ++c)
        values[c] = (*this)(row, c);
    }

    void setRow(long row, T const * values)
    {
      for (long c = 0; c < num_cols; ++c)
        (*this)(row, c) = values[c];
    }

    void getColumn(long col, T * values) const
    {
      for (long r = 0; r < num_rows; ++r)
        values[r] = (*this)(r, col);
    }

    void setColumn(long col, T const * values)
    {
      for (long r = 0; r < num_rows; ++r)
        (*this)(r, col) = values[r];
    }

    // This is needed to avoid a Visual Studio warning about inheriting via dominance
    void makeZero() { AddressableBaseT::makeZero(); }

    void fill(T const & value)
    {
      if (values)
      {
        for (T * vi = values, * e = values + this->numElements(); vi != e; ++vi)
          *vi = value;
      }
    }

    /** Get the storage block of the matrix. */
    T const * data() const { return values; }

    /** Get the storage block of the matrix. */
    T * data() { return values; }

    /**
     * Get the inverse of the matrix (assertion failure if matrix is not square). All computations are done using type T, so do
     * <b>not</b> call this function on integer matrices (built-in (POD) integer types will generate assertion failures).
     * Creates a new matrix for the result, so you might prefer in-place inversion using invert().
     */
    Matrix inverse() const { Matrix result = *this; result.invert(); return result; }

    /**
     * Invert the matrix (assertion failure if matrix is not square). All computations are done using type T, so do <b>not</b>
     * call this function on integer matrices (built-in (POD) integer types will generate assertion failures).
     */
    void invert()
    {
      alwaysAssertM(this->isSquare(), "Matrix: Only square matrices can be inverted");

      TheaArray<long> col_index((size_t)numRows()), row_index((size_t)numRows()), pivot((size_t)numRows());
      Internal::invertMatrix(*this, &col_index[0], &row_index[0], &pivot[0]);
    }

    /** Get the transpose of the matrix. Creates a new matrix for the result, so you might prefer transpose(Matrix &). */
    Matrix transpose() const
    {
      Matrix t(numColumns(), numRows());

      for (long r = 0; r < numRows(); ++r)
        for (long c = 0; c < numColumns(); ++c)
          t(c, r) = (*this)(r, c);

      return t;
    }

    /**
     * Get the transpose of the matrix. No memory allocation/deallocation is required if the input and output matrices have the
     * same number of elements.
     *
     * @note The result matrix cannot be the same object as this matrix.
     */
    template <typename U, MatrixLayout::Value ArgLayout> void transpose(Matrix<U, ArgLayout> & result) const
    {
      result.resize(numColumns(), numRows());

      if (L == MatrixLayout::ROW_MAJOR)
      {
        for (long r = 0; r < numRows(); ++r)
          for (long c = 0; c < numColumns(); ++c)
            result(c, r) = (*this)(r, c);
      }
      else
      {
        for (long c = 0; c < numColumns(); ++c)
          for (long r = 0; r < numRows(); ++r)
            result(c, r) = (*this)(r, c);
      }
    }

    /** Get the negation of this matrix (the matrix equal to this one except with every element negated). */
    Matrix operator-() const
    {
      Matrix neg(numRows(), numColumns());
      T * v2 = neg.values;
      for (T const * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
        *v2 = -(*v1);
    }

    /** Addition (per-component). */
    Matrix operator+(Matrix const & rhs) const
    {
      alwaysAssertM(numRows() == rhs.numRows() && numColumns() == rhs.numColumns(),
                    "Matrix: Can't add matrices of different dimensions");

      Matrix result(numRows(), numColumns());
      T * u = result.values;
      for (T const * v1 = values, * v2 = rhs.values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2, ++u)
        *u = *v1 + *v2;

      return result;
    }

    /** Add-and-assign. */
    Matrix & operator+=(Matrix const & rhs)
    {
      alwaysAssertM(numRows() == rhs.numRows() && numColumns() == rhs.numColumns(),
                    "Matrix: Can't add matrices of different dimensions");

      T const * v2 = rhs.values;
      for (T * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
        *v1 += *v2;

      return *this;
    }

    /** Subtraction (per-component). */
    Matrix operator-(Matrix const & rhs) const
    {
      alwaysAssertM(numRows() == rhs.numRows() && numColumns() == rhs.numColumns(),
                    "Matrix: Can't subtract matrices of different dimensions");

      Matrix result(numRows(), numColumns());
      T * u = result.values;
      for (T const * v1 = values, * v2 = rhs.values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2, ++u)
        *u = *v1 - *v2;

      return result;
    }

    /** Subtract-and-assign. */
    Matrix & operator-=(Matrix const & rhs)
    {
      alwaysAssertM(numRows() == rhs.numRows() && numColumns() == rhs.numColumns(),
                    "Matrix: Can't add matrices of different dimensions");

      T const * v2 = rhs.values;
      for (T * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
        *v1 -= *v2;

      return *this;
    }

    /** Post-multiply by a scalar. */
    Matrix operator*(T const & s) const
    {
      Matrix result(numRows(), numColumns());
      T * u = result.values;
      for (T const * v = values, * e = values + this->numElements(); v != e; ++v, ++u)
        *u = s * (*v);

      return result;
    }

    /** Post-multiply by a scalar and assign. */
    Matrix & operator*=(T const & s)
    {
      for (T * v = values, * e = values + this->numElements(); v != e; ++v)
        (*v) *= s;

      return *this;
    }

    /** Multiplication. */
    Matrix operator*(Matrix const & rhs) const
    {
      alwaysAssertM(numColumns() == rhs.numRows(), "Matrix: Matrices don't have compatible dimensions for multiplication");

      Matrix result(numRows(), rhs.numColumns());
      for (long i = 0; i < numRows(); ++i)
        for (long k = 0; k < rhs.numColumns(); ++k)
        {
          result(i, k) = (*this)(i, 0) * rhs(0, k);

          for (long j = 1; j < numColumns(); ++j)
            result(i, k) += (*this)(i, j) * rhs(j, k);
        }

      return result;
    }

    /**
     * Multiply a vector \a v by this matrix (M), yielding the product \a result = M * \a v.
     *
     * @param v The vector to be multiplied, must contain numColumns() elements.
     * @param result The result vector, must be preallocated to numRows() elements.
     */
    template <typename U> void postmulVector(U const * v, U * result) const
    {
      debugAssertM(v, "Matrix: Vector to be multiplied cannot be null");
      debugAssertM(result, "Matrix: Result vector of multiplication cannot be null");

      if (L == MatrixLayout::ROW_MAJOR)
      {
        for (long r = 0; r < num_rows; ++r)
        {
          result[r] = static_cast<U>(0);
          for (long c = 0; c < num_cols; ++c)
            result[r] += static_cast<U>((*this)(r, c) * v[c]);
        }
      }
      else
      {
        for (long r = 0; r < num_rows; ++r)
          result[r] = static_cast<U>(0);

        for (long c = 0; c < num_cols; ++c)
          for (long r = 0; r < num_rows; ++r)
            result[r] += static_cast<U>((*this)(r, c) * v[c]);
      }
    }

    /**
     * Multiply a vector \a v by this matrix (M), yielding the product \a result = \a v * M.
     *
     * @param v The vector to be multiplied, must contain numRows() elements.
     * @param result The result vector, must be preallocated to numColumns() elements.
     */
    template <typename U> void premulVector(U const * v, U * result) const
    {
      debugAssertM(v, "Matrix: Vector to be multiplied cannot be null");
      debugAssertM(result, "Matrix: Result vector of multiplication cannot be null");

      if (L == MatrixLayout::ROW_MAJOR)
      {
        for (long c = 0; c < num_cols; ++c)
          result[c] = static_cast<U>(0);

        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            result[c] += static_cast<U>(v[r] * (*this)(r, c));
      }
      else
      {
        for (long c = 0; c < num_cols; ++c)
        {
          result[c] = static_cast<U>(0);
          for (long r = 0; r < num_rows; ++r)
            result[c] += static_cast<U>(v[r] * (*this)(r, c));
        }
      }
    }

    /**
     * Utility function for efficiently computing  a matrix-vector product (w = M * v). This lets this matrix be passed directly
     * to ARPACK++. The implementation simply calls postmulVector().
     *
     * @see postmulVector()
     */
    template <typename U> void MultMv(U const * v, U * w) const { postmulVector(v, w); }

    /** Post-divide by a scalar. */
    Matrix operator/(T const & s) const
    {
      Matrix result(numRows(), numColumns());
      T * u = result.values;
      for (T const * v = values, * e = values + this->numElements(); v != e; ++v, ++u)
        *u = s / (*v);

      return result;
    }

    /** Post-divide by a scalar and assign. */
    Matrix & operator/=(T const & s)
    {
      for (T * v = values, * e = values + this->numElements(); v != e; ++v)
        (*v) /= s;

      return *this;
    }

    T const & min() const { return AddressableBaseT::min(); }
    T const & max() const { return AddressableBaseT::max(); }

    /** Return a matrix containing the component-wise minima of this matrix and another. */
    Matrix min(Matrix const & other) const
    {
      Matrix result(num_rows, num_cols);
      for (long r = 0; r < num_rows; ++r)
        for (long c = 0; c < num_cols; ++c)
          result(r, c) = std::min((*this)(r, c), other(r, c));

      return result;
    }

    /** Return a matrix containing the component-wise maxima of this matrix and another. */
    Matrix max(Matrix const & other) const
    {
      Matrix result(num_rows, num_cols);
      for (long r = 0; r < num_rows; ++r)
        for (long c = 0; c < num_cols; ++c)
          result(r, c) = std::max((*this)(r, c), other(r, c));

      return result;
    }

    /** Per-element sign (-1, 0 or 1). */
    Matrix sign() const
    {
      Matrix result(num_rows, num_cols);
      for (long r = 0; r < num_rows; ++r)
        for (long c = 0; c < num_cols; ++c)
          result(r, c) = Math::sign((*this)(r, c));

      return result;
    }

  private:
    long num_rows, num_cols;
    T * values;
    bool owns_memory;

}; // class Matrix

/** Pre-multiply by a scalar. */
template <typename T, MatrixLayout::Value L>
Matrix<T, L>
operator*(T const & s, Matrix<T, L> const & m)
{
  return m * s;
}

} // namespace Thea

#endif
