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
#include <boost/type_traits/conditional.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <algorithm>
#include <memory>

namespace Thea {

// Forward declarations
template <typename T, MatrixLayout::Value L, typename Index2DT, typename Index1DT> class CompressedSparseMatrix;
template <typename T, typename Index2DT, typename Index1DT> class CompressedRowMatrix;
template <typename T, typename Index2DT, typename Index1DT> class CompressedColumnMatrix;
template <typename T, MatrixLayout::Value L, typename AllocT> class Matrix;
template <typename T, MatrixLayout::Value L, typename AllocT> class Vector;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for dynamic-size matrices and vectors, in either row or column-major format.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <typename T, MatrixLayout::Value L, bool IsVector, typename AllocT>
class /* THEA_DLL_LOCAL */ MatrixBase : public AddressableMatrix<T>, public ResizableMatrix<T>
{
  protected:
    typedef AddressableMatrix<T>  AddressableBaseT;
    typedef ResizableMatrix<T>    ResizableBaseT;

  public:
    static MatrixLayout::Value const Layout = L;  ///< The layout of the matrix (row-major or column-major).

    /** The opposite (transposed) matrix layout. */
    static MatrixLayout::Value const TransposedLayout = (L == MatrixLayout::ROW_MAJOR ? MatrixLayout::COLUMN_MAJOR
                                                                                      : MatrixLayout::ROW_MAJOR);

    /** The user-facing matrix type. */
    typedef typename boost::conditional< IsVector, Vector<T, L, AllocT>, Matrix<T, L, AllocT> >::type MatrixT;

    /** Multiplication result type. */
    template <MatrixLayout::Value L1, bool V1, MatrixLayout::Value L2, bool V2>
    struct MultResultT
    {
      typedef typename boost::conditional< V1,
                  typename boost::conditional< V2,
                      Matrix<T, L1, AllocT>,                                             // V * V = M
                      typename boost::conditional< L1 == MatrixLayout::ROW_MAJOR,
                          Vector<T, MatrixLayout::ROW_MAJOR, AllocT>,                    // RV * M = RV
                          Matrix<T, L1, AllocT> >::type                                  // CV * M = M
                      >::type,
                  typename boost::conditional< V2,
                      typename boost::conditional< L2 == MatrixLayout::COLUMN_MAJOR,
                          Vector<T, MatrixLayout::COLUMN_MAJOR, AllocT>,                 // M * CV = CV
                          Matrix<T, L1, AllocT> >::type,                                 // M * RV = M
                      Matrix<T, L1, AllocT>                                              // M * M = M
                      >::type
              >::type
          type;
    };

    typedef T                                 Value;                 ///< Type of values stored in the matrix.
    typedef T                                 value_type;            ///< Type of values stored in the matrix (STL convention).
    typedef T *                               Iterator;              ///< Forward iterator through elements.
    typedef T const *                         ConstIterator;         ///< Forward const iterator through elements.
    typedef std::reverse_iterator<T *>        ReverseIterator;       ///< Reverse iterator through elements.
    typedef std::reverse_iterator<T const *>  ConstReverseIterator;  ///< Reverse const iterator through elements.

    THEA_DEF_POINTER_TYPES(MatrixT, shared_ptr, weak_ptr)

  protected:
    /** Default constructor. */
    MatrixBase()
    : num_rows(Layout == MatrixLayout::ROW_MAJOR ? 1 : 0),
      num_cols(Layout == MatrixLayout::ROW_MAJOR ? 0 : 1),
      values(NULL), owns_memory(true)
    {}

    /** Copy constructor. */
    MatrixBase(MatrixBase const & m) : num_rows(0), num_cols(0), values(NULL), owns_memory(true) { *this = m; }

    /** Copy from any compatible template instantiation. */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixBase(MatrixBase<T, L2, V2, A2> const & m)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    { *this = m; }

    /** Constructs a matrix of a given size, with uninitialized values. */
    MatrixBase(long num_rows_, long num_cols_) : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      resize(num_rows_, num_cols_);
    }

    /** Constructs a matrix of a given size, filled with a given value. */
    MatrixBase(long num_rows_, long num_cols_, T const & fill_value) : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      resize(num_rows_, num_cols_, fill_value);
    }

    /**
     * Constructs a row vector or a column vector, depending on whether the matrix layout is MatrixLayout::ROW_MAJOR or
     * MatrixLayout::COLUMN_MAJOR, with uninitialized values.
     */
    MatrixBase(long vector_size) : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      resize(vector_size);
    }

    /**
     * Constructs this matrix as a wrapper for an existing block of storage. The matrix thus created is <b>not resizable</b>,
     * and the memory block will <b>not be freed</b> when the Matrix object is destroyed.
     */
    MatrixBase(T * existing_data, long num_rows_, long num_cols_)
    : num_rows(num_rows_), num_cols(num_cols_), values(existing_data), owns_memory(false)
    {
      alwaysAssertM(validDims(num_rows_, num_cols_), "MatrixBase: Invalid dimensions");
    }

    /**
     * Constructs a vector that wraps an existing block of storage. The vector thus created is <b>not resizable</b>, and the
     * memory block will <b>not be freed</b> when the Matrix object is destroyed. A row vector or a column vector is created,
     * depending on whether the matrix layout is MatrixLayout::ROW_MAJOR or MatrixLayout::COLUMN_MAJOR.
     */
    MatrixBase(T * existing_data, long vector_size)
    : num_rows(Layout == MatrixLayout::ROW_MAJOR ? 1 : vector_size),
      num_cols(Layout == MatrixLayout::ROW_MAJOR ? vector_size : 1),
      values(existing_data), owns_memory(false)
    {
      alwaysAssertM(vector_size >= 0, "MatrixBase: Vector length must be non-negative");
    }

    /** Initialize from an addressable matrix. */
    template <typename M> explicit
    MatrixBase(M const & src, typename boost::enable_if< boost::is_base_of< AddressableMatrix<typename M::Value>,
                                                                            M > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(src.numRows(), src.numColumns());
      else
      {
        alwaysAssertM(num_rows == src.numRows() && num_cols == src.numColumns(),  // essentially, 0x0
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
    template <typename M> explicit
    MatrixBase(M const & crm,
               typename boost::enable_if< boost::is_base_of< CompressedRowMatrix< typename M::Value,
                                                                                  typename M::Index2D,
                                                                                  typename M::Index1D >,
                                                             M > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(crm.numRows(), crm.numColumns());
      else
      {
        alwaysAssertM(num_rows == crm.numRows() && num_cols == crm.numColumns(),  // essentially, 0x0
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
    template <typename M> explicit
    MatrixBase(M const & ccm,
               typename boost::enable_if< boost::is_base_of< CompressedColumnMatrix< typename M::Value,
                                                                                     typename M::Index2D,
                                                                                     typename M::Index1D >,
                                                             M > >::type * dummy = NULL)
    : num_rows(0), num_cols(0), values(NULL), owns_memory(true)
    {
      if (owns_memory)
        resize(ccm.numRows(), ccm.numColumns());
      else
      {
        alwaysAssertM(num_rows == ccm.numRows() && num_cols == ccm.numColumns(),  // essentially, 0x0
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

  public:
    /** Destructor. */
    ~MatrixBase()
    {
      if (owns_memory) allocator.deallocate(values, (size_t)this->numElements());
    }

    /** Copy assignment operator. */
    MatrixBase & operator=(MatrixBase const & src)
    {
      if (owns_memory)
        resize(src.num_rows, src.num_cols);
      else
        alwaysAssertM(sameDims(src), "MatrixBase: A wrapper matrix cannot be assigned a value of different dimensions");

      Algorithms::fastCopy(src.values, src.values + this->numElements(), values);
      return *this;
    }

    /** Assign from any compatible template instantiation. */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixBase & operator=(MatrixBase<T, L2, V2, A2> const & src)
    {
      if (owns_memory)
        resize(src.numRows(), src.numColumns());
      else
        alwaysAssertM(sameDims(src), "MatrixBase: A wrapper matrix cannot be assigned a value of different dimensions");

      if (L == L2)
        Algorithms::fastCopy(src.data(), src.data() + this->numElements(), values);
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            (*this)(r, c) = src(r, c);
      }

      return *this;
    }

    /** Get the layout of the matrix (row or column major) */
    static MatrixLayout getLayout() { return L; }

    long numRows() const { return num_rows; }
    long numColumns() const { return num_cols; }

    void resize(long num_rows_, long num_cols_)
    {
      if (num_rows != num_rows_ || num_cols != num_cols_)
      {
        alwaysAssertM(owns_memory, "MatrixBase: This matrix does not own its storage block and can't be resized");
        alwaysAssertM(validDims(num_rows_, num_cols_), "MatrixBase: Invalid dimensions");

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
     * Resize a vector. The result is a row vector or a column vector depending on whether the matrix layout is
     * MatrixLayout::ROW_MAJOR or MatrixLayout::COLUMN_MAJOR.
     */
    void resize(long vector_size)
    {
      if (Layout == MatrixLayout::ROW_MAJOR)
        resize(1, vector_size);
      else
        resize(vector_size, 1);
    }

    /** Reset the matrix to zero elements. */
    void clear()
    {
      if (IsVector)
        resize(0);
      else
        resize(0, 0);
    }

    /**
     * Append a single (uninitialized) row to a row-major matrix. If the matrix is not row-major (or a column vector), or if the
     * matrix does not own its memory block, an assertion failure occurs.
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
     * Append one or more (uninitialized) rows to a row-major matrix. If the matrix is not row-major (or a column vector), or if
     * the matrix does not own its memory block, an assertion failure occurs. \a num_rows_to_append must be non-negative.
     *
     * @warning This is very slow if you repeatedly add small numbers of rows to a large matrix. Instead, preallocate memory for
     *   the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to resize().
     *
     * @see appendRow()
     */
    void appendRows(long num_rows_to_append)
    {
      alwaysAssertM((L == MatrixLayout::ROW_MAJOR && !IsVector) || (L == MatrixLayout::COLUMN_MAJOR && IsVector),
                    "MatrixBase: Cannot append row(s) to a matrix that is not row-major");
      alwaysAssertM(owns_memory, "MatrixBase: Cannot append row(s) to a matrix that does not own its memory block");
      alwaysAssertM(num_rows_to_append >= 0, "MatrixBase: Cannot append a negative number of rows to a matrix");

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
     * Append a single (uninitialized) column to a column-major matrix. If the matrix is not column-major (or a row vector), or
     * if the matrix does not own its memory block, an assertion failure occurs.
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
     * Append one or more (uninitialized) columns to a column-major matrix. If the matrix is not column-major (or a row vector),
     * or if the matrix does not own its memory block, an assertion failure occurs. \a num_cols_to_append must be non-negative.
     *
     * @warning This is very slow if you repeatedly add small numbers of columns to a large matrix. Instead, preallocate memory
     *   for the matrix if possible. It is better to overallocate and then free unused memory at the end with a call to
     *   resize().
     *
     * @see appendColumn()
     */
    void appendColumns(long num_cols_to_append)
    {
      alwaysAssertM((L == MatrixLayout::COLUMN_MAJOR && !IsVector) || (L == MatrixLayout::ROW_MAJOR && IsVector),
                    "MatrixBase: Cannot append column(s) to a matrix that is not column-major");
      alwaysAssertM(owns_memory, "MatrixBase: Cannot append column(s) to a matrix that does not own its memory block");
      alwaysAssertM(num_cols >= 0, "MatrixBase: Cannot append a row to a matrix with no columns");
      alwaysAssertM(num_cols_to_append >= 0, "MatrixBase: Cannot append a negative number of columns to a matrix");

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

    /** Get an iterator pointing to the first element of the matrix (in storage order). */
    Iterator begin() { return values; }

    /** Get a const iterator pointing to the first element of the matrix (in storage order). */
    ConstIterator begin() const { return values; }

    /** Get an iterator pointing to the end of the matrix (in storage order). */
    Iterator end() { return values + num_rows * num_cols; }

    /** Get a const iterator pointing to the end of the matrix (in storage order). */
    ConstIterator end() const { return values + num_rows * num_cols; }

    /** Get a reverse iterator pointing to the last element of the matrix (in storage order). */
    ReverseIterator rbegin() { return ReverseIterator(end()); }

    /** Get a const reverse iterator pointing to the last element of the matrix (in storage order). */
    ConstReverseIterator rbegin() const { return ConstReverseIterator(end()); }

    /** Get a reverse iterator pointing to the beginning of the matrix (in storage order). */
    ReverseIterator rend() { return ReverseIterator(begin()); }

    /** Get a const reverse iterator pointing to the beginning of the matrix (in storage order). */
    ConstReverseIterator rend() const { return ConstReverseIterator(begin()); }

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
    MatrixT inverse() const { MatrixT result = *this; result.invert(); return result; }

    /**
     * Invert the matrix (assertion failure if matrix is not square). All computations are done using type T, so do <b>not</b>
     * call this function on integer matrices (built-in (POD) integer types will generate assertion failures).
     */
    void invert()
    {
      alwaysAssertM(this->isSquare(), "MatrixBase: Only square matrices can be inverted");

      TheaArray<long> col_index((size_t)num_rows), row_index((size_t)num_rows), pivot((size_t)num_rows);
      Internal::invertMatrix(*this, &col_index[0], &row_index[0], &pivot[0]);
    }

    /**
     * Get the transpose of the matrix. Creates a new matrix for the result, so you might prefer transpose(MatrixBase &).
     *
     * @note Unlike transposeLayout(), this will not swap the layout of the matrix (row-major to column-major or vice versa) if
     *   the matrix is not a vector. However, the swap is necessary for vectors.
     */
    typename boost::conditional<IsVector, Vector<T, TransposedLayout, AllocT>, MatrixT>::type transpose() const
    {
      typename boost::conditional<IsVector, Vector<T, TransposedLayout, AllocT>, MatrixT>::type t(num_cols, num_rows);

      if (IsVector)
        Algorithms::fastCopy(values, values + this->numElements(), t.data());
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            t(c, r) = (*this)(r, c);
      }

      return t;
    }

    /**
     * Get the transpose of this matrix as a matrix with the opposite layout. Requires a memory allocation (for the new matrix),
     * but the copying is cheap since the underlying elements are stored in the same order and can be copied as a block.
     *
     * @note Unlike transpose(), this <i>always</i> swaps the layout of the matrix (row-major to column-major or vice versa).
     */
    typename boost::conditional< IsVector, Vector<T, TransposedLayout, AllocT>,
                                           Matrix<T, TransposedLayout, AllocT> >::type
    transposeLayout() const
    {
      typename boost::conditional< IsVector, Vector<T, TransposedLayout, AllocT>,
                                             Matrix<T, TransposedLayout, AllocT> >::type t(num_cols, num_rows);
      Algorithms::fastCopy(values, values + this->numElements(), t.data());

      return t;
    }

    /**
     * Get the transpose of the matrix. No memory allocation/deallocation is required if the input and output matrices have the
     * same number of elements.
     *
     * @note The result matrix cannot be the same object as this matrix.
     */
    template <typename U, MatrixLayout::Value L2, bool V2, typename A2> void transpose(MatrixBase<U, L2, V2, A2> & result) const
    {
      result.resize(num_cols, num_rows);
      if (L != L2)  // element storage order remains same
      {
        U * u = result.data();
        for (T const * v = values, * e = values + this->numElements(); v != e; ++v, ++u)
          *u = static_cast<U>(v);
      }
      else if (L == MatrixLayout::ROW_MAJOR)
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            result(c, r) = (*this)(r, c);
      }
      else
      {
        for (long c = 0; c < num_cols; ++c)
          for (long r = 0; r < num_rows; ++r)
            result(c, r) = (*this)(r, c);
      }
    }

    /** Get the negation of this matrix (the matrix equal to this one except with every element negated). */
    MatrixT operator-() const
    {
      MatrixT neg; neg.resize(num_rows, num_cols);
      T * v2 = neg.values;
      for (T const * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
        *v2 = -(*v1);
    }

    /** Addition (per-component). */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixT operator+(MatrixBase<T, L2, V2, A2> const & rhs) const
    {
      alwaysAssertM(sameDims(rhs), "MatrixBase: Can't add matrices of different dimensions");

      MatrixT result; result.resize(num_rows, num_cols);
      if (L == L2)
      {
        T * u = result.values;
        for (T const * v1 = values, * v2 = rhs.data(), * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2, ++u)
          *u = *v1 + *v2;
      }
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            result(r, c) = (*this)(r, c) + rhs(r, c);
      }

      return result;
    }

    /** Add-and-assign. */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixT & operator+=(MatrixBase<T, L2, V2, A2> const & rhs)
    {
      alwaysAssertM(sameDims(rhs), "MatrixBase: Can't add matrices of different dimensions");

      if (L == L2)
      {
        T const * v2 = rhs.data();
        for (T * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
          *v1 += *v2;
      }
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            (*this)(r, c) += rhs(r, c);
      }

      return *static_cast<MatrixT *>(this);
    }

    /** Subtraction (per-component). */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixT operator-(MatrixBase<T, L2, V2, A2> const & rhs) const
    {
      alwaysAssertM(sameDims(rhs), "MatrixBase: Can't subtract matrices of different dimensions");

      MatrixT result; result.resize(num_rows, num_cols);
      if (L == L2)
      {
        T * u = result.values;
        for (T const * v1 = values, * v2 = rhs.data(), * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2, ++u)
          *u = *v1 - *v2;
      }
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            result(r, c) = (*this)(r, c) - rhs(r, c);
      }

      return result;
    }

    /** Subtract-and-assign. */
    template <MatrixLayout::Value L2, bool V2, typename A2> MatrixT & operator-=(MatrixBase<T, L2, V2, A2> const & rhs)
    {
      alwaysAssertM(sameDims(rhs), "MatrixBase: Can't subtract matrices of different dimensions");

      if (L == L2)
      {
        T const * v2 = rhs.data();
        for (T * v1 = values, * e1 = values + this->numElements(); v1 != e1; ++v1, ++v2)
          *v1 -= *v2;
      }
      else
      {
        for (long r = 0; r < num_rows; ++r)
          for (long c = 0; c < num_cols; ++c)
            (*this)(r, c) -= rhs(r, c);
      }

      return *static_cast<MatrixT *>(this);
    }

    /** Post-multiply by a scalar. */
    template <typename S> typename boost::enable_if< IsCompatibleScalar<S, T>, MatrixT >::type operator*(S const & s) const
    {
      MatrixT result; result.resize(num_rows, num_cols);
      T * u = result.values;
      for (T const * v = values, * e = values + this->numElements(); v != e; ++v, ++u)
        *u = static_cast<T>((*v) * s);

      return result;
    }

    /** Post-multiply by a scalar and assign. */
    template <typename S> typename boost::enable_if< IsCompatibleScalar<S, T>, MatrixT >::type & operator*=(S const & s)
    {
      for (T * v = values, * e = values + this->numElements(); v != e; ++v)
        *v = static_cast<T>((*v) * s);

      return *static_cast<MatrixT *>(this);
    }

    /** Multiplication. */
    template <MatrixLayout::Value L2, bool V2, typename A2>
    typename MultResultT<L, IsVector, L2, V2>::type operator*(MatrixBase<T, L2, V2, A2> const & rhs) const
    {
      alwaysAssertM(num_cols == rhs.numRows(), "MatrixBase: Matrices don't have compatible dimensions for multiplication");

      long rhs_ncols = rhs.numColumns();
      typename MultResultT<L, IsVector, L2, V2>::type result;
      result.resize(num_rows, rhs_ncols);  // do this in a separate step to ensure the correct constructor is called above

      for (long i = 0; i < num_rows; ++i)
        for (long k = 0; k < rhs_ncols; ++k)
        {
          result(i, k) = (*this)(i, 0) * rhs(0, k);

          for (long j = 1; j < num_cols; ++j)
            result(i, k) += (*this)(i, j) * rhs(j, k);
        }

      return result;
    }

    /**
     * Multiply a vector \a v by this matrix (M), yielding the product \a result = M * \a v.
     *
     * @param v The vector to be multiplied, must contain num_cols elements.
     * @param result The result vector, must be preallocated to num_rows elements.
     */
    template <typename U> void postmulVector(U const * v, U * result) const
    {
      debugAssertM(v, "MatrixBase: Vector to be multiplied cannot be null");
      debugAssertM(result, "MatrixBase: Result vector of multiplication cannot be null");

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
     * @param v The vector to be multiplied, must contain num_rows elements.
     * @param result The result vector, must be preallocated to num_cols elements.
     */
    template <typename U> void premulVector(U const * v, U * result) const
    {
      debugAssertM(v, "MatrixBase: Vector to be multiplied cannot be null");
      debugAssertM(result, "MatrixBase: Result vector of multiplication cannot be null");

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
    template <typename S> typename boost::enable_if< IsCompatibleScalar<S, T>, MatrixT >::type operator/(S const & s) const
    {
      MatrixT result; result.resize(num_rows, num_cols);
      T * u = result.values;
      for (T const * v = values, * e = values + this->numElements(); v != e; ++v, ++u)
        *u = static_cast<T>((*v) / s);

      return result;
    }

    /** Post-divide by a scalar and assign. */
    template <typename S> typename boost::enable_if< IsCompatibleScalar<S, T>, MatrixT >::type & operator/=(S const & s)
    {
      for (T * v = values, * e = values + this->numElements(); v != e; ++v)
        *v = static_cast<T>((*v) / s);

      return *static_cast<MatrixT *>(this);
    }

  protected:
    AllocT allocator;  // not static because: https://bytes.com/topic/c/answers/131402-could-should-allocator-t-static-member
    long num_rows, num_cols;
    T * values;
    bool owns_memory;

    /** Check if proposed dimensions are valid for this matrix. */
    static bool validDims(long num_rows_, long num_cols_)
    {
      return (IsVector && ((Layout == MatrixLayout::ROW_MAJOR    && num_rows_ == 1 && num_cols_ >= 0)
                        || (Layout == MatrixLayout::COLUMN_MAJOR && num_cols_ == 1 && num_rows_ >= 0)))
          || (!IsVector && num_rows_ >= 0 && num_cols_ >= 0);
    }

    /** Check if two matrices have the same dimensions. */
    template <typename U, MatrixLayout::Value L2, bool V2, typename A2> bool sameDims(MatrixBase<U, L2, V2, A2> const & m) const
    { return num_rows == m.numRows() && num_cols == m.numColumns(); }

}; // class MatrixBase

} // namespace Internal

/** Pre-multiply by a scalar. */
template <typename S, typename T, MatrixLayout::Value L, bool V, typename A>
typename boost::enable_if< IsCompatibleScalar<S, T>, typename Internal::MatrixBase<T, L, V, A>::MatrixT >::type
operator*(S const & s, Internal::MatrixBase<T, L, V, A> const & m)
{
  return m * s;
}

/** A resizable dense matrix, either row-major or column-major. */
template < typename T, MatrixLayout::Value L = MatrixLayout::ROW_MAJOR, typename AllocT = std::allocator<T> >
class /* THEA_API */ Matrix : public Internal::MatrixBase<T, L, false, AllocT>
{
  private:
     typedef Internal::MatrixBase<T, L, false, AllocT> BaseT;

  public:
    /** Constructs an empty matrix. */
    Matrix() {}

    /** Constructs a matrix of a given size, with uninitialized values. */
    Matrix(long num_rows_, long num_cols_) : BaseT(num_rows_, num_cols_) {}

    /** Constructs a matrix of a given size, filled with a given value. */
    Matrix(long num_rows_, long num_cols_, T const & fill_value) : BaseT(num_rows_, num_cols_, fill_value) {}

    /**
     * Constructs this matrix as a wrapper for an existing block of storage. The matrix thus created is <b>not resizable</b>,
     * and the memory block will <b>not be freed</b> when the Matrix object is destroyed.
     */
    Matrix(T * existing_data, long num_rows_, long num_cols_) : BaseT(existing_data, num_rows_, num_cols_) {}

    /** Construct a square diagonal matrix, given the values on the diagonal. */
    template <typename U> static Matrix fromDiagonal(long n, U const * diagonal)
    {
      Matrix m(n, n, static_cast<T>(0));
      for (long i = 0; i < n; ++i)
        m(i, i) = static_cast<T>(diagonal[i]);

      return m;
    }

    /** Construct a square diagonal matrix, given the values on the diagonal as an array object. */
    template <typename U> static Matrix fromDiagonal(TheaArray<U> const & diagonal)
    {
      return fromDiagonal((long)diagonal.size(), !diagonal.empty() ? &diagonal[0] : NULL);
    }

    /** Construct a square diagonal matrix, given the values on the diagonal as a (row or column) vector. */
    template <typename U, MatrixLayout::Value L2, bool V2, typename A2>
    static Matrix fromDiagonal(Internal::MatrixBase<U, L2, V2, A2> const & diagonal)
    {
      alwaysAssertM(diagonal.isVector(), "Matrix: Diagonal values must be supplied as a vector");
      return fromDiagonal(diagonal.numElements(), diagonal.data());
    }

    /** Initialize from an addressable matrix. */
    template <typename M> explicit
    Matrix(M const & src, typename boost::enable_if< boost::is_base_of< AddressableMatrix<typename M::Value>,
                                                                        M > >::type * dummy = NULL)
    : BaseT(src) {}

    /**
     * Initialize from a compressed row matrix. The class T must be <code>static_cast</code>able from an integer, specifically
     * from zero.
     *
     * <b>Implementation note:</b> Since this is a templated function, the forward declaration is enough. The body of
     * CompressedRowMatrix will necessarily be available at the point of instantiation.
     */
    template <typename M> explicit
    Matrix(M const & crm,
           typename boost::enable_if< boost::is_base_of< CompressedRowMatrix< typename M::Value,
                                                                              typename M::Index2D,
                                                                              typename M::Index1D >,
                                                         M > >::type * dummy = NULL)
    : BaseT(crm) {}

    /**
     * Initialize from a compressed column matrix. The class T must be <code>static_cast</code>able from an integer,
     * specifically from zero.
     *
     * <b>Implementation note:</b> Since this is a templated function, the forward declaration is enough. The body of
     * CompressedColumnMatrix will necessarily be available at the point of instantiation.
     */
    template <typename M> explicit
    Matrix(M const & ccm,
           typename boost::enable_if< boost::is_base_of< CompressedColumnMatrix< typename M::Value,
                                                                                 typename M::Index2D,
                                                                                 typename M::Index1D >,
                                                         M > >::type * dummy = NULL)
    : BaseT(ccm) {}

    /** Copy constructor. */
    Matrix(Matrix const & m) : BaseT(m) {}

    /** Construct from any compatible base type. */
    template <MatrixLayout::Value L2, bool V2, typename A2> Matrix(Internal::MatrixBase<T, L2, V2, A2> const & m) : BaseT(m) {}

    /** Copy assignment operator. */
    Matrix & operator=(Matrix const & src)
    {
      BaseT::operator=(src);
      return *this;
    }

    /** Assign from any compatible base type. */
    template <MatrixLayout::Value L2, bool V2, typename A2> Matrix & operator=(Internal::MatrixBase<T, L2, V2, A2> const & src)
    {
      BaseT::operator=(src);
      return *this;
    }

}; // class Matrix

} // namespace Thea

#endif
