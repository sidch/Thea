//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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

#ifndef __Thea_CompressedSparseMatrix_hpp__
#define __Thea_CompressedSparseMatrix_hpp__

#include "Common.hpp"
#include "AddressableMatrix.hpp"
#include "Array.hpp"
#include "BasicMatrix.hpp"
#include "Map.hpp"
#include "ResizableMatrix.hpp"
#include "Algorithms/FastCopy.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace Thea {

namespace CompressedSparseMatrixInternal {

struct Base {};

} // namespace CompressedSparseMatrixInternal

/**
 * A matrix stored in the Compressed Row/Column format (http://www.cs.utk.edu/~dongarra/etemplates/node373.html). Significant
 * space savings are possible if the matrix is sparse, but the representation is not amenable to most matrix operations, except
 * multiplication which is efficient.
 *
 * This class has hidden constructors and cannot be used directly. Instead use the CompressedRowMatrix and
 * CompressedColumnMatrix derivatives.
 *
 * <b>Implementation-specific:</b> The subscripts 1 and 2 stand for row and column respectively if the matrix is interpreted as
 * being in row-major format, and column and row if the matrix is interpreted as being in column-major format. This <b>NOT</b>
 * the same as the Boost.uBLAS convention, where 1 always stands for rows and 2 for columns.
 */
template <typename T, MatrixLayout::Value L, typename Index2DT = int, typename Index1DT = long>
class /* THEA_API */ CompressedSparseMatrix : public virtual BasicMatrix<T>, public CompressedSparseMatrixInternal::Base
{
  public:
    THEA_DEF_POINTER_TYPES(CompressedSparseMatrix, shared_ptr, weak_ptr)

    static MatrixLayout::Value const Layout = L;  ///< The layout of the matrix (row-major or column-major).
    typedef Index2DT  Index2D;                    ///< The type of 2D element indices.
    typedef Index1DT  Index1D;                    ///< The type of 1D (flat) element indices.

    // Since this class can't be directly instantiated (protected constructors), we might as well make these public, to avoid
    // access issues for some templated constructors
    long                 size1;
    long                 size2;
    TheaArray<Index1D>  indices1;
    TheaArray<Index2D>  indices2;
    TheaArray<T>        values;

    /** Assignment operator. */
    CompressedSparseMatrix & operator=(CompressedSparseMatrix const & src)
    {
      size1 = src.size1;
      size2 = src.size2;
      indices1 = src.indices1;
      indices2 = src.indices2;
      values = src.values;

      return *this;
    }

    /** Get the layout of the matrix (row or column major) */
    static MatrixLayout getLayout() { return L; }

    /** Get the set of logically non-zero values in the matrix. */
    TheaArray<T> const & getValues() const { return values; }

    /** Get the set of logically non-zero values in the matrix. Assert isValid() if you modify this array. */
    TheaArray<T> & getValues() { return values; }

    /** Get the number of logically non-zero entries in the matrix. */
    long numSetElements() const { return (long)values.size(); }

    /**
     * Check that the matrix is internally consistent. You should assert this after externally modifying the underlying arrays.
     */
    bool isValid() const
    {
      return (indices2.size() == values.size()
           && indices1.size() == (array_size_t)size1 + 1
           && (long)indices1[indices1.size() - 1] == numSetElements());
    }

    void makeZero()
    {
      indices2.clear();
      values.clear();
      for (typename TheaArray<Index1D>::iterator i = indices1.begin(); i != indices1.end(); ++i)
        *i = static_cast<T>(0);
    }

  protected:
    /**
     * Constructs a zero matrix of the specified size. If the number of columns is omitted or zero, a square matrix is created.
     */
    CompressedSparseMatrix(int size1_ = 0, int size2_ = 0)
    : size1(size1_), size2(size2_ ? size2_ : size1_), indices1((array_size_t)size1_ + 1, 0)
    {
      alwaysAssertM(size1_ >= 0 && size2_ >= 0, "CompressedSparseMatrix: Dimensions must be non-negative");
    }

    /** Copy constructor. */
    CompressedSparseMatrix(CompressedSparseMatrix const & src)
    : size1(src.size1), size2(src.size2), indices1(src.indices1), indices2(src.indices2), values(src.values)
    {}

    /** Initialize from another compressed sparse matrix. */
    template <typename MatrixT> explicit
    CompressedSparseMatrix(MatrixT const & src,
                           typename boost::enable_if< boost::is_base_of< CompressedSparseMatrixInternal::Base,
                                                                         MatrixT > >::type * dummy = NULL)
    : size1(0), size2(0)
    {
      if (src.getLayout() == L)
      {
        size1 = src.size1;
        size2 = src.size2;

        indices1.resize(src.indices1.size());
        Algorithms::fastCopy(src.indices1.begin(), src.indices1.end(), indices1.begin());

        indices2.resize(src.indices2.size());
        Algorithms::fastCopy(src.indices2.begin(), src.indices2.end(), indices2.begin());

        values.resize(src.values.size());
        Algorithms::fastCopy(src.values.begin(), src.values.end(), values.begin());
      }
      else
      {
        size1 = src.size2;
        size2 = src.size1;

        indices1.resize(src.size2 + 1);
        std::fill(indices1.begin(), indices1.end(), 0);

        indices2.resize(src.indices1.size());
        values.resize(src.values.size());

        typedef std::pair<typename MatrixT::Index2D, int> IndexPair;
        typedef TheaMap<IndexPair, typename MatrixT::Value> ValueMap;

        // Put all the source values into a map, with indices swapped so they can be read off using the inverse layout
        ValueMap src_values;
        array_size_t curr_pos = 0;
        for (int i = 0; i < src.size1; ++i)
        {
          int num_elems = (int)(src.indices1[(array_size_t)i + 1] - src.indices1[(array_size_t)i]);
          for (int e = 0; e < num_elems; ++e, ++curr_pos)
          {
            typename MatrixT::Index2D j = src.indices2[curr_pos];
            src_values[IndexPair(j, i)] = src.values[curr_pos];
          }
        }

        // Now read back the values in lexicographic order
        long i, last_index1 = -1;
        curr_pos = 0;
        for (typename ValueMap::const_iterator vi = src_values.begin(); vi != src_values.end(); ++vi, ++curr_pos)
        {
          i = (long)vi->first.first;
          if (i != last_index1)  // change of primary index since last iteration
          {
            for (long k = last_index1 + 1; k <= i; ++k)
              indices1[(array_size_t)k] = curr_pos;

            last_index1 = i;
          }

          indices2[curr_pos] = static_cast<Index2D>(vi->first.second);
          values  [curr_pos] = static_cast<T>(vi->second);
        }
      }
    }

    /** Initialize from an addressable matrix. */
    template <typename MatrixT> explicit
    CompressedSparseMatrix(MatrixT const & src,
                           typename boost::enable_if< boost::is_base_of< AddressableMatrix<typename MatrixT::Value>,
                                                                         MatrixT > >::type * dummy = NULL)
    : size1(src.numRows()), size2(src.numColumns())
    {
      // All the conditionals on layout will hopefully be optimized out by the compiler

      if (L == MatrixLayout::ROW_MAJOR)
      {
        size1 = src.numRows();
        size2 = src.numColumns();
      }
      else
      {
        size1 = src.numColumns();
        size2 = src.numRows();
      }

      indices1.resize((array_size_t)size1 + 1);
      indices2.clear();
      values.clear();

      Index1D curr_pos = 0;  // current position in value array
      for (int i = 0; i < size1; ++i)
      {
        indices1[(array_size_t)i] = curr_pos;

        for (int j = 0; j < size2; ++j)
        {
          T value = static_cast<T>(L == MatrixLayout::ROW_MAJOR ? src.get(i, j) : src.get(j, i));
          if (value != 0)  // only retain non-zero values
          {
            indices2.push_back((Index1D)j);
            values.push_back(value);
            curr_pos++;
          }
        }
      }

      // One beyond last position, stores number of non-zeros
      indices1[(array_size_t)size1] = curr_pos;
    }

    /** Resizes the matrix to the specified dimensions. All existing data is discarded and the matrix is set to zero. */
    void resize(long size1_, long size2_)
    {
      alwaysAssertM(size1_ >= 0 && size2_ >= 0, "CompressedSparseMatrix: Dimensions must be non-negative");

      size1 = size1_;
      size2 = size2_;
      indices1.resize((array_size_t)size1 + 1, 0);
      indices2.clear();
      values.clear();
    }

}; // class CompressedSparseMatrix

/** Row-major sparse matrix. Non-zero values are packed contiguously row-by-row. */
template <typename T, typename Index2DT = int, typename Index1DT = int>
class /* THEA_API */ CompressedRowMatrix : public CompressedSparseMatrix<T, MatrixLayout::ROW_MAJOR, Index2DT, Index1DT>,
                                           public ResizableMatrix<T>
{
  private:
    typedef CompressedSparseMatrix<T, MatrixLayout::ROW_MAJOR, Index2DT, Index1DT> BaseT;

  public:
    THEA_DEF_POINTER_TYPES(CompressedRowMatrix, shared_ptr, weak_ptr)

    /**
     * Default constructor. Initializes a zero matrix of the specified size. If the number of columns is omitted or zero, a
     * square matrix is created.
     */
    CompressedRowMatrix(long num_rows = 0, long num_cols = 0) : BaseT(num_rows, num_cols) {}

    /** Copy contructor. */
    CompressedRowMatrix(CompressedRowMatrix const & src) : BaseT(src) {}

    /** Templated copy constructor. */
    template <typename MatrixT> explicit CompressedRowMatrix(MatrixT const & m) : BaseT(m) {}

    long numRows() const { return BaseT::size1; }
    long numColumns() const { return BaseT::size2; }

    void resize(long num_rows, long num_cols) { BaseT::resize(num_rows, num_cols); }

    // This is needed to avoid a Visual Studio warning about inheriting via dominance
    void makeZero() { BaseT::makeZero(); }

    /** Get the column index of each logically non-zero entry. */
    TheaArray<Index2DT> const & getColumnIndices() const { return BaseT::indices2; }

    /** Get the column index of each logically non-zero entry. Assert isValid() if you modify this array. */
    TheaArray<Index2DT> & getColumnIndices() { return BaseT::indices2; }

    /** Get the index (in the values array) of the starting element of each row. */
    TheaArray<Index1DT> const & getRowIndices() const { return BaseT::indices1; }

    /** Get the index (in the values array) of the starting element of each row. Assert isValid() if you modify this array. */
    TheaArray<Index1DT> & getRowIndices() { return BaseT::indices1; }

    /**
     * Multiply a vector \a v by this matrix (M), yielding the product \a result = M * \a v.
     *
     * @param v The vector to be multiplied, must contain numColumns() elements.
     * @param result The result vector, must be preallocated to numRows() elements.
     */
    template <typename U> void postmulVector(U const * v, U * result) const
    {
      array_size_t curr_pos = 0;
      U * rp = result;
      for (long row = 0; row < BaseT::size1; ++row, ++rp)
      {
        *rp = static_cast<U>(0);
        long num_elems = (long)(BaseT::indices1[row + 1] - BaseT::indices1[row]);
        for (long e = 0; e < num_elems; ++e, ++curr_pos)
        {
          Index2DT col = BaseT::indices2[curr_pos];
          *rp += static_cast<U>(BaseT::values[curr_pos] * v[col]);
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

}; // class CompressedRowMatrix

/** Column-major sparse matrix. Non-zero values are packed contiguously column-by-column. */
template <typename T, typename Index2DT = int, typename Index1DT = int>
class /* THEA_API */ CompressedColumnMatrix : public CompressedSparseMatrix<T, MatrixLayout::COLUMN_MAJOR, Index2DT,
                                                                            Index1DT>,
                                              public ResizableMatrix<T>
{
  private:
    typedef CompressedSparseMatrix<T, MatrixLayout::COLUMN_MAJOR, Index2DT, Index1DT> BaseT;

  public:
    THEA_DEF_POINTER_TYPES(CompressedColumnMatrix, shared_ptr, weak_ptr)

    /**
     * Default constructor. Initializes a zero matrix of the specified size. If the number of columns is omitted or zero, a
     * square matrix is created.
     */
    CompressedColumnMatrix(long num_rows = 0, long num_cols = 0) : BaseT(num_cols, num_rows) {}

    /** Copy contructor. */
    CompressedColumnMatrix(CompressedColumnMatrix const & src) : BaseT(src) {}

    /** Templated copy constructor. */
    template <typename MatrixT> explicit CompressedColumnMatrix(MatrixT const & m) : BaseT(m) {}

    long numRows() const { return BaseT::size2; }
    long numColumns() const { return BaseT::size1; }

    void resize(long num_rows, long num_cols) { BaseT::resize(num_cols, num_rows); }

    // This is needed to avoid a Visual Studio warning about inheriting via dominance
    void makeZero() { BaseT::makeZero(); }

    /** Get the row index of each logically non-zero entry. */
    TheaArray<Index2DT> const & getRowIndices() const { return BaseT::indices2; }

    /** Get the row index of each logically non-zero entry. Assert isValid() if you modify this array. */
    TheaArray<Index2DT> & getRowIndices() { return BaseT::indices2; }

    /** Get the index (in the values array) of the starting element of each column. */
    TheaArray<Index1DT> const & getColumnIndices() const { return BaseT::indices1; }

    /**
     * Get the index (in the values array) of the starting element of each column. Assert isValid() if you modify this array.
     */
    TheaArray<Index1DT> & getColumnIndices() { return BaseT::indices1; }

    /**
     * Multiply a vector \a v by this matrix (M), yielding the product \a result = M * \a v.
     *
     * @param v The vector to be multiplied, must contain numColumns() elements.
     * @param result The result vector, must be preallocated to numRows() elements.
     */
    template <typename U> void postmulVector(U const * v, U * result) const
    {
      // Initialize output to zero
      U * r_end = result + BaseT::size2;
      for (U * rp = result; rp != r_end; ++rp) *rp = static_cast<U>(0);

      array_size_t curr_pos = 0;
      for (long col = 0; col < BaseT::size1; ++col)
      {
        long num_elems = (long)(BaseT::indices1[col + 1] - BaseT::indices1[col]);
        for (long e = 0; e < num_elems; ++e, ++curr_pos)
        {
          Index2DT row = BaseT::indices2[curr_pos];
          result[row] += static_cast<U>(BaseT::values[curr_pos] * v[col]);
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

}; // class CompressedColumnMatrix

} // namespace Thea

#endif
