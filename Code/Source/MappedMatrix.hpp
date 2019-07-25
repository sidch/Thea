//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2019, Siddhartha Chaudhuri
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

#ifndef __Thea_MappedMatrix_hpp__
#define __Thea_MappedMatrix_hpp__

#include "AbstractAddressableMatrix.hpp"
#include "AbstractSparseMatrix.hpp"
#include "Hash.hpp"
#include "UnorderedSet.hpp"
#include <utility>

namespace Thea {

/** Matrix that stores elements as a (row, col) --> value map. */
template <typename T = Real, typename IndexT = intx>
class /* THEA_API */ MappedMatrix : public virtual AbstractAddressableMatrix<T>, public virtual AbstractSparseMatrix<T>
{
  public:
    typedef IndexT Index;  ///< The type of row or column indices.

    /**
     * The stored (%row, %col, %value) triplet type. We use a set of triplets, rather than a (%row, %col) --> %value map, for
     * easy interoperability with the SparseMatrix::setFromTriplets() function. Implementionally, the two are similar since
     * std::unordered_map is a set of (key, %value) pairs.
     */
    class Triplet
    {
      public:
        /** Default constructor. */
        Triplet() {}

        /** Construct from %row, %column and %value fields. */
        Triplet(IndexT r_, IndexT c_, T const & v_ = static_cast<T>(0)) : r(r_), c(c_), v(v_) {}

        /** Get the %row field. */
        Index row() const { return r; }

        /** Get the %column field. */
        Index col() const { return c; }

        /** Get the %value field. */
        T const & value() const { return v; }

        /** Get the %value field. */
        T & value() { return v; }

        /** Get a string representing the triplet. */
        std::string toString() const
        {
          std::ostringstream oss;
          oss << '(' << r << ", " << c << ", " << v << ')';
          return oss.str();
        }

      private:
        IndexT r, c;
        T v;

    }; // class Triplet

  private:

    /** A hash function for the (row, col) key of a triplet. */
    struct TripletHash
    {
      size_t operator()(Triplet const & t) const
      {
        return std::hash< std::pair<Index, Index> >()(std::make_pair(t.row(), t.col()));
      }
    };

    /** Compare the (row, col) keys of two triplets for equality. */
    struct TripletEqualTo
    {
      size_t operator()(Triplet const & a, Triplet const & b) const
      {
        return a.row() == b.row() && a.col() == b.col();
      }
    };

    typedef UnorderedSet<Triplet, TripletHash, TripletEqualTo> TripletSet;  ///< Storage container for triplets.

  public:
    THEA_DECL_SMART_POINTERS(MappedMatrix)

    typedef typename TripletSet::const_iterator TripletConstIterator;  ///< Read-only iterator over triplets.

    /** Constructor. */
    MappedMatrix(intx nrows = 0, intx ncols = 0)
    : num_rows(nrows), num_cols(ncols)
    {
      alwaysAssertM(nrows >= 0 && ncols >= 0,
                    format("MappedMatrix: Cannot resize to negative dimensions %ldx%ld", (long)nrows, (long)ncols));
    }

    int64 rows() const { return num_rows; }
    int64 cols() const { return num_cols; }
    void setZero() { triplets.clear(); }
    int8 isResizable() const { return true; }

    int8 resize(int64 nrows, int64 ncols)
    {
      alwaysAssertM(nrows >= 0 && ncols >= 0,
                    format("MappedMatrix: Cannot resize to negative dimensions %ldx%ld", (long)nrows, (long)ncols));

      triplets.clear();
      num_rows = nrows;
      num_cols = ncols;

      return true;
    }

    int64 numStoredElements() const { return (int64)triplets.size(); }

    /** Get a read-only iterator pointing to the first (row, col, value) triplet in the matrix. */
    TripletConstIterator tripletsBegin() const { return triplets.begin(); }

    /** Get a read-only iterator pointing to one past the last (row, col, value) triplet in the matrix. */
    TripletConstIterator tripletsEnd() const { return triplets.end(); }

    /**
     * Get the element at position (\a row, \a col). If the element is not already stored in the matrix, it will return a zero
     * value.
     */
    T const & operator()(Index row, Index col) const
    {
      static T const ZERO_VAL(0);

      auto loc = triplets.find(Triplet(row, col, ZERO_VAL));  // value is ignored in this lookup
      return (loc != triplets.end() ? loc->value() : ZERO_VAL);
    }

    /** Set the element at position (\a row, \a col). */
    void set(Index row, Index col, T const & value)
    {
      Triplet t(row, col, value);
      auto loc = triplets.insert(t);
      if (!loc.second)
        const_cast<T &>(loc.first->value()) = value;
    }

    /**
     * \copydoc AbstractAddressableMatrix::at()
     *
     * If the element is not already stored in the matrix, it will return a zero value.
     */
    T const & at(int64 row, int64 col) const
    {
      return (*this)((Index)row, (Index)col);
    }

    /**
     * \copydoc AbstractAddressableMatrix::mutableAt()
     *
     * If the element is not already stored in the matrix, <b>it will be inserted</b> with an initial value of zero, increasing
     * the storage size.
     */
    T & mutableAt(int64 row, int64 col)
    {
      auto loc = triplets.insert(Triplet((Index)row, (Index)col, static_cast<T>(0)));  // insert a zero if not already set
      return const_cast<T &>(loc.first->value());  // either the existing element, or the newly inserted one.
    }

    /**
     * \copydoc AbstractAddressableMatrix::getRow()
     *
     * This function is in general <b>very slow</b> because it has to iterate through all the stored elements, and should be
     * avoided.
     */
    void getRow(int64 row, T * values) const
    {
      std::fill(values, values + num_cols, static_cast<T>(0));

      for (auto const & t : triplets)
        if (t.row() == row)
          values[t.col()] = t.value();
    }

    /**
     * \copydoc AbstractAddressableMatrix::setRow()
     *
     * This function compares the elements of \a values to 0 using the not-equals operator <tt>T::operator!=</tt>. This may be
     * prone to numerical (floating-point) error.
     */
    void setRow(int64 row, T const * values)
    {
      for (intx c = 0; c < num_cols; ++c)
        if (values[c] != static_cast<T>(0))
          set((Index)row, (Index)c, values[c]);
    }

    /**
     * \copydoc AbstractAddressableMatrix::getColumn()
     *
     * This function is in general <b>very slow</b> because it has to iterate through all the stored elements, and should be
     * avoided.
     */
    void getColumn(int64 col, T * values) const
    {
      std::fill(values, values + num_rows, static_cast<T>(0));

      for (auto const & t : triplets)
        if (t.col() == col)
          values[t.row()] = t.value();
    }

    /**
     * \copydoc AbstractAddressableMatrix::setColumn()
     *
     * This function compares the elements of \a values to 0 using the not-equals operator <tt>T::operator!=</tt>. This may be
     * prone to numerical (floating-point) error.
     */
    void setColumn(int64 col, T const * values)
    {
      for (intx r = 0; r < num_rows; ++r)
        if (values[r] != static_cast<T>(0))
          set((Index)r, (Index)col, values[r]);
    }

    // Type-casting functions
    AbstractAddressableMatrix<T> const * asAddressable() const { return this; }
    AbstractAddressableMatrix<T> * asAddressable() { return this; }
    AbstractSparseMatrix<T> const * asSparse() const { return this; }
    AbstractSparseMatrix<T> * asSparse() { return this; }
    AbstractDenseMatrix<T> const * asDense() const { return NULL; }
    AbstractDenseMatrix<T> * asDense() { return NULL; }
    AbstractCompressedSparseMatrix<T> const * asCompressed() const { return NULL; }
    AbstractCompressedSparseMatrix<T> * asCompressed() { return NULL; }

  private:
    intx num_rows, num_cols;
    TripletSet triplets;  ///< The set of stored triplets.

}; // class MappedMatrix

} // namespace Thea

#endif
