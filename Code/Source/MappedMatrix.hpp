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
// First version: 2019
//
//============================================================================

#ifndef __Thea_MappedMatrix_hpp__
#define __Thea_MappedMatrix_hpp__

#include "IAddressableMatrix.hpp"
#include "ISparseMatrix.hpp"
#include "Hash.hpp"
#include "UnorderedSet.hpp"
#include <utility>

namespace Thea {

/** Matrix that stores elements as a (row, col) --> value map. */
template <typename T = Real, typename IndexT = intx>
class /* THEA_API */ MappedMatrix : public virtual IAddressableMatrix<T>, public virtual ISparseMatrix<T>
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

    int64 THEA_ICALL rows() const { return num_rows; }
    int64 THEA_ICALL cols() const { return num_cols; }
    void THEA_ICALL setZero() { triplets.clear(); }
    int8 THEA_ICALL isResizable() const { return true; }

    int8 THEA_ICALL resize(int64 nrows, int64 ncols)
    {
      alwaysAssertM(nrows >= 0 && ncols >= 0,
                    format("MappedMatrix: Cannot resize to negative dimensions %ldx%ld", (long)nrows, (long)ncols));

      triplets.clear();
      num_rows = nrows;
      num_cols = ncols;

      return true;
    }

    int64 THEA_ICALL numStoredElements() const { return (int64)triplets.size(); }

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
     * \copydoc IAddressableMatrix::at()
     *
     * If the element is not already stored in the matrix, it will return a zero value.
     */
    T const & THEA_ICALL at(int64 row, int64 col) const
    {
      return (*this)((Index)row, (Index)col);
    }

    /**
     * \copydoc IAddressableMatrix::mutableAt()
     *
     * If the element is not already stored in the matrix, <b>it will be inserted</b> with an initial value of zero, increasing
     * the storage size.
     */
    T & THEA_ICALL mutableAt(int64 row, int64 col)
    {
      auto loc = triplets.insert(Triplet((Index)row, (Index)col, static_cast<T>(0)));  // insert a zero if not already set
      return const_cast<T &>(loc.first->value());  // either the existing element, or the newly inserted one.
    }

    /**
     * \copydoc IAddressableMatrix::getRow()
     *
     * This function is in general <b>very slow</b> because it has to iterate through all the stored elements, and should be
     * avoided.
     */
    void THEA_ICALL getRow(int64 row, T * values) const
    {
      std::fill(values, values + num_cols, static_cast<T>(0));

      for (auto const & t : triplets)
        if (t.row() == row)
          values[t.col()] = t.value();
    }

    /**
     * \copydoc IAddressableMatrix::setRow()
     *
     * This function compares the elements of \a values to 0 using the not-equals operator <tt>T::operator!=</tt>. This may be
     * prone to numerical (floating-point) error.
     */
    void THEA_ICALL setRow(int64 row, T const * values)
    {
      for (intx c = 0; c < num_cols; ++c)
        if (values[c] != static_cast<T>(0))
          set((Index)row, (Index)c, values[c]);
    }

    /**
     * \copydoc IAddressableMatrix::getColumn()
     *
     * This function is in general <b>very slow</b> because it has to iterate through all the stored elements, and should be
     * avoided.
     */
    void THEA_ICALL getColumn(int64 col, T * values) const
    {
      std::fill(values, values + num_rows, static_cast<T>(0));

      for (auto const & t : triplets)
        if (t.col() == col)
          values[t.row()] = t.value();
    }

    /**
     * \copydoc IAddressableMatrix::setColumn()
     *
     * This function compares the elements of \a values to 0 using the not-equals operator <tt>T::operator!=</tt>. This may be
     * prone to numerical (floating-point) error.
     */
    void THEA_ICALL setColumn(int64 col, T const * values)
    {
      for (intx r = 0; r < num_rows; ++r)
        if (values[r] != static_cast<T>(0))
          set((Index)r, (Index)col, values[r]);
    }

    // Type-casting functions
    IAddressableMatrix<T> const * THEA_ICALL asAddressable() const { return this; }
    IAddressableMatrix<T> * THEA_ICALL asAddressable() { return this; }
    ISparseMatrix<T> const * THEA_ICALL asSparse() const { return this; }
    ISparseMatrix<T> * THEA_ICALL asSparse() { return this; }
    IDenseMatrix<T> const * THEA_ICALL asDense() const { return nullptr; }
    IDenseMatrix<T> * THEA_ICALL asDense() { return nullptr; }
    ICompressedSparseMatrix<T> const * THEA_ICALL asCompressed() const { return nullptr; }
    ICompressedSparseMatrix<T> * THEA_ICALL asCompressed() { return nullptr; }

  private:
    intx num_rows, num_cols;
    TripletSet triplets;  ///< The set of stored triplets.

}; // class MappedMatrix

} // namespace Thea

#endif
