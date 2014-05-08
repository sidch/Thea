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

#ifndef __Thea_AddressableMatrix_hpp__
#define __Thea_AddressableMatrix_hpp__

#include "IteratableMatrix.hpp"

namespace Thea {

/**
 * Interface for a 2D matrix with elements that may be directly accessed via a (row, column) pair. The type T should be a
 * well-ordered field, pipeable to a <code>std::ostream</code>.
 */
template <typename T>
class /* THEA_API */ AddressableMatrix : public virtual IteratableMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(AddressableMatrix, shared_ptr, weak_ptr)

    typedef std::pair<long, long> IndexPair;  ///< A (row, column) index pair.
    typedef std::pair<IndexPair, T> Entry;    ///< An entry in the matrix, mapping a (row, column) pair to a value.

    /** Generic read-only iterator for an addressable matrix. */
    class ConstIterator
    {
      public:
        /** Constructor. */
        ConstIterator(AddressableMatrix const & m_, long r = 0, long c = 0)
        : m(m_), nrows(m_.numRows()), ncols(m_.numColumns()), entry(IndexPair(r, c), 0)
        {}

        /** Get the row of the current element. */
        long getRow() const { return entry.first.first; }

        /** Get the column of the current element. */
        long getColumn() const { return entry.first.second; }

        /** Dereference to the current element. */
        Entry const & operator*() const
        {
          entry.second = m.get(entry.first.first, entry.first.second);
          return entry;
        }

        /** Arrow operator for dereferencing. */
        Entry const * operator->() const
        {
          return &(this->operator*());
        }

        /** Pre-increment. */
        ConstIterator & operator++()
        {
          if (++entry.first.second >= ncols)
          {
            ++entry.first.first;
            entry.first.second = 0;
          }

          return *this;
        }

        /** Post-increment. */
        ConstIterator operator++(int)
        {
          ConstIterator old = *this;
          this->operator++();
          return old;
        }

        /** Test for equality. */
        bool operator==(ConstIterator const & rhs) const
        {
          return !(*this != rhs);
        }

        /** Test for inequality. */
        bool operator!=(ConstIterator const & rhs) const
        {
          debugAssertM(&m == &rhs.m, "Matrix: Comparing iterators from different matrices for equality");
          return (entry.first != rhs.entry.first);
        }

      private:
        AddressableMatrix const & m;
        long nrows, ncols;
        mutable Entry entry;

    }; // class ConstIterator

    /** Get an iterator pointing to the beginning of the matrix. */
    ConstIterator begin() const { return ConstIterator(*this, 0, 0); }

    /** Get an iterator pointing to the end of the matrix. */
    ConstIterator end() const { return ConstIterator(*this, (this->numColumns() > 0 ? this->numRows() : 0), 0); }

    /**
     * Get element. Most derived classes define operator() to access an element quicker, without the virtual function overhead.
     * Use this function only in generic algorithms that need polymorphic access to matrices without using templates.
     */
    virtual T const & get(long row, long col) const = 0;

    /**
     * Get an element that can be directly modified. Most derived classes define operator() to access an element quicker,
     * without the virtual function overhead. Use this function only in generic algorithms that need polymorphic access to
     * matrices without using templates.
     */
    virtual T & getMutable(long row, long col) = 0;

    /**
     * Set element. Most derived classes define operator() to access an element quicker, without the virtual function overhead.
     * Use this function only in generic algorithms that need polymorphic access to matrices without using templates.
     */
    virtual void set(long row, long col, T const & value) = 0;

    void makeZero() { fill(static_cast<T>(0)); }

    /** Set all elements of the matrix to a given value. */
    virtual void fill(T const & value)
    {
      long nr = this->numRows(), nc = this->numColumns();

      for (long r = 0; r < nr; ++r)
        for (long c = 0; c < nc; ++c)
          set(r, c, value);
    }

    /** Get the minimum element (according to signed comparison) in the matrix. */
    virtual T const & min() const
    {
      long nr = this->numRows(), nc = this->numColumns();
      T const * m = NULL;

      for (long r = 0; r < nr; ++r)
        for (long c = 0; c < nc; ++c)
        {
          T const & e = get(r, c);
          if (!m || e < *m)
            m = &e;
        }

      debugAssertM(m, "AddressableMatrix: No minimum element found");
      return *m;
    }

    /** Get the maximum element (according to signed comparison) in the matrix. */
    virtual T const & max() const
    {
      long nr = this->numRows(), nc = this->numColumns();
      T const * m = NULL;

      for (long r = 0; r < nr; ++r)
        for (long c = 0; c < nc; ++c)
        {
          T const & e = get(r, c);
          if (!m || e > *m)
            m = &e;
        }

      debugAssertM(m, "AddressableMatrix: No maximum element found");
      return *m;
    }

    /** Get a row of the matrix. \a values must be preallocated with numColumns() elements. */
    virtual void getRow(long row, T * values) const
    {
      long nc = this->numColumns();
      for (long c = 0; c < nc; ++c)
        values[c] = get(row, c);
    }

    /** Set a row of the matrix. \a values must contain numColumns() elements. */
    virtual void setRow(long row, T const * values)
    {
      long nc = this->numColumns();
      for (long c = 0; c < nc; ++c)
        set(row, c, values[c]);
    }

    /** Get a column of the matrix. \a values must be preallocated with numRows() elements. */
    virtual void getColumn(long col, T * values) const
    {
      long nr = this->numRows();
      for (long r = 0; r < nr; ++r)
        values[r] = get(r, col);
    }

    /** Set a column of the matrix. \a values must contain numRows() elements. */
    virtual void setColumn(long col, T const * values)
    {
      long nr = this->numRows();
      for (long r = 0; r < nr; ++r)
        set(r, col, values[r]);
    }

    /** Get a string representing the matrix. */
    virtual std::string toString() const
    {
      std::ostringstream oss;
      oss << '[';

      long nr = this->numRows(), nc = this->numColumns();
      for (long r = 0; r < nr; ++r)
      {
        for (long c = 0; c < nc; ++c)
        {
          oss << get(r, c);
          if (c < nc - 1) oss << ", ";
        }

        if (r < nr - 1) oss << "; ";
      }

      oss << ']';
      return oss.str();
    }

    /**
     * Copy the contents of the matrix to another addressable matrix. The matrices must have the same dimensions, else an
     * assertion failure is generated.
     */
    template <typename U> void copyTo(AddressableMatrix<U> & dst) const
    {
      long nr = this->numRows(), nc = this->numColumns();
      alwaysAssertM(nr == dst.numRows() && nc == dst.numColumns(),
                    "AddressableMatrix: Copy destination has different dimensions");

      for (long r = 0; r < nr; ++r)
        for (long c = 0; c < nc; ++c)
          dst.set(r, c, get(r, c));
    }

}; // class AddressableMatrix

/** Pipe a textual representation of an addressable matrix to a <code>std::ostream</code>. */
template <typename T>
std::ostream &
operator<<(std::ostream & os, AddressableMatrix<T> const & m)
{
  return os << m.toString();
}

} // namespace Thea

#endif
