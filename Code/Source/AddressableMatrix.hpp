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

#include "BasicMatrix.hpp"

namespace Thea {

/**
 * Interface for a 2D matrix with elements that may be directly accessed via a (row, column) pair. The type T should be a
 * well-ordered field, pipeable to a <code>std::ostream</code>.
 */
template <typename T>
class /* THEA_API */ AddressableMatrix : public virtual BasicMatrix<T>
{
  private:
    typedef BasicMatrix<T> BaseType;

  public:
    THEA_DEF_POINTER_TYPES(AddressableMatrix, shared_ptr, weak_ptr)

    /** Generic iterator for an addressable matrix. */
    class Iterator
    {
      public:
        /** Constructor. */
        Iterator(AddressableMatrix const & m_, long r = 0, long c = 0)
        : m(m_), row(r), col(c), nrows(m_.numRows()), ncols(m_.numColumns())
        {}

        /** Get the row of the current element. */
        long getRow() const { return row; }

        /** Get the column of the current element. */
        long getColumn() const { return col; }

        /** Get the value of the current element. */
        T const & operator*() const { return m.get(row, col); }

        /** Pre-increment. */
        Iterator & operator++()
        {
          if (++col >= ncols)
          {
            ++row;
            col = 0;
          }

          return *this;
        }

        /** Post-increment. */
        Iterator operator++(int)
        {
          Iterator old = *this;
          this->operator++();
          return old;
        }

      private:
        AddressableMatrix const & m;
        long row, col;
        long nrows, ncols;

    }; // class Iterator

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

    /** Get the minimum element in the matrix. */
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

    /** Get the maximum element in the matrix. */
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
