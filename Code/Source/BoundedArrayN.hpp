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
// First version: 2011
//
//============================================================================

#ifndef __Thea_BoundedArrayN_hpp__
#define __Thea_BoundedArrayN_hpp__

#include "Common.hpp"
#include <algorithm>
#include <iterator>

namespace Thea {

/**
 * An array of bounded maximum size N, stored on the stack. No new elements can be added to the end of the array once it reaches
 * its maximum capacity, until some are removed via erase(). Elements can be inserted in the middle of the array even when it is
 * full: in this case the last element is dropped to make space. This class is useful for very fast allocation of space for a
 * few elements, where the exact number of elements is not known but is guaranteed to have an upper limit.
 *
 * The implementation always allocates enough space to store the maximum number of instances of T. The capacity N should be
 * <b>positive</b> (non-zero).
 *
 * For a bounded stack-based array that also guarantees that its elements are in sorted order, see BoundedSortedArrayN.
 */
template <size_t N, typename T>
class BoundedArrayN
{
  private:
    static_assert(N > 0, "BoundedArrayN: Capacity must be positive");

    size_t num_elems;
    T values[N];

  public:
    typedef T *                                    iterator;                ///< Iterator over elements.
    typedef T const *                              const_iterator;          ///< Iterator over immutable elements.
    typedef std::reverse_iterator<iterator>        reverse_iterator;        ///< Reverse iterator over elements.
    typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;  ///< Reverse iterator over immutable elements.

    /** Constructor. */
    BoundedArrayN() : num_elems(0) {}

    /** Get the maximum number of elements the array can hold. */
    static constexpr size_t getCapacity() { return N; }

    /** Get the number of elements in the array. */
    size_t size() const { return num_elems; }

    /** Check if the array is empty or not. */
    bool empty() const { return num_elems <= 0; }

    /** Check if the array has reached maximum capacity or not. */
    bool isFull() const { return num_elems >= N; }

    /** Get a pointer to the data buffer. */
    T const * data() const { return values; }

    /** Get a pointer to the data buffer. */
    T * data() { return values; }

    iterator begin() noexcept               { return values; }  ///< Iterator to the beginning of the array.
    const_iterator begin() const noexcept   { return values; }  ///< Iterator to the beginning of the array.
    const_iterator cbegin() const noexcept  { return values; }  ///< Iterator to the beginning of the array.
    iterator end() noexcept                 { return values + num_elems; }  ///< Iterator to the end of the array.
    const_iterator end() const noexcept     { return values + num_elems; }  ///< Iterator to the end of the array.
    const_iterator cend() const noexcept    { return values + num_elems; }  ///< Iterator to the end of the array.

    /**< Reverse iterator to the end of the array. */
    reverse_iterator rbegin() noexcept               { return reverse_iterator(end()); }

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator rbegin() const noexcept   { return reverse_iterator(end()); }

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator crbegin() const noexcept  { return reverse_iterator(end()); }

    /**< Reverse iterator to the beginning of the array. */
    reverse_iterator rend() noexcept                 { return reverse_iterator(begin()); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator rend() const noexcept     { return reverse_iterator(begin()); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator crend() const noexcept    { return reverse_iterator(begin()); }

    /** Get the first element in the array. */
    T const & front() const
    {
      theaAssertM(num_elems > 0, "BoundedArrayN: Can't get first element of empty array");
      return values[0];
    }

    /** Get the last element in the array. */
    T const & back() const
    {
      theaAssertM(num_elems > 0, "BoundedArrayN: Can't get last element of empty array");
      return values[num_elems - 1];
    }

    /** Get the element at a given position in the array. Bounds checks are only performed in debug mode. */
    T const & operator[](size_t i) const
    {
      theaAssertM(i < num_elems, format("BoundedArrayN: Index %ld out of bounds [0, %ld)", (long)i, (long)num_elems));
      return values[i];
    }

    /** Get the element at a given position in the array. Bounds checks are only performed in debug mode. */
    T & operator[](size_t i)
    {
      theaAssertM(i < num_elems, format("BoundedArrayN: Index %ld out of bounds [0, %ld)", (long)i, (long)num_elems));
      return values[i];
    }

    /** Add a new element to the end of array, if it is not full. If the array is full, the operation fails silently. */
    void append(T const & t)
    {
      if (num_elems < N)
        values[num_elems++] = t;
#ifdef THEA_DEBUG_BUILD
      else
        THEA_WARNING << "BoundedArrayN: Cannot append element to a full array of size " << N;
#endif
    }

    /**
     * Add a new element to the end of array, if it is not full (STL-style syntax). If the array is full, the operation fails
     * silently.
     */
    void push_back(T const & t) { append(t); }

    /**
     * Insert an element at a given position in the array, shifting all existing elements at or after this position up by one
     * position to make space. If the array is full, the last element is dropped.
     */
    void insert(size_t i, T const & t)
    {
      theaAssertM(i <= std::min(num_elems, N - 1),
                  format("BoundedArrayN: Index %ld out of bounds [0, %ld]", (long)i, (long)std::min(num_elems, N - 1)));

      if (isFull())
        std::copy_backward(values + i, values + N - 1, values + N);
      else if (i < num_elems)
      {
        std::copy_backward(values + i, values + num_elems, values + num_elems);
        ++num_elems;
      }

      values[i] = t;
    }

    /** Remove the element at the given position from the array. */
    void erase(size_t i)
    {
      if (i < num_elems)
      {
        std::copy(values + i + 1, values + num_elems, values + i);
        --num_elems;
      }
    }

    /** Remove all elements from the array. */
    void clear()
    {
      num_elems = 0;
    }

}; // class BoundedArrayN

} // namespace Thea

#endif
