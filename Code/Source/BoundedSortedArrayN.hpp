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
// First version: 2010
//
//============================================================================

#ifndef __Thea_BoundedSortedArrayN_hpp__
#define __Thea_BoundedSortedArrayN_hpp__

#include "Common.hpp"
#include "Algorithms/FastCopy.hpp"
#include <functional>
#include <iterator>

namespace Thea {

/**
 * A sorted array of bounded maximum size N, ordered in ascending order according to a comparator. If the array is full and a
 * new element is added, the last element is dropped. The capacity is set as a template parameter N, so the array can be stored
 * entirely on the stack. If the array size is known at compile-time, this class is usually a more efficient alternative to
 * BoundedSortedArray.
 *
 * To get some extra speed when T has a trivial (bit-copy) assignment operator, make sure that
 * <tt>std::is_trivially_copyable</tt> is true for T.
 *
 * The implementation always allocates enough space to store the maximum number of instances of T. The capacity N should be
 * <b>positive</b> (non-zero).
 *
 * @see BoundedSortedArray
 */
template < size_t N, typename T, typename Compare = std::less<T> >
class BoundedSortedArrayN
{
  private:
    static_assert(N > 0, "BoundedSortedArrayN: Capacity must be positive");

    Compare compare;
    size_t num_elems;
    T values[N];

  public:
    typedef T const *                              const_iterator;          ///< Iterator over immutable elements.
    typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;  ///< Reverse iterator over immutable elements.

    /**
     * Constructor.
     *
     * @param compare_ The comparator that evaluates the "less-than" operator on objects of type T.
     */
    BoundedSortedArrayN(Compare compare_ = Compare()) : compare(compare_), num_elems(0) {}

    /** Get the maximum number of elements the array can hold. */
    static constexpr size_t getCapacity() { return N; }

    /** Get the number of elements in the array. */
    size_t size() const { return num_elems; }

    /** Check if the array is empty or not. */
    bool empty() const { return num_elems <= 0; }

    /** Get a pointer to the data buffer. */
    T const * data() const { return values; }

    const_iterator begin() const noexcept   { return values; }  ///< Iterator to the beginning of the array.
    const_iterator cbegin() const noexcept  { return values; }  ///< Iterator to the beginning of the array.
    const_iterator end() const noexcept     { return values + num_elems; }  ///< Iterator to the end of the array.
    const_iterator cend() const noexcept    { return values + num_elems; }  ///< Iterator to the end of the array.

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator rbegin() const noexcept   { return reverse_iterator(end()); }

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator crbegin() const noexcept  { return reverse_iterator(end()); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator rend() const noexcept     { return reverse_iterator(begin()); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator crend() const noexcept    { return reverse_iterator(begin()); }

    /** Get the first element in the sorted sequence. */
    T const & first() const
    {
      debugAssertM(num_elems > 0, "BoundedSortedArrayN: Can't get first element of empty array");
      return values[0];
    }

    /** Get the last element in the sorted sequence. */
    T const & last() const
    {
      debugAssertM(num_elems > 0, "BoundedSortedArrayN: Can't get last element of empty array");
      return values[num_elems - 1];
    }

    /** Get the element at a given position in the sorted sequence. */
    T const & operator[](size_t i) const
    {
      debugAssertM(i < num_elems, format("BoundedSortedArrayN: Index %ld out of bounds [0, %ld)", (long)i, (long)num_elems));
      return values[i];
    }

    /** Check if the array contains an element with a given value. */
    bool contains(T const & t) const { return find(t) < N; }

    /**
     * Check if the array already contains an element with a given value, by testing every element in the set for equality with
     * the query. This is useful when searching with other notions of equality than that defined by the ordering comparator.
     */
    template <typename EqualityComparatorT> bool contains(T const & t, EqualityComparatorT const & comp) const
    {
      for (size_t i = 0; i < num_elems; ++i)
        if (comp(values[i], t))
          return true;

      return false;
    }

    /**
     * Get the index of a given value in the array. If the value is not present in the array, the capacity of the array is
     * returned. If the value occurs multiple times, the index of any one occurrence is returned.
     */
    size_t find(T const & t) const
    {
      size_t lb = lowerBound(t);
      return (lb < num_elems && !compare(t, values[lb])) ? lb : N;
    }

    /**
     * Get the index of the first element strictly greater than \a t, or return the capacity of the array (N) if no such element
     * is present.
     */
    size_t upperBound(T const & t) const
    {
      size_t first = 0, mid, step;
      size_t count = num_elems;
      while (count > 0)
      {
        step = count >> 1;
        mid = first + step;
        if (!compare(t, values[mid]))
        {
          first = mid + 1;
          count -= (step + 1);
        }
        else
          count = step;
      }

      return first;
    }

    /**
     * Get the index of the first element equal to or greater than \a t, or return the capacity of the array (N) if no such
     * element is present.
     */
    size_t lowerBound(T const & t) const
    {
      size_t first = 0, mid, step;
      size_t count = num_elems;
      while (count > 0)
      {
        step = count >> 1;
        mid = first + step;
        if (compare(values[mid], t))
        {
          first = mid + 1;
          count -= (step + 1);
        }
        else
          count = step;
      }

      return first;
    }

    /**
     * Check if a value can be inserted in the array. This requires that either the array has fewer elements than its capacity,
     * or the value is "less than" the last element in the array.
     */
    bool isInsertable(T const & t) const
    {
      return num_elems < N || compare(t, last());
    }

    /**
     * Insert a value into the array.
     *
     * @return The index of the newly inserted element, or the capacity of the array (N) if the value could not be inserted.
     */
    size_t insert(T const & t)
    {
      if (num_elems <= 0)
      {
        values[0] = t;
        ++num_elems;
        return 0;
      }
      else if (isInsertable(t))
      {
        size_t ub = upperBound(t);
        T * end = values + (num_elems < N ? num_elems : N - 1);
        Algorithms::fastCopyBackward(values + ub, end, end + 1);
        values[ub] = t;
        if (num_elems < N) ++num_elems;
        return ub;
      }

      return N;
    }

    /**
     * Insert a value into the array only if it does not already exist.
     *
     * @return The index of the newly inserted element, or the capacity of the array (N) if the value could not be inserted.
     *
     * @todo Make this faster by merging the containment test with the lookup for the insertion position.
     */
    size_t insertUnique(T const & t)
    {
      if (contains(t))
        return N;

      return insert(t);
    }

    /** Remove the element at the given position from the array. */
    void erase(size_t i)
    {
      if (i < num_elems)
      {
        Algorithms::fastCopy(values + i + 1, values + num_elems, values + i);
        --num_elems;
      }
    }

    /** Remove (one occurrence of) the given value from the array, if it is present. */
    void erase(T const & t)
    {
      erase(find(t));
    }

    /** Remove all elements from the array. */
    void clear()
    {
      num_elems = 0;
    }

}; // class BoundedSortedArrayN

} // namespace Thea

#endif
