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

#ifndef __Thea_BoundedSortedArray_hpp__
#define __Thea_BoundedSortedArray_hpp__

#include "Common.hpp"
#include <algorithm>
#include <functional>
#include <iterator>

namespace Thea {

/**
 * A sorted array of bounded maximum size, ordered in ascending order according to a comparator. If the array is full and a new
 * element is added, the last element is dropped.
 *
 * The implementation always allocates enough space to store the maximum number of instances of T. This space is allocated on
 * the heap: if the capacity is known at compile-time, BoundedSortedArrayN, which can use stack storage, is usually preferred.
 *
 * @see BoundedSortedArrayN
 */
template < typename T, typename Compare = std::less<T> >
class BoundedSortedArray
{
  private:
    Compare compare;
    size_t capacity, num_elems;
    T * values;

  public:
    typedef T const *                              const_iterator;          ///< Iterator over immutable elements.
    typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;  ///< Reverse iterator over immutable elements.

    /**
     * Constructor. Allocates memory for \a capacity_ elements.
     *
     * @param capacity_ The maximum number of elements the array can hold. Must be a non-negative integer.
     * @param compare_ The comparator that evaluates the "less-than" operator on objects of type T.
     */
    BoundedSortedArray(size_t capacity_ = 0, Compare compare_ = Compare())
    : compare(compare_), capacity(capacity_), num_elems(0), values(nullptr)
    {
      alwaysAssertM(capacity >= 0, "BoundedSortedArray: Capacity must be non-negative");

      if (capacity > 0)
        values = new T[capacity];
    }

    /** Copy constructor. */
    BoundedSortedArray(BoundedSortedArray const & src)
    : compare(src.compare), capacity(src.capacity), num_elems(src.num_elems),
      values(src.capacity > 0 ? new T[src.capacity] : nullptr)
    {
      if (src.num_elems > 0)
        std::copy(src.values, src.values + src.num_elems, values);
    }

    /** Assignment operator. */
    BoundedSortedArray & operator=(BoundedSortedArray const & src)
    {
      compare = src.compare;
      num_elems = src.num_elems;

      if (capacity != src.capacity)
      {
        delete [] values;
        values = (src.capacity > 0 ? new T[src.capacity] : nullptr);
        capacity = src.capacity;
      }

      if (src.num_elems > 0)
        std::copy(src.values, src.values + src.num_elems, values);

      return *this;
    }

    /** Destructor. */
    ~BoundedSortedArray() { delete [] values; }

    /** Set the maximum number of elements the array can hold. The array is cleared and all prior data is discarded. */
    void setCapacity(size_t capacity_)
    {
      if (capacity != capacity_)
      {
        delete [] values;
        values = (capacity_ > 0 ? new T[capacity_] : nullptr);
        capacity = capacity_;
      }

      num_elems = 0;
    }

    /** Get the maximum number of elements the array can hold. */
    size_t getCapacity() const { return capacity; }

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
    T const & front() const
    {
      debugAssertM(num_elems > 0, "BoundedSortedArray: Can't get first element of empty array");
      return values[0];
    }

    /** Get the last element in the sorted sequence. */
    T const & back() const
    {
      debugAssertM(num_elems > 0, "BoundedSortedArray: Can't get last element of empty array");
      return values[num_elems - 1];
    }

    /** Get the element at a given position in the sorted sequence. */
    T const & operator[](size_t i) const
    {
      debugAssertM(i < num_elems, format("BoundedSortedArray: Index %ld out of bounds [0, %ld)", (long)i, (long)num_elems));
      return values[i];
    }

    /** Check if the array contains an element with a given value. */
    bool contains(T const & t) const { return find(t) < capacity; }

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
      return (lb < num_elems && !compare(t, values[lb])) ? lb : capacity;
    }

    /**
     * Get the index of the first element strictly greater than \a t, or return the capacity of the array if no such element is
     * present.
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
     * Get the index of the first element equal to or greater than \a t, or return the capacity of the array if no such element
     * is present.
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
      return capacity > 0 && (num_elems < capacity || compare(t, back()));
    }

    /**
     * Insert a value into the array.
     *
     * @return The index of the newly inserted element, or the capacity of the array if the value could not be inserted.
     */
    size_t insert(T const & t)
    {
      if (capacity <= 0)
        return capacity;

      if (num_elems <= 0)
      {
        values[0] = t;
        ++num_elems;
        return 0;
      }
      else if (isInsertable(t))
      {
        size_t ub = upperBound(t);
        T * end = values + (num_elems < capacity ? num_elems : capacity - 1);
        std::copy_backward(values + ub, end, end + 1);
        values[ub] = t;
        if (num_elems < capacity) ++num_elems;
        return ub;
      }

      return capacity;
    }

    /**
     * Insert a value into the array only if it does not already exist.
     *
     * @return The index of the newly inserted element, or the capacity of the array if the value could not be inserted.
     *
     * @todo Make this faster by merging the containment test with the lookup for the insertion position.
     */
    size_t insertUnique(T const & t)
    {
      if (contains(t))
        return capacity;

      return insert(t);
    }

    /** Remove the element at the given position from the array. */
    void erase(size_t i)
    {
      if (i >= 0 && i < num_elems)
      {
        std::copy(values + i + 1, values + num_elems, values + i);
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

}; // class BoundedSortedArray

} // namespace Thea

#endif
