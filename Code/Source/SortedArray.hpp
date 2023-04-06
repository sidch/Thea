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
// First version: 2012
//
//============================================================================

#ifndef __Thea_SortedArray_hpp__
#define __Thea_SortedArray_hpp__

#include "Common.hpp"
#include <functional>

namespace Thea {

/** An array sorted in ascending order according to a comparator. */
template < typename T, typename Compare = std::less<T> >
class SortedArray
{
  private:
    Array<T> values;
    Compare compare;

  public:
    typedef typename Array<T>::const_iterator          const_iterator;          ///< Iterator over immutable elements.
    typedef typename Array<T>::const_reverse_iterator  const_reverse_iterator;  ///< Reverse iterator over immutable elements.

    /**
     * Constructor.
     *
     * @param compare_ The comparator that evaluates the "less-than" operator on objects of type T.
     */
    SortedArray(Compare compare_ = Compare()) : compare(compare_) {}

    /** Get the number of elements in the array. */
    intx size() const { return (intx)values.size(); }

    /** Check if the array is empty or not. */
    bool empty() const { return values.empty(); }

    /** Get the first element in the sorted sequence. */
    T const & front() const
    {
      theaAssertM(!values.empty(), "SortedArray: Can't get first element of empty array");
      return values[0];
    }

    /** Get the last element in the sorted sequence. */
    T const & back() const
    {
      theaAssertM(!values.empty(), "SortedArray: Can't get last element of empty array");
      return values[size() - 1];
    }

    /** Get the element at a given position in the sorted sequence. */
    T const & operator[](intx i) const
    {
      theaAssertM(i >= 0 && i < size(), format("SortedArray: Index %d out of bounds [0, %ld)", i, size()));
      return values[(size_t)i];
    }

    const_iterator begin() const noexcept   { return values.begin();  }  ///< Iterator to the beginning of the array.
    const_iterator cbegin() const noexcept  { return values.cbegin(); }  ///< Iterator to the beginning of the array.
    const_iterator end() const noexcept     { return values.end();    }  ///< Iterator to the end of the array.
    const_iterator cend() const noexcept    { return values.cend();   }  ///< Iterator to the end of the array.

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator rbegin() const noexcept   { return values.rbegin(); }

    /**< Reverse iterator to the end of the array. */
    const_reverse_iterator crbegin() const noexcept  { return values.crbegin(); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator rend() const noexcept     { return values.rend(); }

    /**< Reverse iterator to the beginning of the array. */
    const_reverse_iterator crend() const noexcept    { return values.crend(); }

    /** Check if the array contains an element with a given value. */
    bool contains(T const & t) const { return find(t) >= 0; }

    /**
     * Check if the array already contains an element with a given value, by testing every element in the set for equality with
     * the query. This is useful when searching with other notions of equality than that defined by the ordering comparator.
     */
    template <typename EqualityComparatorT> bool contains(T const & t, EqualityComparatorT const & comp) const
    {
      for (int i = 0; i < size(); ++i)
        if (comp(values[i], t))
          return true;

      return false;
    }

    /**
     * Get the index of a given value, or negative if it is not present in the array. If the value occurs multiple times, the
     * index of any one occurrence is returned.
     */
    intx find(T const & t) const
    {
      intx lb = lowerBound(t);
      return (lb < size() && !compare(t, values[lb])) ? lb : -1;
    }

    /**
     * Get the index of the first element strictly greater than \a t, or return the size of the array if no such element is
     * present.
     */
    intx upperBound(T const & t) const
    {
      intx first = 0, mid, step;
      intx count = size();
      while (count > 0)
      {
        step = count >> 1;
        mid = first + step;
        if (!compare(t, values[(size_t)mid]))
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
     * Get the index of the first element equal to or greater than \a t, or return the size of the array if no such element is
     * present.
     */
    intx lowerBound(T const & t) const
    {
      intx first = 0, mid, step;
      intx count = size();
      while (count > 0)
      {
        step = count >> 1;
        mid = first + step;
        if (compare(values[(size_t)mid], t))
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
     * Insert a value into the array.
     *
     * @return The index of the newly inserted element, or negative if the value could not be inserted.
     */
    intx insert(T const & t)
    {
      if (values.empty())
      {
        values.push_back(t);
        return 0;
      }
      else
      {
        intx ub = upperBound(t);
        values.insert(values.begin() + (size_t)ub, t);
        return ub;
      }
    }

    /**
     * Insert a value into the array only if it does not already exist.
     *
     * @return The index of the newly inserted element, or negative if the value could not be inserted.
     *
     * @todo Make this faster by merging the containment test with the lookup for the insertion position.
     */
    intx insertUnique(T const & t)
    {
      if (contains(t))
        return -1;

      return insert(t);
    }

    /** Remove the element at the given position from the array. */
    void erase(intx i)
    {
      if (i >= 0 && i < (intx)values.size())
        values.erase(values.begin() + (size_t)i);
    }

    /** Remove (one occurrence of) the given value from the array, if it is present. */
    void erase(T const & t)
    {
      erase(find(t));
    }

    /** Remove all elements from the array. */
    void clear()
    {
      values.clear();
    }

}; // class SortedArray

} // namespace Thea

#endif
