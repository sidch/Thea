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

#ifndef __Thea_BoundedArrayN_hpp__
#define __Thea_BoundedArrayN_hpp__

#include "Common.hpp"
#include "Algorithms/FastCopy.hpp"

namespace Thea {

/**
 * An array of bounded maximum size N, stored on the stack. No new elements can be added to the end of the array once it reaches
 * its maximum capacity, until some are removed via erase(). Elements can be inserted in the middle of the array even when it is
 * full: in this case the last element is dropped to make space. This class is useful for very fast allocation of space for a
 * few elements, where the exact number of elements is not known but is guaranteed to have an upper limit.
 *
 * To get some extra speed when T has a trivial (bit-copy) assignment operator, make sure that
 * <tt>boost::has_trivial_assign</tt> is true for T.
 *
 * The implementation always allocates enough space to store the maximum number of instances of T. The capacity N should be
 * <b>positive</b> (non-zero).
 *
 * For a bounded stack-based array that also guarantees that its elements are in sorted order, see BoundedSortedArrayN.
 */
template <int N, typename T>
class BoundedArrayN
{
  private:
    int num_elems;
    T values[N];

  public:
    /** Constructor. */
    BoundedArrayN() : num_elems(0) {}

    /** Get the maximum number of elements the array can hold. */
    static int getCapacity() { return N; }

    /** Get the number of elements in the array. */
    int size() const { return num_elems; }

    /** Check if the array is empty or not. */
    bool isEmpty() const { return num_elems <= 0; }

    /** Check if the array has reached maximum capacity or not. */
    bool isFull() const { return num_elems >= N; }

    /** Get the first element in the array. */
    T const & first() const
    {
      debugAssertM(num_elems > 0, "BoundedArrayN: Can't get first element of empty array");
      return values[0];
    }

    /** Get the last element in the array. */
    T const & last() const
    {
      debugAssertM(num_elems > 0, "BoundedArrayN: Can't get last element of empty array");
      return values[num_elems - 1];
    }

    /** Get the element at a given position in the array. Bounds checks are only performed in debug mode. */
    T const & operator[](int i) const
    {
      debugAssertM(i >= 0 && i < num_elems, format("BoundedArrayN: Index %d out of bounds [0, %d)", i, num_elems));
      return values[i];
    }

    /** Get the element at a given position in the array. Bounds checks are only performed in debug mode. */
    T & operator[](int i)
    {
      debugAssertM(i >= 0 && i < num_elems, format("BoundedArrayN: Index %d out of bounds [0, %d)", i, num_elems));
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
    void insert(int i, T const & t)
    {
      debugAssertM(i >= 0 && i <= std::min(num_elems, N - 1), format("BoundedArrayN: Index %d out of bounds [0, %d]", i,
                                                                     std::min(num_elems, N - 1)));

      if (isFull())
        Algorithms::fastCopyBackward(values + i, values + N - 1, values + i + 1);
      else if (i < num_elems)
        Algorithms::fastCopyBackward(values + i, values + num_elems, values + i + 1);

      values[i] = t;
    }

    /** Remove the element at the given position from the array. */
    void erase(int i)
    {
      if (i >= 0 && i < num_elems)
      {
        Algorithms::fastCopy(values + i + 1, values + num_elems, values + i);
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
