//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, * except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Cornell University
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

#ifndef __Thea_AlignedAllocator_hpp__
#define __Thea_AlignedAllocator_hpp__

#include "Common.hpp"

#ifdef THEA_WINDOWS
#  include <malloc.h>
#else
#  include <stdlib.h>
#endif

namespace Thea {

/**
 * Allocates aligned memory blocks.
 *
 * From http://stackoverflow.com/a/8545389
 */
template <typename T, size_t N = 16>
class AlignedAllocator
{
  public:
    typedef T value_type;                    ///< Type of allocated objects.
    typedef std::size_t size_type;           ///< Type of block sizes.
    typedef std::ptrdiff_t difference_type;  ///< Type of difference of two pointers.

    typedef T * pointer;                     ///< Pointer to an allocated object/block.
    typedef T const * const_pointer;         ///< Const pointer to an allocated object/block.

    typedef T & reference;                   ///< Reference to an allocated object.
    typedef T const & const_reference;       ///< Const reference to an allocated object.

  public:
    /** Default Constructor. */
    inline AlignedAllocator() throw () {}

    /** Copy constructor. */
    template <typename T2> inline AlignedAllocator(AlignedAllocator<T2, N> const &) throw () {}

    /** Destructor. */
    inline ~AlignedAllocator() throw () {}

    /** Get the address of a referenced object. */
    inline pointer address(reference r) { return &r; }

    /** Get the address of a referenced object. */
    inline const_pointer address(const_reference r) const { return &r; }

    /** Allocate an aligned block of \a n bytes. */
    inline pointer allocate(size_type n)
    {
#ifdef THEA_WINDOWS
      return (pointer)_aligned_malloc(n * sizeof(value_type), N);
#else
      pointer p;
      if (posix_memalign((void **)(&p), N, n * sizeof(value_type)) == 0)
        return p;
      else
        return NULL;
#endif
    }

    /** Deallocate an aligned block. */
    inline void deallocate(pointer p, size_type n = 0)
    {
#ifdef THEA_WINDOWS
      _aligned_free(p);
#else
      free(p);
#endif
    }

    /** Construct an object at a memory location. */
    inline void construct (pointer p, const value_type & wert) { new (p) value_type (wert); }

    /** Destroy an object at a memory location. */
    inline void destroy (pointer p) { p->~value_type (); }

    /** Get the maximum number of elements that can theoretically be allocated. */
    inline size_type max_size () const throw () { return size_type (-1) / sizeof (value_type); }

    /** A structure that enables this allocator to allocate storage for objects of another type. */
    template <typename T2>
    struct rebind
    {
      typedef AlignedAllocator<T2, N> other;
    };

    /** Check if two allocators are different. */
    bool operator!=(const AlignedAllocator<T, N> & other) const  { return !(*this == other); }

    /**
     * Check if two allocators are the same. Returns true if and only if storage allocated from *this can be deallocated from
     * \a other, and vice versa. Always returns true for stateless allocators like this one.
     */
    bool operator==(const AlignedAllocator<T, N> & other) const { return true; }

}; // class AlignedAllocator

} // namespace Thea

#endif
