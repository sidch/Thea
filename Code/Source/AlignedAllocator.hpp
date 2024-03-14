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
// First version: 2014
//
//============================================================================

#ifndef __Thea_AlignedAllocator_hpp__
#define __Thea_AlignedAllocator_hpp__

#include "Platform.hpp"
#include <cstdlib>

namespace Thea {

/**
 * Allocates aligned memory blocks for one or more objects of type `T`. `Alignment` should be the desired alignment in bytes,
 * e.g. 1, 4, 8 or 16.
 *
 * Originally from http://stackoverflow.com/a/8545389 but now updated to use C++17 std::aligned_alloc.
 */
template <typename T, size_t Alignment = 16>
class AlignedAllocator
{
  private:
    static_assert(Alignment >= 1, "AlignedAllocator: Alignment must be non-negative");

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
    AlignedAllocator() throw () {}

    /** Copy constructor. */
    template <typename T2> AlignedAllocator(AlignedAllocator<T2, Alignment> const &) throw () {}

    /** Destructor. */
    ~AlignedAllocator() throw () {}

    /** Get the address of a referenced object. */
    pointer address(reference r) { return &r; }

    /** Get the address of a referenced object. */
    const_pointer address(const_reference r) const { return &r; }

    /** Allocate an aligned block of \a count objects. */
    pointer allocate(size_type count) { return (pointer)std::aligned_alloc(Alignment, count * sizeof(value_type)); }

    /** Deallocate an aligned block. */
    void deallocate(pointer p, size_type count) { (void)count; std::free(p); }

    /** Construct an object at a memory location. */
    void construct(pointer p, const value_type & val) { new (p) value_type(val); }

    /** Destroy an object at a memory location. */
    void destroy(pointer p) { p->~value_type (); }

    /** Get the maximum number of elements that can theoretically be allocated. */
    size_type max_size() const throw () { return size_type(-1) / sizeof(value_type); }

    /** A structure that enables this allocator to allocate storage for objects of another type. */
    template <typename T2>
    struct rebind
    {
      typedef AlignedAllocator<T2, Alignment> other;
    };

    /** Check if two allocators are different. */
    bool operator!=(const AlignedAllocator<T, Alignment> & other) const { return !(*this == other); }

    /**
     * Check if two allocators are the same. Returns true if and only if storage allocated from *this can be deallocated from
     * \a other, and vice versa. Always returns true for stateless allocators like this one.
     */
    bool operator==(const AlignedAllocator<T, Alignment> & other) const { return true; }

}; // class AlignedAllocator

} // namespace Thea

#endif
