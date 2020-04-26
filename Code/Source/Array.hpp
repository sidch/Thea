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

#ifndef __Thea_Array_hpp__
#define __Thea_Array_hpp__

#include "Platform.hpp"
#include "AlignedAllocator.hpp"
#include <vector>

namespace Thea {

/** Dynamically resizable array. */
template < typename T, typename Alloc = std::allocator<T> > using Array = std::vector<T, Alloc>;

/** Dynamically resizable array with aligned memory allocation. */
template < typename T, size_t N = 16 > using AlignedArray = std::vector< T, AlignedAllocator<T, N> >;

} // namespace Thea

#endif
