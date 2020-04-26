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

#ifndef __Thea_UnorderedSet_hpp__
#define __Thea_UnorderedSet_hpp__

#include "Platform.hpp"
#include <unordered_set>

namespace Thea {

/** Hash table-based set of objects. */
template < typename T,
           typename Hash = std::hash<T>,
           typename Pred = std::equal_to<T>,
           typename Alloc = std::allocator<T>
         > using UnorderedSet = std::unordered_set<T, Hash, Pred, Alloc>;

/** Hash table-based set of objects, with possible duplication. */
template < typename T,
           typename Hash = std::hash<T>,
           typename Pred = std::equal_to<T>,
           typename Alloc = std::allocator<T>
         > using UnorderedMultiSet = std::unordered_multiset<T, Hash, Pred, Alloc>;

} // namespace Thea

#endif
