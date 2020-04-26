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

#ifndef __Thea_UnorderedMap_hpp__
#define __Thea_UnorderedMap_hpp__

#include "Platform.hpp"
#include <unordered_map>

namespace Thea {

/** Hash table-based mapping from keys to values. */
template < typename Key,
           typename T,
           typename Hash = std::hash<Key>,
           typename Pred = std::equal_to<Key>,
           typename Alloc = std::allocator< std::pair<Key const, T> >
         > using UnorderedMap = std::unordered_map<Key, T, Hash, Pred, Alloc>;

/** Hash table-based mapping from (possibly duplicate) keys to values. */
template < typename Key,
           typename T,
           typename Hash = std::hash<Key>,
           typename Pred = std::equal_to<Key>,
           typename Alloc = std::allocator< std::pair<Key const, T> >
         > using UnorderedMultiMap = std::unordered_multimap<Key, T, Hash, Pred, Alloc>;

} // namespace Thea

#endif
