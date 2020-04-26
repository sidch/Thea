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

#ifndef __Thea_Map_hpp__
#define __Thea_Map_hpp__

#include "Platform.hpp"
#include <map>

namespace Thea {

/** Maps keys to values. Requires an ordering on the keys. */
template < typename Key,
           typename T,
           typename Compare = std::less<Key>,
           typename Alloc = std::allocator< std::pair<Key const, T> >
         > using Map = std::map<Key, T, Compare, Alloc>;

/** Maps (possibly duplicate) keys to values. Requires an ordering on the keys. */
template < typename Key,
           typename T,
           typename Compare = std::less<Key>,
           typename Alloc = std::allocator< std::pair<Key const, T> >
         > using MultiMap = std::multimap<Key, T, Compare, Alloc>;

} // namespace Thea

#endif
