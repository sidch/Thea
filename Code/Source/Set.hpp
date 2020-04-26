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

#ifndef __Thea_Set_hpp__
#define __Thea_Set_hpp__

#include "Platform.hpp"
#include <set>

namespace Thea {

/** Set of objects. Requires an ordering on the objects. */
template < typename T,
           typename Compare = std::less<T>,
           typename Alloc = std::allocator<T>
         > using Set = std::set<T, Compare, Alloc>;

/** Set of objects, with possible duplication. Requires an ordering on the objects. */
template < typename T,
           typename Compare = std::less<T>,
           typename Alloc = std::allocator<T>
         > using MultiSet = std::multiset<T, Compare, Alloc>;

} // namespace Thea

#endif
