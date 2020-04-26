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

#ifndef __Thea_List_hpp__
#define __Thea_List_hpp__

#include "Platform.hpp"
#include <list>

namespace Thea {

/** Linked list. */
template < typename T, typename Alloc = std::allocator<T> > using List = std::list<T, Alloc>;

} // namespace Thea

#endif
