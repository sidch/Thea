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

#ifndef __Thea_Queue_hpp__
#define __Thea_Queue_hpp__

#include "Platform.hpp"
#include <deque>
#include <queue>

/** First-in-first-out queue of objects. */
template < typename T, typename Container = std::deque<T> > using Queue = std::queue<T, Container>;

#endif
