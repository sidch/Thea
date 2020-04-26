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
// First version: 2019
//
//============================================================================

#ifndef __Thea_Hash_hpp__
#define __Thea_Hash_hpp__

#include <boost/functional/hash.hpp>
#include <utility>

namespace std {

// Specialize std::hash for std::pair types. For some inexplicable reason this is not present in the standard library.
template <typename U, typename V>
struct hash< pair<U, V> >
{
  size_t operator()(pair<U, V> const & p) const
  {
    return boost::hash_value(p);
  }
};

} // namespace std

#endif
