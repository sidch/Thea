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

#include "Platform.hpp"
#include <functional>
#include <utility>

namespace Thea {

/** Custom hasher that extends `std::hash`. The default specialization just duplicates `std::hash`. */
template <typename T, typename Enable = void>
struct Hasher : public std::hash<T> {};

/** Convenience function that returns the hash of any object supported by Hasher. */
template <typename T>
std::size_t
hashValue(T const & v)
{
  return Hasher<T>()(v);
}

/** Hash a value \a v and combine it with a preexisting hash \a seed. The combination is <b>order-dependent</b>. */
template <class T>
void
hashCombine(std::size_t & seed, T const & v)
{
  // Copied from boost::hash_combine
  seed ^= hashValue(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); // magic random number ensures spreading of hashes
}

/** Hash a range of values, optionally starting with a pre-existing hash \a seed. */
template <typename Iterator>
std::size_t
hashRange(Iterator first, Iterator last, std::size_t seed = 0)
{
  // Copied from boost::hash_range
  std::size_t h = seed;
  for(; first != last; ++first)
    hashCombine(seed, *first);

  return h;
}

/** Specialization of Hasher for `std::pair` types. For some inexplicable reason this is not present in the standard library. */
template <typename U, typename V>
struct Hasher< std::pair<U, V> >
{
  std::size_t operator()(std::pair<U, V> const & p) const
  {
    std::size_t h = hashValue(p.first);
    hashCombine(h, p.second);
    return h;
  }
};

} // namespace Thea

#endif
