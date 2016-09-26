//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2010, Siddhartha Chaudhuri/Stanford University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#ifndef __Thea_Algorithms_FastCopy_hpp__
#define __Thea_Algorithms_FastCopy_hpp__

#include "../Common.hpp"
#include <boost/type_traits.hpp>
#include <cstring>
#include <iterator>

namespace Thea {
namespace Algorithms {

// Adapted from the Boost type_traits examples.

namespace FastCopyInternal {

template <typename I1, typename I2, bool b>
I2
fastCopyImpl(I1 first, I1 last, I2 out, boost::integral_constant<bool, b> const &)
{
  typedef typename std::iterator_traits<I2>::value_type value_type;
  while (first != last) *(out++) = static_cast<value_type>(*(first++));
  return out;
}

template <typename T>
T *
fastCopyImpl(T const * first, T const * last, T * out, boost::true_type const &)
{
   memcpy(out, first, (last - first) * sizeof(T));
   return out + (last - first);
}

template <typename I1, typename I2, bool b>
I2
fastCopyBackwardImpl(I1 first, I1 last, I2 out, boost::integral_constant<bool, b> const &)
{
  typedef typename std::iterator_traits<I2>::value_type value_type;
  while (last != first) *(--out) = static_cast<value_type>(*(--last));
  return out;
}

template <typename T>
T *
fastCopyBackwardImpl(T const * first, T const * last, T * out, boost::true_type const &)
{
   memmove(out, first, (last - first) * sizeof(T));
   return out;
}

} // namespace FastCopyInternal

/**
 * A version of <tt>std::copy</tt> that calls <tt>memcpy</tt> where appropriate (if the class has a trivial assignment operator
 * and the iterators are raw pointers) for speed.
 *
 * To take advantage of fast copying, specialize boost::has_trivial_assign to return true for the value type.
 */
template <typename I1, typename I2>
inline I2
fastCopy(I1 first, I1 last, I2 out)
{
  //
  // We can copy with memcpy if T has a trivial assignment operator,
  // and if the iterator arguments are actually pointers (this last
  // requirement we detect with overload resolution):
  //
  typedef typename std::iterator_traits<I1>::value_type value_type;
  return FastCopyInternal::fastCopyImpl(first, last, out, boost::has_trivial_assign<value_type>());
}

/**
 * A version of <tt>std::copy_backward</tt> that calls <tt>memmove</tt> where appropriate (if the class has a trivial assignment
 * operator) for speed.
 *
 * To take advantage of fast copying, specialize boost::has_trivial_assign to return true for the value type.
 */
template <typename I1, typename I2>
inline I2
fastCopyBackward(I1 first, I1 last, I2 out)
{
  typedef typename std::iterator_traits<I1>::value_type value_type;
  return FastCopyInternal::fastCopyBackwardImpl(first, last, out, boost::has_trivial_assign<value_type>());
}

} // namespace Algorithms
} // namespace Thea

#endif
