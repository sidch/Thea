//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_IteratorModifiers_hpp__
#define __Thea_Algorithms_IteratorModifiers_hpp__

#include "../Common.hpp"
#include <iterator>

namespace Thea {
namespace Algorithms {

/**
 * Converts an iterator dereferencing to a pointer to T, to an iterator dereferencing to T. The new iterator supports the
 * operations (e.g. increment and comparison) of the original iterator.
 */
template <typename T, typename PtrIterator>
class PtrToRefIterator : public PtrIterator
{
  public:
    /** Constructor. */
    PtrToRefIterator(PtrIterator const & ii = PtrIterator()) : PtrIterator(ii) {}

    /** Dereferences the iterator to an object of type T. */
    T & operator*() const { return *(this->PtrIterator::operator*()); }

}; // class PtrToRefIterator

#define THEA_PTR_PTR_TO_REF_ITERATOR_BODY(T)                                                                                  \
    PtrToRefIterator(PtrToRefIterator const & src) : ii(src.ii) {}                                                            \
                                                                                                                              \
    T & operator*() const { return **ii; }                                                                                    \
    template <typename IntegerT> T & operator[](IntegerT n) { return *ii[n]; }                                                \
                                                                                                                              \
    PtrToRefIterator & operator++() { ++ii; return *this; }                                                                   \
    PtrToRefIterator   operator++(int) { PtrToRefIterator tmp(*this); operator++(); return tmp; }                             \
    PtrToRefIterator & operator--() { --ii; return *this; }                                                                   \
    PtrToRefIterator   operator--(int) { PtrToRefIterator tmp(*this); operator--(); return tmp; }                             \
                                                                                                                              \
    PtrToRefIterator & operator=(PtrToRefIterator const & src) { ii = src.ii; return *this; }                                 \
    template <typename IntegerT> PtrToRefIterator & operator+=(IntegerT n) { ii += n; return *this; }                         \
    template <typename IntegerT> PtrToRefIterator & operator-=(IntegerT n) { ii -= n; return *this; }                         \
                                                                                                                              \
    template <typename IntegerT> PtrToRefIterator operator+(IntegerT n) const { return PtrToRefIterator(ii + n); }            \
    template <typename IntegerT> PtrToRefIterator operator-(IntegerT n) const { return PtrToRefIterator(ii - n); }            \
                                                                                                                              \
    bool operator==(PtrToRefIterator const & rhs) const { return ii == rhs.ii; }                                              \
    bool operator!=(PtrToRefIterator const & rhs) const { return ii != rhs.ii; }                                              \
                                                                                                                              \
    bool operator< (PtrToRefIterator const & rhs) const { return ii <  rhs.ii; }                                              \
    bool operator> (PtrToRefIterator const & rhs) const { return ii >  rhs.ii; }                                              \
    bool operator<=(PtrToRefIterator const & rhs) const { return ii <= rhs.ii; }                                              \
    bool operator>=(PtrToRefIterator const & rhs) const { return ii >= rhs.ii; }

// Specialization when the iterator is a pointer to a pointer to T.
template <typename T>
class PtrToRefIterator<T, T **> : public std::iterator<std::random_access_iterator_tag, T, std::ptrdiff_t, T *, T &>
{
  public:
    explicit PtrToRefIterator(T ** ii_ = NULL) : ii(ii_) {}

    THEA_PTR_PTR_TO_REF_ITERATOR_BODY(T)

  private:
    T ** ii;

}; // class PtrToRefIterator<T, T **>

// Specialization when the iterator is a const pointer to a const pointer to T.
template <typename T>
class PtrToRefIterator<T const, T const * const *> : public std::iterator<std::random_access_iterator_tag,
                                                                          T const, std::ptrdiff_t, T const *, T const &>
{
  public:
    explicit PtrToRefIterator(T const * const * ii_ = NULL) : ii(ii_) {}

    THEA_PTR_PTR_TO_REF_ITERATOR_BODY(T const)

  private:
    T const * const * ii;

}; // class PtrToRefIterator<T const, T const * const *>

#undef THEA_PTR_PTR_TO_REF_ITERATOR_BODY

} // namespace Algorithms
} // namespace Thea

#endif
