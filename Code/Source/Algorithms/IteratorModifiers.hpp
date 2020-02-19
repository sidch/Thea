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
#include <type_traits>

namespace Thea {
namespace Algorithms {

/**
 * Converts an iterator dereferencing to a pointer to T, to an iterator dereferencing to T. The new iterator supports the
 * operations (e.g. increment and comparison) of the original iterator. If the iterator does not dereference to a pointer, the
 * wrapper exactly mimics the behavior of the original iterator.
 */
template <typename IteratorT, typename Enable = void>
class RefIterator : public IteratorT
{
  private:
    /** The underlying value type (no top-level pointer). */
    typedef typename std::remove_pointer<
                typename std::remove_cv< typename std::iterator_traits<IteratorT>::value_type >::type >::type T;

  public:
    typedef T    value_type;  ///< Value type.
    typedef T *  pointer;     ///< Pointer to value.
    typedef T &  reference;   ///< Reference to value.

    /** Constructor. */
    RefIterator(IteratorT const & ii = IteratorT()) : IteratorT(ii) {}

    /** Dereferences the iterator to an object of type T. */
    T & operator*() const { return *(this->IteratorT::operator*()); }

    /** Arrow operator. */
    T * operator->() const { return this->IteratorT::operator*(); }

}; // class RefIterator

#define THEA_RANDOM_ACCESS_ITERATOR_BODY(class_name)                                                                          \
    class_name(class_name const & src) : ii(src.ii) {}                                                                        \
                                                                                                                              \
    class_name & operator++() { ++ii; return *this; }                                                                         \
    class_name   operator++(int) { class_name tmp(*this); operator++(); return tmp; }                                         \
    class_name & operator--() { --ii; return *this; }                                                                         \
    class_name   operator--(int) { class_name tmp(*this); operator--(); return tmp; }                                         \
                                                                                                                              \
    class_name & operator=(class_name const & src) { ii = src.ii; return *this; }                                             \
    template <typename IntegerT> class_name & operator+=(IntegerT n) { ii += n; return *this; }                               \
    template <typename IntegerT> class_name & operator-=(IntegerT n) { ii -= n; return *this; }                               \
                                                                                                                              \
    template <typename IntegerT> class_name operator+(IntegerT n) const { return class_name(ii + n); }                        \
    template <typename IntegerT> class_name operator-(IntegerT n) const { return class_name(ii - n); }                        \
                                                                                                                              \
    bool operator==(class_name const & rhs) const { return ii == rhs.ii; }                                                    \
    bool operator!=(class_name const & rhs) const { return ii != rhs.ii; }                                                    \
                                                                                                                              \
    bool operator< (class_name const & rhs) const { return ii <  rhs.ii; }                                                    \
    bool operator> (class_name const & rhs) const { return ii >  rhs.ii; }                                                    \
    bool operator<=(class_name const & rhs) const { return ii <= rhs.ii; }                                                    \
    bool operator>=(class_name const & rhs) const { return ii >= rhs.ii; }

// Specialization when the iterator does *NOT* dereference to a pointer.
template <typename IteratorT>
class RefIterator< IteratorT, typename std::enable_if<
                                  ! std::is_pointer<
                                        typename std::remove_cv< typename std::iterator_traits<IteratorT>::value_type >::type
                                    >::value
                               && ! std::is_pointer< typename std::remove_cv<IteratorT>::type >::value
                              >::type >
: public IteratorT
{
  public:
    RefIterator(IteratorT const & ii = IteratorT()) : IteratorT(ii) {}
};

// Specialization when the iterator is a pointer to T.
template <typename T>
class RefIterator< T *, typename std::enable_if< ! std::is_pointer< typename std::remove_cv<T>::type >::value >::type >
: public std::iterator<std::random_access_iterator_tag, T, std::ptrdiff_t, T *, T &>
{
  public:
    RefIterator(T * ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(RefIterator)

    T & operator*() const { return *ii; }
    T * operator->() const { return ii; }
    template <typename IntegerT> T & operator[](IntegerT n) { return ii[n]; }

  private:
    T * ii;
};

// Specialization when the iterator is a pointer to a pointer to T.
template <typename T>
class RefIterator<T **> : public std::iterator<std::random_access_iterator_tag, T, std::ptrdiff_t, T *, T &>
{
  public:
    explicit RefIterator(T ** ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(RefIterator)

    T & operator*() const { return **ii; }
    T * operator->() const { return *ii; }
    template <typename IntegerT> T & operator[](IntegerT n) { return *ii[n]; }

  private:
    T ** ii;

}; // class RefIterator<T, T **>

// Specialization when the iterator is a const pointer to a const pointer to T.
template <typename T>
class RefIterator<T const * const *>
: public std::iterator<std::random_access_iterator_tag, T const, std::ptrdiff_t, T const *, T const &>
{
  public:
    explicit RefIterator(T const * const * ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(RefIterator)

    T const & operator*() const { return **ii; }
    T const * operator->() const { return *ii; }
    template <typename IntegerT> T const & operator[](IntegerT n) { return *ii[n]; }

  private:
    T const * const * ii;

}; // class RefIterator<T const, T const * const *>

/**
 * Converts an iterator dereferencing to T, to an iterator dereferencing to T *. The new iterator supports the operations (e.g.
 * increment and comparison) of the original iterator. If the iterator already dereferences to a pointer, the wrapper exactly
 * mimics the behavior of the original iterator.
 *
 * @note Currently the implementation produces a read-only iterator because an iterator's <tt>value_type</tt> may not be
 * const-qualified even if it is a <tt>const_iterator</tt>, but a read-write version would probably not be very useful anyway.
 */
template <typename IteratorT, typename Enable = void>
class PtrIterator : public IteratorT
{
  private:
    /** The underlying value type. */
    typedef typename std::remove_cv< typename std::iterator_traits<IteratorT>::value_type >::type T;

  public:
    typedef T const *          value_type;  ///< Value type.
    typedef T const * const *  pointer;     ///< Pointer to value.
    typedef T const * const &  reference;   ///< Reference to value.

    /** Constructor. */
    PtrIterator(IteratorT const & ii = IteratorT()) : IteratorT(ii) {}

    /** Dereferences the iterator to an (immutable) object of type T. */
    T const * operator*() const { return &(this->IteratorT::operator*()); }

    /** Arrow operator is <b>DISABLED</b> (throws assertion error) since it cannot be called on an iterator-over-pointers. */
    T const * const * operator->() const { alwaysAssertM(false, "PtrIterator: Can't call '->' on iterator-over-pointers"); }

}; // class PtrIterator

// Specialization when the iterator already dereferences to a pointer.
template <typename IteratorT>
class PtrIterator< IteratorT, typename std::enable_if<
                                  std::is_pointer<
                                      typename std::remove_cv< typename std::iterator_traits<IteratorT>::value_type >::type
                                  >::value
                               && ! std::is_pointer< typename std::remove_cv<IteratorT>::type >::value
                              >::type >
: public IteratorT
{
  public:
    PtrIterator(IteratorT const & ii = IteratorT()) : IteratorT(ii) {}
};

// Specialization when the iterator is a pointer to T.
template <typename T>
class PtrIterator<T *>
: public std::iterator<std::random_access_iterator_tag, T, std::ptrdiff_t, T *, T &>
{
  public:
    PtrIterator(T * ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(PtrIterator)

    T * operator*() const { return ii; }
    template <typename IntegerT> T * operator[](IntegerT n) { return ii + n; }

  private:
    T * ii;
};

// Specialization when the iterator is a pointer to a pointer to T.
template <typename T>
class PtrIterator<T **>
: public std::iterator<std::random_access_iterator_tag, T *, std::ptrdiff_t, T **, T * &>
{
  public:
    PtrIterator(T ** ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(PtrIterator)

    T * operator*() const { return *ii; }
    T ** operator->() const { alwaysAssertM(false, "PtrIterator: Can't call '->' on iterator-over-pointers"); }
    template <typename IntegerT> T * operator[](IntegerT n) { return ii[n]; }

  private:
    T ** ii;
};

// Specialization when the iterator is a const pointer to a const pointer to T.
template <typename T>
class PtrIterator<T const * const *>
: public std::iterator<std::random_access_iterator_tag, T const *, std::ptrdiff_t, T const * const *, T const * const &>
{
  public:
    PtrIterator(T const * const * ii_ = NULL) : ii(ii_) {}

    THEA_RANDOM_ACCESS_ITERATOR_BODY(PtrIterator)

    T const * operator*() const { return *ii; }
    T const * const * operator->() const { alwaysAssertM(false, "PtrIterator: Can't call '->' on iterator-over-pointers"); }
    template <typename IntegerT> T const * operator[](IntegerT n) { return ii[n]; }

  private:
    T const * const * ii;
};

/** Utility function to create a specialization of RefIterator of the correct type wrapping a given iterator. */
template <typename IteratorT>
RefIterator<IteratorT>
makeRefIterator(IteratorT const & ii)
{
  return RefIterator<IteratorT>(ii);
}

/** Utility function to create a specialization of PtrIterator of the correct type wrapping a given iterator. */
template <typename IteratorT>
PtrIterator<IteratorT>
makePtrIterator(IteratorT const & ii)
{
  return PtrIterator<IteratorT>(ii);
}

#undef THEA_RANDOM_ACCESS_ITERATOR_BODY

} // namespace Algorithms
} // namespace Thea

#endif
