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

#ifndef __Thea_Memory_hpp__
#define __Thea_Memory_hpp__

#include "Platform.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/functional/hash.hpp>
#include <boost/version.hpp>

namespace Thea {

/**
 * Smart pointer types. These <b><em>MUST</em> be thread-safe</b> -- we assume thread-safety throughout.
 */
using boost::shared_ptr;
using boost::shared_array;
using boost::weak_ptr;
using boost::static_pointer_cast;
using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;
using boost::enable_shared_from_this;
using boost::swap;

/**
 * Generic constructors similar to those of shared_ptr. Saves typing when you subclass shared_ptr. E.g. if you create a subclass
 * CPtr from shared_ptr<%C> to add extra functions specific to class C, use
 * <code>THEA_DEF_SHARED_PTR_CONSTRUCTORS(CPtr, shared_ptr<%C>)</code> in the class declaration to generate the standard
 * constructors.
 */
#define THEA_DEF_SHARED_PTR_CONSTRUCTORS(type, base_type) \
    type() : base_type() {} \
    template<class Y> explicit type(Y * p) : base_type(p) {} \
    template<class Y, class D> type(Y * p, D d) : base_type(p, d) {} \
    template<class Y, class D, class A> type(Y * p, D d, A a) \
    : base_type(p, d, a) {} \
    template<class Y> type(shared_ptr<Y> const & r) : base_type(r) {} \
    template<class Y> explicit type(weak_ptr<Y> const & r) \
    : base_type(r) {} \
    template<class Y> explicit type(std::auto_ptr<Y> & r) : base_type(r) {}

/**
 * Generic constructors similar to those of shared_array. Saves typing when you subclass shared_array.
 *
 * @see THEA_DEF_SHARED_PTR_CONSTRUCTORS
 */
#define THEA_DEF_SHARED_ARRAY_CONSTRUCTORS(type, base_type, obj_type) \
    explicit type(obj_type * p = 0) : base_type(p) {} \
    template<class D> type(obj_type * p, D d) : base_type(p, d) {} \
    type(type const & r) : base_type(r) {}

/**
 * Generic constructors similar to those of weak_ptr. Saves typing when you subclass weak_ptr.
 *
 * @see THEA_DEF_SHARED_PTR_CONSTRUCTORS
 */
#define THEA_DEF_WEAK_PTR_CONSTRUCTORS(type, base_type) \
    type() {} \
    template<class Y> type(shared_ptr<Y> const & r) : base_type(r) {} \
    type(type const & r) : base_type(r) {} \
    template<class Y> type(weak_ptr<Y> const & r) : base_type(r) {}

/**
 * Define smart pointer types for a class.
 */
#define THEA_DEF_POINTER_TYPES(type, shared_pointer, weak_pointer) \
    typedef shared_pointer< type > Ptr; \
    typedef shared_pointer< type const > ConstPtr; \
    typedef weak_pointer< type > WeakPtr; \
    typedef weak_pointer< type const > ConstWeakPtr;

/**
 * Macro that declares the standard smart pointer types as 'extern'.
 */
#ifdef THEA_EXTERN_TEMPLATES
#  define THEA_DECL_EXTERN_SMART_POINTERS(class_name)                                                                         \
     namespace boost {                                                                                                        \
       extern template class shared_ptr<class_name>;                                                                          \
       extern template class shared_ptr<class_name const>;                                                                    \
       extern template class weak_ptr<class_name>;                                                                            \
       extern template class weak_ptr<class_name const>;                                                                      \
     }
#else
#  define THEA_DECL_EXTERN_SMART_POINTERS(class_name)  /* blank */
#endif // THEA_EXTERN_TEMPLATES

/**
 * Macro that explicitly instantiates the standard smart pointer types
 */
#ifdef THEA_EXTERN_TEMPLATES
#  define THEA_INSTANTIATE_SMART_POINTERS(class_name)                                                                         \
     namespace boost {                                                                                                        \
       template class shared_ptr<class_name>;                                                                                 \
       template class shared_ptr<class_name const>;                                                                           \
       template class weak_ptr<class_name>;                                                                                   \
       template class weak_ptr<class_name const>;                                                                             \
     }
#else
#  define THEA_INSTANTIATE_SMART_POINTERS(class_name)  /* blank */
#endif // THEA_EXTERN_TEMPLATES

} // namespace Thea

// Boost 1.47.0 first defined hash_value for shared pointers
#if BOOST_VERSION < 104700

namespace boost {

template <class T>
std::size_t
hash_value(shared_ptr<T> const & ptr)
{
  return hash_value(ptr.get());
}

template <class T>
std::size_t
hash_value(shared_array<T> const & arr)
{
  return hash_value(arr.get());
}

} // namespace boost

#endif

#endif
