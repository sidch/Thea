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
#include <memory>

namespace Thea {

/**
 * Define smart pointer types for a class.
 */
#define THEA_DECL_SMART_POINTERS(type)                                                                                        \
    typedef std::shared_ptr< type > Ptr;                                                                                      \
    typedef std::shared_ptr< type const > ConstPtr;                                                                           \
    typedef std::weak_ptr< type > WeakPtr;                                                                                    \
    typedef std::weak_ptr< type const > ConstWeakPtr;

/**
 * Macro that declares the standard smart pointer types as 'extern'.
 */
#ifdef THEA_EXTERN_TEMPLATES
#  define THEA_DECL_EXTERN_SMART_POINTERS(class_name)                                                                         \
     namespace std {                                                                                                          \
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
     namespace std {                                                                                                          \
       template class shared_ptr<class_name>;                                                                                 \
       template class shared_ptr<class_name const>;                                                                           \
       template class weak_ptr<class_name>;                                                                                   \
       template class weak_ptr<class_name const>;                                                                             \
     }
#else
#  define THEA_INSTANTIATE_SMART_POINTERS(class_name)  /* blank */
#endif // THEA_EXTERN_TEMPLATES

} // namespace Thea

#endif
