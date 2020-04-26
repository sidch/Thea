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
// First version: 2009
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
