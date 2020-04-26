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
// First version: 2012
//
//============================================================================

#ifndef __Thea_Concept_hpp__
#define __Thea_Concept_hpp__

#include "Common.hpp"

namespace Thea {

/** Define a template class that checks if the class supplied as the template argument defines a type with the given name. */
#define THEA_HAS_TYPE(check_name, type_name)                                                                                  \
template <typename check_name##T_>                                                                                            \
class check_name                                                                                                              \
{                                                                                                                             \
  private:                                                                                                                    \
    typedef char Yes;                                                                                                         \
    typedef struct { char a[2]; } No;                                                                                         \
    template <typename check_name##U_> static Yes test(typename check_name##U_::type_name *);                                 \
    template <typename check_name##U_> static No test(...);                                                                   \
                                                                                                                              \
  public:                                                                                                                     \
    static bool const value = (sizeof(test<check_name##T_>(0)) == sizeof(Yes));                                               \
};

#if __GNUC__ == 4 && ((__GNUC_MINOR__ == 4 && __GNUC_PATCHLEVEL__ > 2) || __GNUC_MINOR__ == 5)

#warning "Concept-checking for members disabled on 4.4.2 < GCC < 4.6"
#warning "    cf. http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46170"

#define THEA_HAS_MEMBER(check_name, member_name)                                                                              \
template <typename check_name##T_>                                                                                            \
class check_name                                                                                                              \
{                                                                                                                             \
  public:                                                                                                                     \
    static bool const value = true;                                                                                           \
};

#else

/**
 * Define a template class that checks if the class supplied as the template argument has a (data or function) member with the
 * given name.
 */
#define THEA_HAS_MEMBER(check_name, member_name)                                                                              \
template <typename check_name##T_>                                                                                            \
class check_name                                                                                                              \
{                                                                                                                             \
  private:                                                                                                                    \
    typedef char Yes;                                                                                                         \
    typedef struct { char a[2]; } No;                                                                                         \
    struct BaseMixin { int member_name; };                                                                                    \
    struct Base : public check_name##T_, public BaseMixin {};                                                                 \
    template <typename check_name##U_, check_name##U_> class Helper {};                                                       \
    template <typename check_name##U_> static No test(check_name##U_ *, Helper<int BaseMixin::*,                              \
                                                     & check_name##U_::member_name> * = 0);                                   \
    template <typename check_name##U_> static Yes test(...);                                                                  \
                                                                                                                              \
  public:                                                                                                                     \
    static bool const value = (sizeof(test<check_name##T_>((Base *)0)) == sizeof(Yes));                                       \
};

#endif

/** Utility function that directly tests the value of concept_name::value using THEA_STATIC_ASSERT. */
#define THEA_CONCEPT_CHECK(concept_name) static_assert(concept_name::value, "Assertion failed")

} // namespace Thea

#endif
