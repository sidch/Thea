//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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
#define THEA_CONCEPT_CHECK(concept_name) THEA_STATIC_ASSERT(concept_name::value)

} // namespace Thea

#endif
