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

#ifndef __Thea_Common_hpp__
#define __Thea_Common_hpp__

#include "Platform.hpp"
#include "AtomicInt32.hpp"
#include "CommonEnums.hpp"
#include "EnumClass.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "Memory.hpp"
#include "Noncopyable.hpp"
#include "NumericTypes.hpp"
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#if defined(__GNUC__) && __GNUC__ >= 3
#  include <cxxabi.h>
#  define THEA_HAVE_CXA_DEMANGLE
#endif

#undef debugAssertM
#ifdef THEA_DEBUG_BUILD
#  define debugAssertM(test, msg)                                                                                             \
   {                                                                                                                          \
     if (!(test))                                                                                                             \
     {                                                                                                                        \
       THEA_ERROR << "Debug-mode assertion failed: " << (msg);                                                                \
       throw Thea::FatalError(msg);                                                                                           \
     }                                                                                                                        \
   }
#else
#  define debugAssertM(test, msg) {}
#endif

#define alwaysAssertM(test, msg)                                                                                              \
   {                                                                                                                          \
     if (!(test))                                                                                                             \
     {                                                                                                                        \
       THEA_ERROR << "Assertion failed: " << (msg);                                                                           \
       throw Thea::FatalError(msg);                                                                                           \
     }                                                                                                                        \
   }

// Visual Studio requires templates instantiations to be explicitly imported from DLL's to avoid conflicts like
// http://www.codesynthesis.com/~boris/blog/2010/01/18/dll-export-cxx-templates/
#ifdef _MSC_VER
#  define THEA_EXPORT_INSTANTIATION
#endif

/** Root namespace for the %Thea library. */
namespace Thea {

/**
 * Require an expression to evaluate to true at compile-time. Example usage:
 *
 * \code
 *   THEA_STATIC_ASSERT(sizeof(int) == 2)
 * \endcode
 *
 * From Ralf Holly, http://www.drdobbs.com/compile-time-assertions/184401873
 */
#define THEA_STATIC_ASSERT(e) \
do \
{ \
  enum { assert_static__ = 1/((int)(e)) }; \
} while (0)

/** Get the class of an object. */
template <typename T>
/* THEA_API */
std::string
getClass(T const & obj)
{
  // GCC doesn't demangle type_info::name(). New versions include a demangling function in the ABI.
#ifdef THEA_HAVE_CXA_DEMANGLE

  size_t length = 1024;
  int status;

  char const * class_name = typeid(obj).name();
  char * ret = abi::__cxa_demangle(class_name, NULL, &length, &status);
  if (ret)
  {
    std::string result(ret);
    std::free(ret);
    return result;
  }
  else
  {
    std::free(ret);
    return class_name;
  }

#else

  return typeid(obj).name();

#endif
}

} // namespace Thea

// Commonly-used but requires the stuff in Common.hpp to be declared first, so is included here.
#include "StringAlg.hpp"

#endif
