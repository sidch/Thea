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

#ifndef __Thea_Common_hpp__
#define __Thea_Common_hpp__

#include "Platform.hpp"
#include "CommonEnums.hpp"
#include "EnumClass.hpp"
#include "Error.hpp"
#include "Hash.hpp"
#include "Memory.hpp"
#include "Noncopyable.hpp"
#include "NumericType.hpp"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#if defined(__GNUC__) && __GNUC__ >= 3
#  include <cxxabi.h>
#  define THEA_HAVE_CXA_DEMANGLE
#endif

// Visual Studio requires templates instantiations to be explicitly imported from DLL's to avoid conflicts like
// http://www.codesynthesis.com/~boris/blog/2010/01/18/dll-export-cxx-templates/
#ifdef _MSC_VER
#  define THEA_EXPORT_INSTANTIATION
#endif

/** Root namespace for the %Thea library. */
namespace Thea {

/** Check if a test condition is true, and immediately abort the program with an error code if not, in debug mode only. */
#ifdef THEA_DEBUG_BUILD

template <typename CondT, typename MessageT>
inline void
debugAssertM(CondT const & test, MessageT const & msg)
{
  if (!test)
  {
    std::cerr << "!!! Debug-mode assertion failed !!! " << msg << std::endl;
    std::exit(-1);
  }
}

#else
#  define debugAssertM(test, msg) (void)0
#endif

/** Check if a test condition is true, and immediately abort the program with an error code if not. */
template <typename CondT, typename MessageT>
inline void
alwaysAssertM(CondT const & test, MessageT const & msg)
{
  if (!test)
  {
    std::cerr << "!!! Assertion failed !!! " << msg << std::endl;
    std::exit(-1);
  }
}

// GCC doesn't demangle type_info::name(). New versions include a demangling function in the ABI.
#ifdef THEA_HAVE_CXA_DEMANGLE

namespace CommonInternal {

inline std::string
demangle(char const * class_name)
{
  size_t length = 1024;
  int status;

  char * ret = abi::__cxa_demangle(class_name, nullptr, &length, &status);
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
}

} // namespace CommonInternal

#endif // THEA_HAVE_CXA_DEMANGLE

/** Get the name of a type, ignoring reference and cv-qualifiers. */
template <typename T>
std::string
getTypeName()
{
#ifdef THEA_HAVE_CXA_DEMANGLE
  return CommonInternal::demangle(typeid(T).name());
#else
  return typeid(T).name();
#endif
}

/** Get the name of the dynamic type of an object, ignoring reference and cv-qualifiers. */
template <typename T>
std::string
getTypeName(T const & obj)
{
#ifdef THEA_HAVE_CXA_DEMANGLE
  return CommonInternal::demangle(typeid(obj).name());
#else
  return typeid(obj).name();
#endif
}

} // namespace Thea

// Commonly-used but requires the stuff in Common.hpp to be declared first, so is included here.
#include "StringAlg.hpp"

#endif
