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

/**
 * \file
 * This file <b>must</b> be included in <em>every</em> file in the project, <em>before</em> any other include.
 */

#ifndef __Thea_Platform_hpp__
#define __Thea_Platform_hpp__

#include <cstddef>

#define WIN32_LEAN_AND_MEAN

#if !defined(THEA_DEBUG_BUILD) && (!defined(NDEBUG) || defined(_DEBUG))
#  define THEA_DEBUG_BUILD
#endif

// NDEBUG needed for assert to be deactivated on release builds
#if !defined(THEA_DEBUG_BUILD) && !defined(NDEBUG)
#  define NDEBUG
#endif

#ifdef _MSC_VER
#  define THEA_WINDOWS 1
   // NOMINMAX required to stop windows.h messing up std::min
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif

   // Disable pesky warnings of the type: "class <member-classname> needs to have dll-interface to be used by clients of class"
   // and "non dll-interface class <subclass> used as base for dll-interface class <superclass>."
#  pragma warning( disable: 4251 )
#  pragma warning( disable: 4275 )

   // Disable numerical type conversion "possible loss of precision" warnings
#  pragma warning( disable: 4244 )

#elif defined(__FreeBSD__)
#  define THEA_FREEBSD 1
#elif defined(__OpenBSD__)
#  define THEA_OPENBSD 1
#elif defined(__linux__)
#  define THEA_LINUX 1
#elif defined(__APPLE__)
   // Retain old OSX tag for backwards-compatibility
#  define THEA_MAC 1
#  define THEA_OSX 1
   // Prevent OS X fp.h header from being included; it defines pi as a constant, which creates a conflict with G3D
#  define __FP__
#else
#  error Unknown platform
#endif

#if defined(THEA_FREEBSD) || defined(THEA_OPENBSD)
#  define THEA_BSD 1
#endif

// http://nadeausoftware.com/articles/2012/02/c_c_tip_how_detect_processor_type_using_compiler_predefined_macros
// http://sourceforge.net/p/predef/wiki/Architectures/
#if defined(_M_X64) || defined(__x86_64__) || defined(__amd_64__)
#  define THEA_64BIT 1
#  define THEA_X64 1
#elif defined(_M_IA64) || defined(__ia64) || defined(__ia64__)
#  define THEA_64BIT 1
#  define THEA_ITANIUM 1
#elif defined(__i386) || defined(_M_IX86)
#  define THEA_32BIT 1
#  define THEA_X86 1
#elif defined(__aarch64__) || defined(__aarch64)
#  define THEA_64BIT 1
#  define THEA_ARM 1
#elif defined(_M_ARM) || defined(__arm) || defined(__arm__)
#  define THEA_32BIT 1
#  define THEA_ARM 1
#elif defined(__powerpc__) || defined(__PPC__) || defined(__POWERPC__)
#  if defined(__ppc64__) || defined(__powerpc64__) || defined(__64BIT__) || defined(_LP64)
#    define THEA_64BIT 1
#  else
#    define THEA_32BIT 1
#  endif
#  define THEA_POWERPC 1
#endif

// Symbol visibility flags
#include "SymbolVisibility.hpp"

// Calling conventions
#if defined(_MSC_VER) && defined(THEA_32BIT)  // also clang and icc on Windows?
#  if defined(__GNUC__)
#    define THEA_CDECL __attribute__((cdecl))
#    define THEA_STDCALL __attribute__((stdcall))
#  else
#    define THEA_CDECL __cdecl
#    define THEA_STDCALL __stdcall
#  endif
#else
#  define THEA_CDECL
#  define THEA_STDCALL
#endif

/** Calling convention for member functions of abstract interface classes (for talking to shared libraries). */
#define THEA_ICALL THEA_STDCALL

/**
 * @def THEA_PRAGMA(expression)
 * \#pragma may not appear inside a macro, so this uses the pragma operator to create an equivalent statement.
 */
#ifdef _MSC_VER
// Microsoft's version http://msdn.microsoft.com/en-us/library/d9x1s805.aspx
#  define THEA_PRAGMA(x) __pragma(x)
#else
// C99 standard http://www.delorie.com/gnu/docs/gcc/cpp_45.html
#  define THEA_PRAGMA(x) _Pragma(#x)
#endif

/** Use tight alignment for a class. Useful for small classes that must be packed into an array. */
#ifdef _MSC_VER
#  define THEA_BEGIN_PACKED_CLASS(byte_align)  THEA_PRAGMA( pack(push, byte_align) )
#else
#  define THEA_BEGIN_PACKED_CLASS(byte_align)
#endif

/** Mark end of class that uses tight alignment. */
#ifdef _MSC_VER
#  define THEA_END_PACKED_CLASS(byte_align)  ; THEA_PRAGMA( pack(pop) )
#elif defined(__GNUC__)
#  define THEA_END_PACKED_CLASS(byte_align)  __attribute__((aligned(byte_align))) ;
#else
#  define THEA_END_PACKED_CLASS(byte_align)  ;
#endif

namespace Thea {

// The following typedefs are here and not in NumericType.hpp since they are needed by headers that NumericType.hpp itself
// depends on.

/**
 * A signed integer suitable for indexing a structure held in memory. 32-bit on 32-bit systems, 64-bit on 64-bit systems. This
 * is expected to be fast and can usually be used for any application involving integers, unless smaller or larger storage is
 * explicitly required. It also matches the current default type of Eigen::Index. Note that it is <b>not</b> equivalent to int,
 * long or long long, whose sizes depend on LP32/64, ILP32/64 or LLP64 conventions.
 *
 * The phonetic similarity to "index" is intentional.
 *
 * @see uintx
 * @see https://rust-lang.github.io/rfcs/0544-rename-int-uint.html
 */
typedef std::ptrdiff_t intx;

/**
 * An unsigned integer suitable for indexing a structure held in memory. 32-bit on 32-bit systems, 64-bit on 64-bit systems.
 * Note that it is <b>not</b> equivalent to unsigned int, unsigned long or unsigned long long, whose sizes depend on LP32/64,
 * ILP32/64 or LLP64 conventions.
 *
 * The phonetic similarity to "(unsigned) index" is intentional.
 *
 * @see intx
 * @see https://rust-lang.github.io/rfcs/0544-rename-int-uint.html
 */
typedef std::size_t uintx;

} // namespace Thea

#endif
