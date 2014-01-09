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

/**
 * \file
 * This file <b>must</b> be included in <em>every</em> file in the project, <em>before</em> any other include.
 */

#ifndef __Thea_Platform_hpp__
#define __Thea_Platform_hpp__

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

#ifndef _MSC_VER
#  define __cdecl
#endif

/** Use tight alignment for a class. Useful for small classes that must be packed into an array. */
#ifdef _MSC_VER
#  define THEA_BEGIN_PACKED_CLASS(byte_align)  PRAGMA( pack(push, byte_align) )
#else
#  define THEA_BEGIN_PACKED_CLASS(byte_align)
#endif

/**
 * @def PRAGMA(expression)
 * \#pragma may not appear inside a macro, so this uses the pragma operator to create an equivalent statement.
 */
#ifdef _MSC_VER
// Microsoft's version http://msdn.microsoft.com/en-us/library/d9x1s805.aspx
#  define PRAGMA(x) __pragma(x)
#else
// C99 standard http://www.delorie.com/gnu/docs/gcc/cpp_45.html
#  define PRAGMA(x) _Pragma(#x)
#endif

/** Mark end of class that uses tight alignment. */
#ifdef _MSC_VER
#  define THEA_END_PACKED_CLASS(byte_align)  ; PRAGMA( pack(pop) )
#elif defined(__GNUC__)
#  define THEA_END_PACKED_CLASS(byte_align)  __attribute((aligned(byte_align))) ;
#else
#  define THEA_END_PACKED_CLASS(byte_align)  ;
#endif

#endif
