//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

/*
 ORIGINAL HEADER

 @file stringutils.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @author  2000-09-09
 @edited  2010-03-05
*/

#ifndef __Thea_BasicStringAlg_hpp_
#define __Thea_BasicStringAlg_hpp_

#include "Platform.hpp"
#include <cctype>
#include <cstdarg>
#include <string>

namespace Thea {

/** The newline character sequence for the current platform. LF for Unix (Linux, Mac OS X), CR-LF for Windows. */
extern THEA_API char const * NEWLINE;

/**
 * Finds the index of the first '\\' or '/' character, starting at index \a start.
 *
 * @return The index of the first slash if one is found, else a negative number.
 *
 * @see findLastSlash, isSlash
 */
inline long
findFirstSlash(std::string const & f, long start = 0)
{
  size_t pos = f.find_first_of("/\\", (size_t)start);
  return (pos == std::string::npos ? -1 : (long)pos);
}

/**
 * Finds the index of the last '\\' or '/' character, starting at index \a start (starts at the end of the string if \a start is
 * negative).
 *
 * @return The index of the last slash if one is found, else a negative number.
 *
 * @see findFirstSlash, isSlash
 */
inline long
findLastSlash(std::string const & f, long start = -1)
{
  if (start == -1)
    start = (long)f.length() - 1;

  size_t pos = f.find_last_of("/\\", (size_t)start);
  return (pos == std::string::npos ? -1 : (long)pos);
}

/** Check if the test string begins with the pattern string. */
THEA_API bool beginsWith(std::string const & test, std::string const & pattern);

/** Check if the test string ends with the pattern string. */
THEA_API bool endsWith(std::string const & test, std::string const & pattern);

/**
 * Produces a new string that is the input string wrapped at a certain number of columns (where the line is broken at the latest
 * space before the column limit). Platform specific newlines are inserted to wrap, or a specific "newline" character may be
 * specified.
 */
THEA_API std::string wordWrap(std::string const & input, long num_cols, char const * newline = NEWLINE);

/** Get the uppercase version of a string. */
THEA_API std::string toUpper(std::string const & s);

/** Get the lowercase version of a string. */
THEA_API std::string toLower(std::string const & s);

/** Strips whitespace from both ends of the string. */
THEA_API std::string trimWhitespace(std::string const & s);

/** Check if a character is a whitespace character. */
inline bool
isWhitespace(char c)
{
  return std::isspace(c) != 0;
}

/** Check if a character is a newline character. */
inline bool
isNewline(char c)
{
  return (c == '\n') || (c == '\r');
}

/** Check if a character is a digit. */
inline bool
isDigit(char c)
{
  return std::isdigit(c) != 0;
}

/** Check if a character is a letter of the alphabet. */
inline bool
isAlpha(char c)
{
  return std::isalpha(c) != 0;
}

/** Check if a character is a slash (forward or backward). */
inline bool
isSlash(char c)
{
  return (c == '\\') || (c == '/');
}

/** Check if a character is a quote (single or double). */
inline bool
isQuote(char c)
{
  return (c == '\'') || (c == '\"');
}

#ifdef __GNUC__
#  define THEA_CHECK_PRINTF_ARGS   __attribute__((__format__(__printf__, 1, 2)))
#  define THEA_CHECK_VPRINTF_ARGS  __attribute__((__format__(__printf__, 1, 0)))
#else
#  define THEA_CHECK_PRINTF_ARGS
#  define THEA_CHECK_VPRINTF_ARGS
#endif

/**
 * Produces a string from arguments in the style of printf. This avoids problems with buffer overflows when using sprintf and
 * makes it easy to use the result functionally.  This function is fast when the resulting string is under 160 characters (not
 * including terminator) and slower when the string is longer.
 */
THEA_API std::string __cdecl format(char const * fmt, ...) THEA_CHECK_PRINTF_ARGS;

/**
 * Produces a string from arguments in the style of printf, can be called with the argument list from a function that itself
 * takes variable arguments. This avoids problems with buffer overflows when using sprintf and makes it easy to use the result
 * functionally.  This function is fast when the resulting string is under 160 characters (not including terminator) and slower
 * when the string is longer.
 */
THEA_API std::string vformat(char const * fmt, va_list arg_ptr) THEA_CHECK_VPRINTF_ARGS;

#undef THEA_CHECK_PRINTF_ARGS
#undef THEA_CHECK_VPRINTF_ARGS

} // namespace Thea

#endif
