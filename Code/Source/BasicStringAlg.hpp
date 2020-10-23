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
// First version: 2013
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
inline intx
findFirstSlash(std::string const & f, intx start = 0)
{
  size_t pos = f.find_first_of("/\\", (size_t)start);
  return (pos == std::string::npos ? -1 : (intx)pos);
}

/**
 * Finds the index of the last '\\' or '/' character, starting at index \a start (starts at the end of the string if \a start is
 * negative).
 *
 * @return The index of the last slash if one is found, else a negative number.
 *
 * @see findFirstSlash, isSlash
 */
inline intx
findLastSlash(std::string const & f, intx start = -1)
{
  if (start == -1)
    start = (intx)f.length() - 1;

  size_t pos = f.find_last_of("/\\", (size_t)start);
  return (pos == std::string::npos ? -1 : (intx)pos);
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
THEA_API std::string wordWrap(std::string const & input, intx num_cols, char const * newline = NEWLINE);

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
#  define THEA_CHECK_PRINTF_ARGS          __attribute__((__format__(__printf__, 1, 2)))
#  define THEA_CHECK_VPRINTF_ARGS         __attribute__((__format__(__printf__, 1, 0)))
#  define THEA_CHECK_MEMBER_PRINTF_ARGS   __attribute__((__format__(__printf__, 2, 3)))
#  define THEA_CHECK_MEMBER_VPRINTF_ARGS  __attribute__((__format__(__printf__, 2, 0)))
#else
#  define THEA_CHECK_PRINTF_ARGS
#  define THEA_CHECK_VPRINTF_ARGS
#  define THEA_CHECK_MEMBER_PRINTF_ARGS
#  define THEA_CHECK_MEMBER_VPRINTF_ARGS
#endif

/**
 * Produces a string from arguments in the style of printf. This avoids problems with buffer overflows when using sprintf and
 * makes it easy to use the result functionally.  This function is fast when the resulting string is under 160 characters (not
 * including terminator) and slower when the string is longer.
 */
THEA_API std::string THEA_CDECL format(char const * fmt, ...) THEA_CHECK_PRINTF_ARGS;

/**
 * Produces a string from arguments in the style of printf, can be called with the argument list from a function that itself
 * takes variable arguments. This avoids problems with buffer overflows when using sprintf and makes it easy to use the result
 * functionally.  This function is fast when the resulting string is under 160 characters (not including terminator) and slower
 * when the string is longer.
 */
THEA_API std::string vformat(char const * fmt, va_list arg_ptr) THEA_CHECK_VPRINTF_ARGS;

} // namespace Thea

#endif
