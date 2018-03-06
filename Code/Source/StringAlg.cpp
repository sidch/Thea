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

 @file stringutils.cpp

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2000-09-09
 @edited  2008-01-10
*/

#include "StringAlg.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstring>

// Not cstdarg, to make sure we pull in the vsnprintf etc functions
#include <stdarg.h>

#ifdef _WIN32
  #include <windows.h>
#endif

namespace Thea {

#ifdef _MSC_VER
// disable: "C++ exception handler used"
#  pragma warning (push)
#  pragma warning (disable : 4530)
#endif

#ifdef THEA_WINDOWS

char const * NEWLINE = "\r\n";

#else

char const * NEWLINE = "\n";

#endif

void
parseCommaSeparated(std::string const & s, TheaArray<std::string> & array, bool strip_quotes)
{
  array.clear();

  if (s == "")
  {
    return;
  }

  size_t begin = 0;
  char const delimiter = ',';
  char const quote = '\"';

  do
  {
    size_t end = begin;
    // Find the next comma, or the end of the string
    bool in_quotes = false;

    while ((end < s.length()) && (in_quotes || (s[end] != delimiter)))
    {
      if (s[end] == quote)
      {
        if ((end < s.length() - 2) && (s[end + 1] == quote) && (s[end + 2]) == quote)
        {
          // Skip over the superquote
          end += 2;
        }

        in_quotes = ! in_quotes;
      }

      ++end;
    }

    array.push_back(s.substr(begin, end - begin));
    begin = end + 1;

  } while (begin < s.length());

  if (strip_quotes)
  {
    for (size_t i = 0; i < array.size(); ++i)
    {
      std::string & t = array[i];
      size_t L = t.length();

      if ((L > 1) && (t[0] == quote) && (t[L - 1] == quote))
      {
        if ((L > 6)  && (t[1] == quote) && (t[2] == quote) && (t[L - 3] == quote) && (t[L - 2] == quote))
        {
          // Triple-quote
          t = t.substr(3, L - 6);
        }
        else
        {
          // Double-quote
          t = t.substr(1, L - 2);
        }
      }
    }
  }
}

bool
beginsWith(std::string const & test, std::string const & pattern)
{
  if (test.size() >= pattern.size())
  {
    for (size_t i = 0; i < pattern.size(); ++i)
    {
      if (pattern[i] != test[i])
      {
        return false;
      }
    }

    return true;
  }
  else
  {
    return false;
  }
}

bool
endsWith(std::string const & test, std::string const & pattern)
{
  if (test.size() >= pattern.size())
  {
    long te = (long)test.size() - 1;
    long pe = (long)pattern.size() - 1;

    for (long i = (long)pattern.size() - 1; i >= 0; --i)
    {
      if (pattern[(size_t)(pe - i)] != test[(size_t)(te - i)])
      {
        return false;
      }
    }

    return true;
  }
  else
  {
    return false;
  }
}

std::string
wordWrap(std::string const & input, long num_cols, char const * newline)
{
  std::ostringstream output;

  size_t  c = 0;
  long    len;
  // Don't make lines less than this length
  long    min_length = num_cols / 4;
  size_t  in_len = input.size();
  bool    first = true;

  while (c < in_len)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      output << newline;
    }

    if ((long)in_len - (long)c - 1 < num_cols)
    {
      // The end
      output << input.substr(c, in_len - c);
      break;
    }

    len = num_cols;

    // Look at character c + num_cols, see if it is a space.
    while ((len > min_length) &&
           (input[(size_t)(c + len)] != ' '))
    {
      len--;
    }

    if (len == min_length)
    {
      // Just crop
      len = num_cols;
    }

    output << input.substr(c, len);
    c += len;

    if (c < input.size())
    {
      // Collapse multiple spaces.
      while ((input[c] == ' ') && (c < input.size()))
      {
        c++;
      }
    }
  }

  return output.str();
}

std::string
toUpper(std::string const & s)
{
  std::string result = s;
  std::transform(result.begin(), result.end(), result.begin(), ::toupper);
  return result;
}

std::string
toLower(std::string const & s)
{
  std::string result = s;
  std::transform(result.begin(), result.end(), result.begin(), ::tolower);
  return result;
}

long
stringSplit(std::string const & s, char split_char, TheaArray<std::string> & result, bool skip_empty_fields)
{
  return stringSplit(s, std::string(1, split_char), result, skip_empty_fields);
}

long
stringSplit(std::string const & s, std::string const & split_chars, TheaArray<std::string> & result, bool skip_empty_fields)
{
  result.clear();

  // Pointer to delimiters
  char const * delims = split_chars.c_str();

  // Pointers to the beginning and end of the substring
  char const * start = s.c_str();
  char const * stop = start;

  while ((stop = std::strpbrk(start, delims)))
  {
    if (!skip_empty_fields || stop != start)
      result.push_back(std::string(start, stop - start));

    start = stop + 1;
  }

  // Append the last one
  std::string last = std::string(start);
  if (!skip_empty_fields || !last.empty())
    result.push_back(last);

  return (long)result.size();
}

std::string
trimWhitespace(std::string const & s)
{
  size_t left = 0;

  // Trim from left
  while ((left < s.length()) && isWhitespace(s[left]))
  {
    ++left;
  }

  long right = (long)s.length() - 1;

  // Trim from right
  while ((right > (long)left) && isWhitespace(s[(size_t)right]))
  {
    --right;
  }

  return s.substr(left, (size_t)(right - (long)left + 1));
}

// If your platform does not have vsnprintf, you can find a
// implementation at http://www.ijs.si/software/snprintf/

std::string __cdecl
format(char const * fmt, ...)
{
  va_list arg_list;
  va_start(arg_list, fmt);
  std::string result = vformat(fmt, arg_list);
  va_end(arg_list);

  return result;
}

#if defined(_MSC_VER) && (_MSC_VER >= 1300)
// Both MSVC seems to use the non-standard vsnprintf
// so we are using vscprintf to determine buffer size, however
// only MSVC7 and up headers include vscprintf for some reason.
std::string
vformat(char const * fmt, va_list arg_ptr)
{
  // We draw the line at a 1MB string.
  long const MAX_SIZE = 1000000;

  // If the string is less than 161 characters, allocate it on the stack because this saves the malloc/free time.
  long const BUF_SIZE = 161;

  // MSVC does not support va_copy
  long actual_size = _vscprintf(fmt, arg_ptr) + 1;

  if (actual_size > BUF_SIZE)
  {
    // Now use the heap.
    char * heap_buffer = NULL;

    if (actual_size < MAX_SIZE)
    {
      heap_buffer = (char *)std::malloc(MAX_SIZE + 1);
      _vsnprintf(heap_buffer, MAX_SIZE, fmt, arg_ptr);
      heap_buffer[MAX_SIZE] = '\0';
    }
    else
    {
      heap_buffer = (char *)std::malloc(actual_size);
      vsprintf(heap_buffer, fmt, arg_ptr);
    }

    std::string formatted_string(heap_buffer);
    std::free(heap_buffer);
    return formatted_string;
  }
  else
  {
    char stack_buffer[BUF_SIZE];
    vsprintf(stack_buffer, fmt, arg_ptr);
    return std::string(stack_buffer);
  }
}

#elif defined(_MSC_VER) && (_MSC_VER < 1300)

std::string
vformat(char const * fmt, va_list arg_ptr)
{
  // We draw the line at a 1MB string.
  long const MAX_SIZE = 1000000;

  // If the string is less than 161 characters,
  // allocate it on the stack because this saves
  // the malloc/free time.
  long const BUF_SIZE = 161;
  char stack_buffer[BUF_SIZE];

  // MSVC6 doesn't support va_copy, however it also seems to compile
  // correctly if we just pass our argument list along.  Note that
  // this whole code block is only compiled if we're on MSVC6 anyway
  long actual_written = _vsnprintf(stack_buffer, BUF_SIZE, fmt, arg_ptr);

  // Not a big enough buffer, BUF_SIZE characters written
  if (actual_written == -1)
  {
    long heap_size = 512;
    double pow_size = 1.0;
    char * heap_buffer = (char *)std::malloc(heap_size);

    while ((_vsnprintf(heap_buffer, heap_size, fmt, arg_ptr) == -1) && (heap_size < MAX_SIZE))
    {
      heap_size = std::ceil(heap_size * std::pow(2.0, pow_size++));
      heap_buffer = (char *)std::realloc(heap_buffer, heap_size);
    }

    heap_buffer[heap_size - 1] = '\0';

    std::string heap_string(heap_buffer);
    std::free(heap_buffer);

    return heap_string;
  }
  else
  {
    return std::string(stack_buffer);
  }
}

#else

// glibc 2.1 has been updated to the C99 standard
std::string
vformat(char const * fmt, va_list arg_ptr)
{
  // If the string is less than 161 characters,
  // allocate it on the stack because this saves
  // the malloc/free time.  The number 161 is chosen
  // to support two lines of text on an 80 character
  // console (plus the null terminator).
  long const BUF_SIZE = 161;
  char stack_buffer[BUF_SIZE];

  va_list arg_ptr_copy;
  va_copy(arg_ptr_copy, arg_ptr);
  long num_chars = vsnprintf(stack_buffer, BUF_SIZE, fmt, arg_ptr_copy);
  va_end(arg_ptr_copy);

  if (num_chars >= BUF_SIZE)
  {
    // We didn't allocate a big enough string.
    char * heap_buffer = (char *)std::malloc(num_chars + 1);

    debugAssertM(heap_buffer, "vformat: Heap buffer allocation failed");
    long num_chars2 = vsnprintf(heap_buffer, num_chars + 1, fmt, arg_ptr);
    debugAssertM(num_chars2 == num_chars, "vformat: Number of characters written does not match expected value");
    (void)num_chars2;  // avoid unused variable warnings

    std::string result(heap_buffer);
    std::free(heap_buffer);

    return result;
  }
  else
  {
    return std::string(stack_buffer);
  }
}

#endif

#define THEA_EOS            '\0'
#define THEA_RANGE_MATCH      1
#define THEA_RANGE_NOMATCH    0
#define THEA_RANGE_ERROR    (-1)

namespace StringAlgInternal {

static int
rangeMatch(char const * pattern, char test, int flags, char ** newp)
{
  int negate, ok;
  char c, c2;

  /*
   * A bracket expression starting with an unquoted circumflex
   * character produces unspecified results (IEEE 1003.2-1992,
   * 3.13.2).  This implementation treats it like '!', for
   * consistency with the regular expression syntax.
   * J.T. Conklin (conklin@ngai.kaleida.com)
   */
  if ((negate = (*pattern == '!' || *pattern == '^')))
    ++pattern;

  if (flags & FNM_CASEFOLD)
    test = ::tolower((unsigned char)test);

  /*
   * A right bracket shall lose its special meaning and represent
   * itself in a bracket expression if it occurs first in the list.
   * -- POSIX.2 2.8.3.2
   */
  ok = 0;
  c = *pattern++;

  do
  {
    if (c == '\\' && !(flags & FNM_NOESCAPE))
      c = *pattern++;

    if (c == THEA_EOS)
      return (THEA_RANGE_ERROR);

    if (c == '/' && (flags & FNM_PATHNAME))
      return (THEA_RANGE_NOMATCH);

    if ((flags & FNM_CASEFOLD))
      c = ::tolower((unsigned char)c);

    if (*pattern == '-'
        && (c2 = *(pattern + 1)) != THEA_EOS && c2 != ']')
    {
      pattern += 2;

      if (c2 == '\\' && !(flags & FNM_NOESCAPE))
        c2 = *pattern++;

      if (c2 == THEA_EOS)
        return (THEA_RANGE_ERROR);

      if (flags & FNM_CASEFOLD)
        c2 = ::tolower((unsigned char)c2);

      if (c <= test && test <= c2)
        ok = 1;
    }
    else if (c == test)
      ok = 1;
  }
  while ((c = *pattern++) != ']');

  *newp = (char *)pattern;
  return (ok == negate ? THEA_RANGE_NOMATCH : THEA_RANGE_MATCH);
}

int
fnmatch(char const * pattern, char const * query, int flags)
{
  char const * stringstart;
  char * newp;
  char c, test;

  for (stringstart = query;;)
  {
    switch (c = *pattern++)
    {
      case THEA_EOS:
        if ((flags & FNM_LEADING_DIR) && *query == '/')
          return (0);

        return (*query == THEA_EOS ? 0 : FNM_NOMATCH);

      case '?':
        if (*query == THEA_EOS)
          return (FNM_NOMATCH);

        if (*query == '/' && (flags & FNM_PATHNAME))
          return (FNM_NOMATCH);

        if (*query == '.' && (flags & FNM_PERIOD) &&
            (query == stringstart ||
             ((flags & FNM_PATHNAME) && *(query - 1) == '/')))
          return (FNM_NOMATCH);

        ++query;
        break;

      case '*':
        c = *pattern;

        /* Collapse multiple stars. */
        while (c == '*')
          c = *++pattern;

        if (*query == '.' && (flags & FNM_PERIOD) &&
            (query == stringstart ||
             ((flags & FNM_PATHNAME) && *(query - 1) == '/')))
          return (FNM_NOMATCH);

        /* Optimize for pattern with * at end or before /. */
        if (c == THEA_EOS)
        {
          if (flags & FNM_PATHNAME)
            return ((flags & FNM_LEADING_DIR) ||
                    std::strchr(query, '/') == NULL ?
                    0 : FNM_NOMATCH);
          else
            return (0);
        }
        else if (c == '/' && (flags & FNM_PATHNAME))
        {
          if ((query = std::strchr(query, '/')) == NULL)
            return (FNM_NOMATCH);

          break;
        }

        /* General case, use recursion. */
        while ((test = *query) != THEA_EOS)
        {
          if (fnmatch(pattern, query, flags & ~FNM_PERIOD) == 0)
            return (0);

          if (test == '/' && (flags & FNM_PATHNAME))
            break;

          ++query;
        }

        return (FNM_NOMATCH);

      case '[':
        if (*query == THEA_EOS)
          return (FNM_NOMATCH);

        if (*query == '/' && (flags & FNM_PATHNAME))
          return (FNM_NOMATCH);

        if (*query == '.' && (flags & FNM_PERIOD) &&
            (query == stringstart ||
             ((flags & FNM_PATHNAME) && *(query - 1) == '/')))
          return (FNM_NOMATCH);

        switch (rangeMatch(pattern, *query, flags, &newp))
        {
          case THEA_RANGE_ERROR:
            /* not a good range, treat as normal text */
            goto fnmatch_normal;

          case THEA_RANGE_MATCH:
            pattern = newp;
            break;

          case THEA_RANGE_NOMATCH:
            return (FNM_NOMATCH);
        }

        ++query;
        break;

      case '\\':
        if (!(flags & FNM_NOESCAPE))
        {
          if ((c = *pattern++) == THEA_EOS)
          {
            c = '\\';
            --pattern;
          }
        }

        /* FALLTHROUGH */
      default:
fnmatch_normal:
        if (c != *query && !((flags & FNM_CASEFOLD) &&
                              (::tolower((unsigned char)c) ==
                               ::tolower((unsigned char)*query))))
          return (FNM_NOMATCH);

        ++query;
        break;
    }
  }

  /* NOTREACHED */
  return -1;
}

} // namespace StringAlgInternal

bool
patternMatch(std::string const & pattern, std::string const & query, int flags)
{
  int status = StringAlgInternal::fnmatch(pattern.c_str(), query.c_str(), flags);

  switch (status)
  {
    case 0: return true;
    case FNM_NOMATCH: return false;
    default: throw Error("patternMatch: Invalid pattern");
  }
}

} // namespace Thea

#ifdef _MSC_VER
#  pragma warning (pop)
#endif
