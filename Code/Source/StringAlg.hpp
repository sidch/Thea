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

#ifndef __Thea_StringAlg_hpp_
#define __Thea_StringAlg_hpp_

#include "Common.hpp"
#include "BasicStringAlg.hpp"
#include "Array.hpp"
#include <sstream>

#ifdef THEA_WINDOWS
#  ifndef FNM_NOMATCH
#    define FNM_NOMATCH         1      /* Match failed. */
#    define FNM_NOESCAPE     0x01      /* Disable backslash escaping. */
#    define FNM_PATHNAME     0x02      /* Slash must be matched by slash. */
#    define FNM_PERIOD       0x04      /* Period must be matched by period. */
#    define FNM_LEADING_DIR  0x08      /* Ignore /<tail> after Imatch. */
#    define FNM_CASEFOLD     0x10      /* Case insensitive search. */
#  endif
#else  // On non-windows systems, include fnmatch directly. Both GNU/Linux and OS X support all the flags above.
#  include <fnmatch.h>
#endif

namespace Thea {

/**
 * Separates a comma-separated line, properly escaping commas within double quotes (") and super quotes ("""). This matches
 * Microsoft Excel's CSV output.
 *
 * @param s The string to parse.
 * @param array The comma-separated fields are stored here.
 * @param strip_quotes If true, strips leading and trailing " and """.
 *
 * @see stringSplit, TextInput
 */
THEA_API void parseCommaSeparated(std::string const & s, TheaArray<std::string> & array, bool strip_quotes = true);

/**
 * Split a string at each occurrence of a splitting character.
 *
 * @param s The string to split.
 * @param split_char The delimiting character.
 * @param result Used to return the sequence of fields found.
 * @param skip_empty_fields If true, a sequence of delimiters is treated as a single delimiter.
 *
 * @return The number of fields found.
 */
THEA_API long stringSplit(std::string const & s, char split_char, TheaArray<std::string> & result,
                          bool skip_empty_fields = false);

/**
 * Split a string at each occurrence of any splitting character from a provided set.
 *
 * @param s The string to split.
 * @param split_chars The set of delimiting characters. E.g. to split a string on whitespace, use
 *   <tt>split_chars = " \t\n\f\r"</tt>.
 * @param result Used to return the sequence of fields found.
 * @param skip_empty_fields If true, a sequence of delimiters is treated as a single delimiter.
 *
 * @return The number of fields found.
 */
THEA_API long stringSplit(std::string const & s, std::string const & split_chars, TheaArray<std::string> & result,
                          bool skip_empty_fields = false);

/**
 * Concatenate a sequence of serializable objects (typically strings), separated by a padding character, string or other
 * serializable object.
 */
template <typename IteratorT, typename C>
std::string
stringJoin(IteratorT begin, IteratorT end, C const & padding)
{
  std::ostringstream out;
  for (IteratorT iter = begin; iter != end; ++iter)
  {
    if (iter != begin) out << padding;
    out << *iter;
  }

  return out.str();
}

/**
 * Concatenate a sequence of serializable objects (typically strings), separated by a padding character, string or other
 * serializable object. This is shorthand for stringJoin(a.begin(), a.end(), padding).
 */
template <typename ContainerT, typename C>
std::string
stringJoin(ContainerT const & a, C const & padding)
{
  return stringJoin(a.begin(), a.end(), padding);
}

/** Pattern matching flags (enum class). Correspond to flags for POSIX fnmatch. */
struct THEA_API Match
{
  /** Supported values. */
  enum Value
  {
    /** Treat backslash as an ordinary character, instead of an escape character. */
    NOESCAPE     =  FNM_NOESCAPE,

    /**
     * Match a slash in the query only with a slash in the pattern and not by an asterisk (*) or a question mark (?)
     * metacharacter, nor by a bracket expression ([]) containing a slash.
     */
    PATHNAME     =  FNM_PATHNAME,

    /**
     * A leading period in the query has to be matched exactly by a period in the pattern. A period is considered to be leading
     * if it is the first character in \a query, or if both FNM_PATHNAME is set and the period immediately follows a slash.
     */
    PERIOD       =  FNM_PERIOD,

    /**
     * If this flag (a GNU extension) is set, the pattern is considered to be matched if it matches an initial segment of the
     * query which is followed by a slash. This flag is mainly for the internal use of glibc and is only implemented in certain
     * cases.
     */
    LEADING_DIR  =  FNM_LEADING_DIR,

    /** If this flag (a GNU extension) is set, the pattern is matched case-insensitively. */
    CASEFOLD     =  FNM_CASEFOLD
  };

  THEA_ENUM_CLASS_BODY(Match)

}; // struct Match

/**
 * Compares a filename or pathname to a pattern. Equivalent (except for boolean instead of integer return value) to function
 * fnmatch() as specified in POSIX 1003.2-1992, section B.6.
 *
 * The function checks whether the \a query argument matches the \a pattern argument, which is a shell wildcard pattern. The
 * \a flags argument modifies the behaviour; it is the bitwise OR of zero or more Match flags.
 *
 * @return True if \a query matches \a pattern, false otherwise. Throws an error if the pattern is malformed.
 */
THEA_API bool patternMatch(std::string const & pattern, std::string const & query, int flags = 0);

} // namespace Thea

#endif
