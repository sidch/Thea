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

 @file TextOutputStream.cpp

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu
 @created 2004-06-21
 @edited  2010-03-14

 Copyright 2000-2010, Morgan McGuire.
 All rights reserved.
*/

#include "TextOutputStream.hpp"
#include "FilePath.hpp"
#include "FileSystem.hpp"
#include "Math.hpp"
#include <cstdarg>
#include <stdio.h>

namespace Thea {

TextOutputStream::Settings::Settings()
: wordWrap(WordWrap::WITHOUT_BREAKING),
  allowWordWrapInsideDoubleQuotes(false),
  numColumns(80),
  spacesPerIndent(4),
  convertNewlines(true),
  trueSymbol("true"),
  falseSymbol("false")
{
#ifdef THEA_WINDOWS
  newlineStyle = NewlineStyle::WINDOWS;
#else
  newlineStyle = NewlineStyle::UNIX;
#endif
}

TextOutputStream::TextOutputStream(TextOutputStream::Settings const & opt)
: NamedObject("<memory>"),
  startingNewLine(true),
  currentColumn(0),
  inDQuote(false),
  path("<memory>"),
  indentLevel(0),
  m_ok(true),
  changed(false)
{
  setOptions(opt);
}

TextOutputStream::TextOutputStream(std::string const & path_, TextOutputStream::Settings const & opt)
: NamedObject(FilePath::objectName(path_)),
  startingNewLine(true),
  currentColumn(0),
  inDQuote(false),
  path(FileSystem::resolve(path_)),
  indentLevel(0),
  m_ok(true),
  changed(false)
{
  setOptions(opt);

  // Verify ability to write to disk
  _commit(false, true);
}

TextOutputStream::~TextOutputStream()
{
  if (path != "<memory>")
    commit(true);
}

bool
TextOutputStream::ok() const
{
  return m_ok;
}

bool
TextOutputStream::commit(bool flush)
{
  return _commit(flush, false);
}

bool
TextOutputStream::_commit(bool flush, bool force)
{
  // If there is already an error, the commit fails
  if (!m_ok)
    return false;

  // Is there anything new to write?
  if (!force && !changed)
    return true;

  // Nothing to commit for memory streams
  if (path == "<memory>")
    return true;

  // Make sure the directory exists
  std::string dir = FilePath::parent(path);
  if (!FileSystem::exists(dir))
    if (!FileSystem::createDirectory(dir))
    {
      THEA_ERROR << "TextOutputStream: Could not create parent directory of '" << path << "'";
      m_ok = false;
    }

  FILE * f = NULL;
  if (m_ok)
  {
    f = fopen(path.c_str(), "wb");
    if (!f)
    {
      THEA_ERROR << "TextOutputStream: Could not open \'" + path + "\' for writing";
      m_ok = false;
    }
  }

  if (m_ok)
  {
    fwrite(&data[0], 1, data.size(), f);

    if (flush)
      fflush(f);

    fclose(f);

    changed = false;
  }

  return m_ok;
}

void
TextOutputStream::commitToString(std::string & out)
{
  out = std::string(&data[0], data.size());
}

std::string
TextOutputStream::commitToString()
{
  std::string str;
  commitToString(str);
  return str;
}

void
TextOutputStream::setIndentLevel(int i)
{
  indentLevel = i;

  // If there were more pops than pushes, don't let that take us below 0 indent.
  // Don't ever indent more than the number of columns.
  indentSpaces = Math::clamp(options.spacesPerIndent * indentLevel,
                             0,
                             options.numColumns - 1);
}

void
TextOutputStream::setOptions(Settings const & _opt)
{
  if (options.numColumns < 2)
    throw Error(getNameStr() + ": Must specify at least 2 columns in options");

  options = _opt;

  setIndentLevel(indentLevel);
  newline = (options.newlineStyle == NewlineStyle::WINDOWS ? "\r\n" : "\n");
}

void
TextOutputStream::pushIndent()
{
  setIndentLevel(indentLevel + 1);
}

void
TextOutputStream::popIndent()
{
  setIndentLevel(indentLevel - 1);
}

namespace TextOutputStreamInternal {

static std::string
escape(std::string const & str)
{
  std::string result = "";

  for (size_t i = 0; i < str.length(); ++i)
  {
    char c = str.at(i);

    switch (c)
    {
      case '\0':
        result += "\\0";
        break;

      case '\r':
        result += "\\r";
        break;

      case '\n':
        result += "\\n";
        break;

      case '\t':
        result += "\\t";
        break;

      case '\\':
        result += "\\\\";
        break;

      default:
        result += c;
    }
  }

  return result;
}

} // namespace TextOutputStreamInternal

void
TextOutputStream::writeString(std::string const & str)
{
  // Convert special characters to escape sequences
  this->printf("\"%s\"", TextOutputStreamInternal::escape(str).c_str());
}

void
TextOutputStream::writeBoolean(bool b)
{
  this->printf("%s ", b ? options.trueSymbol.c_str() : options.falseSymbol.c_str());
}

void
TextOutputStream::writeNumber(double n)
{
  this->printf("%f ", n);
}

void
TextOutputStream::writeSymbol(std::string const & str)
{
  if (!str.empty())
  {
    // TODO: check for legal symbols?
    this->printf("%s ", str.c_str());
  }
}

void TextOutputStream::writeSymbols(
  std::string const & a,
  std::string const & b,
  std::string const & c,
  std::string const & d,
  std::string const & e,
  std::string const & f)
{
  writeSymbol(a);
  writeSymbol(b);
  writeSymbol(c);
  writeSymbol(d);
  writeSymbol(e);
  writeSymbol(f);
}

void
TextOutputStream::printf(std::string const formatString, ...)
{
  va_list argList;
  va_start(argList, formatString);
  this->vprintf(formatString.c_str(), argList);
  va_end(argList);
}

void
TextOutputStream::printf(char const * formatString, ...)
{
  va_list argList;
  va_start(argList, formatString);
  this->vprintf(formatString, argList);
  va_end(argList);
}

void
TextOutputStream::convertNewlines(std::string const & in, std::string & out)
{
  // TODO: can be significantly optimized in cases where
  // single characters are copied in order by walking through
  // the array and copying substrings as needed.
  if (options.convertNewlines)
  {
    out = "";

    for (size_t i = 0; i < in.size(); ++i)
    {
      if (in[i] == '\n')
      {
        // Unix newline
        out += newline;
      }
      else if ((in[i] == '\r') && (i + 1 < in.size()) && (in[i + 1] == '\n'))
      {
        // Windows newline
        out += newline;
        ++i;
      }
      else
      {
        out += in[i];
      }
    }
  }
  else
  {
    out = in;
  }
}

void
TextOutputStream::writeNewline()
{
  for (size_t i = 0; i < newline.size(); ++i)
  {
    indentAppend(newline[i]);
  }
}

void
TextOutputStream::writeNewlines(int numLines)
{
  for (int i = 0; i < numLines; ++i)
  {
    writeNewline();
  }
}

void
TextOutputStream::wordWrapIndentAppend(std::string const & str)
{
  // TODO: keep track of the last space character we saw so we don't
  // have to always search.
  if ((options.wordWrap == WordWrap::NONE) ||
      (currentColumn + (int)str.size() <= options.numColumns))
  {
    // No word-wrapping is needed

    // Add one character at a time.
    // TODO: optimize for strings without newlines to add multiple
    // characters.
    for (size_t i = 0; i < str.size(); ++i)
    {
      indentAppend(str[i]);
    }

    return;
  }

  // Number of columns to wrap against
  int cols = options.numColumns - indentSpaces;

  // Copy forward until we exceed the column size,
  // and then back up and try to insert newlines as needed.
  for (size_t i = 0; i < str.size(); ++i)
  {
    indentAppend(str[i]);

    if ((str[i] == '\r') && (i + 1 < str.size()) && (str[i + 1] == '\n'))
    {
      // \r\n, we need to hit the \n to enter word wrapping.
      ++i;
      indentAppend(str[i]);
    }

    if (currentColumn >= cols)
    {
      debugAssertM(str[i] != '\n' && str[i] != '\r',
                   "Should never enter word-wrapping on a newline character");
      // True when we're allowed to treat a space as a space.
      bool unquotedSpace = options.allowWordWrapInsideDoubleQuotes || ! inDQuote;
      // Cases:
      //
      // 1. Currently in a series of spaces that ends with a newline
      //     strip all spaces and let the newline
      //     flow through.
      //
      // 2. Currently in a series of spaces that does not end with a newline
      //     strip all spaces and replace them with single newline
      //
      // 3. Not in a series of spaces
      //     search backwards for a space, then execute case 2.
      // Index of most recent space
      size_t lastSpace = data.size() - 1;
      // How far back we had to look for a space
      size_t k = 0;
      size_t maxLookBackward = (size_t)(currentColumn - indentSpaces);

      // Search backwards (from current character), looking for a space.
      while ((k < maxLookBackward) &&
             (lastSpace > 0) &&
             (! ((data[lastSpace] == ' ') && unquotedSpace)))
      {
        --lastSpace;
        ++k;

        if ((data[lastSpace] == '\"') && !options.allowWordWrapInsideDoubleQuotes)
        {
          unquotedSpace = ! unquotedSpace;
        }
      }

      if (k == maxLookBackward)
      {
        // We couldn't find a series of spaces
        if (options.wordWrap == WordWrap::ALWAYS)
        {
          // Strip the last character we wrote, force a newline,
          // and replace the last character;
          data.pop_back();
          writeNewline();
          indentAppend(str[i]);
        }
        else
        {
          // Must be Settings::WRAP_WITHOUT_BREAKING
          //
          // Don't write the newline; we'll come back to
          // the word wrap code after writing another character
        }
      }
      else
      {
        // We found a series of spaces.  If they continue
        // to the new string, strip spaces off both.  Otherwise
        // strip spaces from data only and insert a newline.
        // Find the start of the spaces.  firstSpace is the index of the
        // first non-space, looking backwards from lastSpace.
        size_t firstSpace = lastSpace;

        while ((k < maxLookBackward) &&
               (firstSpace > 0) &&
               (data[firstSpace] == ' '))
        {
          --firstSpace;
          ++k;
        }

        if (k == maxLookBackward)
        {
          ++firstSpace;
        }

        if (lastSpace == data.size() - 1)
        {
          // Spaces continued up to the new string
          data.resize(firstSpace + 1);
          writeNewline();

          // Delete the spaces from the new string
          while ((i < str.size() - 1) && (str[i + 1] == ' '))
          {
            ++i;
          }
        }
        else
        {
          // Spaces were somewhere in the middle of the old string.
          // replace them with a newline.
          // Copy over the characters that should be saved
          TheaArray<char> temp;

          for (size_t j = lastSpace + 1; j < data.size(); ++j)
          {
            char c = data[j];

            if (c == '\"')
            {
              // Undo changes to quoting (they will be re-done
              // when we paste these characters back on).
              inDQuote = !inDQuote;
            }

            temp.push_back(c);
          }

          // Remove those characters and replace with a newline.
          data.resize(firstSpace + 1);
          writeNewline();

          // Write them back
          for (size_t j = 0; j < temp.size(); ++j)
          {
            indentAppend(temp[j]);
          }

          // We are now free to continue adding from the
          // new string, which may or may not begin with spaces.
        } // if spaces included new string
      } // if hit indent
    } // if line exceeded
  } // iterate over str

  changed = true;
}

void
TextOutputStream::indentAppend(char c)
{
  if (startingNewLine)
  {
    for (int j = 0; j < indentSpaces; ++j)
    {
      data.push_back(' ');
    }

    startingNewLine = false;
    currentColumn = indentSpaces;
  }

  data.push_back(c);

  // Don't increment the column count on return character
  // newline is taken care of below.
  if (c != '\r')
  {
    ++currentColumn;
  }

  if (c == '\"')
  {
    inDQuote = ! inDQuote;
  }

  startingNewLine = (c == '\n');

  if (startingNewLine)
  {
    currentColumn = 0;
  }

  changed = true;
}

void
TextOutputStream::vprintf(char const * formatString, va_list argPtr)
{
  std::string str = vformat(formatString, argPtr);
  std::string clean;
  convertNewlines(str, clean);
  wordWrapIndentAppend(clean);
}

} // namespace Thea
