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

 @file TextInputStream.cpp

 @author Morgan McGuire, graphics3d.com

 @cite Based on a lexer written by Aaron Orenstein.

 @created 2001-11-27
 @edited  2010-07-03
*/

#include "TextInputStream.hpp"
#include "FilePath.hpp"
#include "FileSystem.hpp"
#include "Math.hpp"
#include <cstdio>
#include <cstring>

#ifdef _MSC_VER
#   pragma warning (push)
// conversion from 'int' to 'char', possible loss of data (TODO: fix underlying problems)
#   pragma warning (disable: 4244)
#endif

namespace Thea {

Token
TextInputStream::readSignificant()
{
  Token t;

  do
  {
    t = read();
  }
  while ((t.type() == Token::Type::COMMENT) || (t.type() == Token::Type::NEWLINE));

  return t;
}

double
Token::number() const
{
  if (_type == Type::NUMBER)
  {
    return TextInputStream::parseNumber(_string);
  }
  else
  {
    return 0.0;
  }
}

bool
TextInputStream::parseBoolean(std::string const & _string)
{
  return toLower(_string) == "true";
}

double
TextInputStream::parseNumber(std::string const & _string)
{
  std::string s = toLower(_string);

  if (s == "-1.#ind00" || s == "nan")
  {
    return Math::nan<double>();
  }

  if (s == "1.#inf00" || s == "inf" || s == "+inf")
  {
    return Math::inf<double>();
  }

  if (s == "-1.#inf00" || s == "-inf")
  {
    return -Math::inf<double>();
  }

  double n;

  if ((_string.length() > 2) &&
      (_string[0] == '0') &&
      (_string[1] == 'x'))
  {
    // Hex
    uint32 i;
    std::sscanf(_string.c_str(), "%x", &i);
    n = i;
  }
  else
  {
    std::sscanf(_string.c_str(), "%lg", &n);
  }

  return n;
}

TextInputStream::Settings::Settings()
: cppBlockComments(true),
  cppLineComments(true),
  otherLineComments(true),
  escapeSequencesInStrings(true),
  otherCommentCharacter('\0'),
  otherCommentCharacter2('\0'),
  generateCommentTokens(false),
  generateNewlineTokens(false),
  signedNumbers(true),
  singleQuotedStrings(true),
  singleQuoteCharacter('\''),
  startingLineNumberOffset(0),
  msvcFloatSpecials(true),
  simpleFloatSpecials(true),
  proofSymbols(false),
  caseSensitive(true)
{
  trueSymbols.insert("true");
  falseSymbols.insert("false");
}

Token
TextInputStream::peek()
{
  if (stack.size() == 0)
  {
    Token t = nextToken();
    push(t);
  }

  return stack.front();
}

int
TextInputStream::peekLineNumber()
{
  return peek().line();
}

int
TextInputStream::peekCharacterNumber()
{
  return peek().character();
}

Token
TextInputStream::read()
{
  if (stack.size() > 0)
  {
    Token t = stack.front();
    stack.pop_front();
    return t;
  }
  else
  {
    return nextToken();
  }
}

std::string
TextInputStream::readLine()
{
  // Go to the front of the next token
  Token t = read();

  // Reset the position to the start of this token
  currentCharOffset = (size_t)t.bytePosition();
  stack.clear();

  if (currentCharOffset == buffer.size())
  {
    // End of file
    return "";
  }

  std::string s;

  // Read until newline or eof
  char c = '\0';
  do
  {
    c = buffer[currentCharOffset];

    if (c == '\r' || c == '\n')
    {
      // Done
      break;
    }
    else
    {
      s += c;
      ++currentCharOffset;
    }
  }
  while (currentCharOffset < buffer.size());

  return s;
}

static void
toUpper(TheaUnorderedSet<std::string> & set)
{
  TheaArray<std::string> symbols(set.begin(), set.end());
  set.clear();

  for (size_t i = 0; i < symbols.size(); ++i)
  {
    set.insert(toUpper(symbols[i]));
  }
}

void
TextInputStream::init()
{
  currentCharOffset = 0;
  charNumber = 1;
  lineNumber = 1 + options.startingLineNumberOffset;

  if (! options.caseSensitive)
  {
    // Convert true and false symbols to all uppercase for fast comparisons
    toUpper(options.trueSymbols);
    toUpper(options.falseSymbols);
  }
}

void
TextInputStream::push(Token const & t)
{
  stack.push_front(t);
}

bool
TextInputStream::hasMore()
{
  return (peek()._type != Token::Type::END);
}

int
TextInputStream::eatInputChar()
{
  // Don't go off the end
  if (currentCharOffset >= buffer.size())
  {
    return EOF;
  }

  char c = buffer[currentCharOffset];
  ++currentCharOffset;

  // update lineNumber and charNumber to reflect the location of the *next*
  // character which will be read.

  // increment line number for \r, \n and \r\n which matches Token::Type::NEWLINE parsing
  if (c == '\r')
  {
    ++lineNumber;
    charNumber = 1;

    // check for \r\n
    if (currentCharOffset < buffer.size())
    {
      char c2 = buffer[currentCharOffset];

      if (c2 == '\n')
      {
        c = c2;
        ++currentCharOffset;
      }
    }
  }
  else if (c == '\n')
  {
    ++lineNumber;
    charNumber = 1;
  }
  else
  {
    ++charNumber;
  }

  return c;
}

int
TextInputStream::peekInputChar(int distance)
{
  // Don't go off the end
  if ((size_t)(currentCharOffset + distance) >= buffer.size())
  {
    return EOF;
  }

  char c = buffer[currentCharOffset + distance];
  return (int)c;
}

Token
TextInputStream::nextToken()
{
  Token t;
  t._bytePosition = currentCharOffset;
  t._line         = lineNumber;
  t._character    = charNumber;
  t._type         = Token::Type::END;
  t._extendedType = Token::ExtendedType::END;
  int c = peekInputChar();

  if (c == EOF)
  {
    return t;
  }

  // loop through white space, newlines and comments
  // found before other tokens
  bool whitespaceDone = false;

  while (! whitespaceDone)
  {
    whitespaceDone = true;

    // generate newlines tokens for '\n' and '\r' and '\r\n'
    while (isWhitespace(c))
    {
      if (options.generateNewlineTokens && isNewline(c))
      {
        t._type         = Token::Type::NEWLINE;
        t._extendedType = Token::ExtendedType::NEWLINE;
        t._bytePosition = currentCharOffset;
        t._line         = lineNumber;
        t._character    = charNumber;
        t._string       = c;
        int c2 = peekInputChar(1);

        if (c == '\r' && c2 == '\n')
        {
          t._string  += c2;
        }

        eatInputChar();
        return t;
      }
      else
      {
        // Consume the single whitespace
        c = eatAndPeekInputChar();
      }
    }

    // update line and character number to include discarded whitespace
    t._line         = lineNumber;
    t._character    = charNumber;
    t._bytePosition = (uint64)currentCharOffset;
    int c2 = peekInputChar(1);

    // parse comments and generate tokens if enabled
    std::string commentString;

    // check for line comments first
    bool isLineComment = false;

    if (options.cppLineComments && (c == '/' && c2 == '/'))
    {
      // set start of line comment and eat markers
      isLineComment = true;
      eatInputChar();
      eatInputChar();
    }
    else if ( options.otherCommentCharacter &&
              (options.otherCommentCharacter != '\0' && c == options.otherCommentCharacter) )
    {
      // set start of line comment and eat markers
      isLineComment = true;
      eatInputChar();
    }
    else if ( options.otherCommentCharacter &&
              (options.otherCommentCharacter2 != '\0' && c == options.otherCommentCharacter2) )
    {
      // set start of line comment and eat markers
      isLineComment = true;
      eatInputChar();
    }

    if (isLineComment)
    {
      // consume line comment to newline or EOF
      c = peekInputChar();

      while (! isNewline(c) && c != EOF)
      {
        // build comment string for token
        commentString += c;
        c = eatAndPeekInputChar();
      }

      if (options.generateCommentTokens)
      {
        t._type         = Token::Type::COMMENT;
        t._extendedType = Token::ExtendedType::LINE_COMMENT;
        t._string       = commentString;
        return t;
      }
      else
      {
        // There is whitespace after the comment (in particular, the
        // newline that terminates the comment).  There might also be
        // whitespace at the start of the next line.
        whitespaceDone = false;
      }
    }
    else if (options.cppBlockComments && (c == '/' && c2 == '*'))
    {
      // consume block comment to end-marker or EOF
      // consume both start-comment chars, can't let the trailing one
      // help close the comment.
      eatInputChar();
      eatInputChar();
      // c is the next character we'll read, c2 is the one after *that*
      c = peekInputChar();
      c2 = peekInputChar(1);

      while (! ((c == '*') && (c2 == '/')) && (c != EOF))
      {
        commentString += c;
        // Eat input char may consume more than one character if there is a newline
        eatInputChar();
        c = peekInputChar();
        c2 = peekInputChar(1);
      }

      eatInputChar();      // eat closing '*'
      eatInputChar();      // eat closing '/'
      c = peekInputChar();

      if (options.generateCommentTokens)
      {
        t._type         = Token::Type::COMMENT;
        t._extendedType = Token::ExtendedType::BLOCK_COMMENT;
        t._string       = commentString;
        return t;
      }
      else
      {
        // There is whitespace after the comment (in particular, the
        // newline that terminates the comment).  There might also be
        // whitespace at the start of the next line.
        whitespaceDone = false;
      }
    }
  } // while (! whitespaceDone)

  t._line      = lineNumber;
  t._character = charNumber;
  t._bytePosition = (uint64)currentCharOffset;

  // handle EOF
  if (c == EOF)
  {
    return t;
  }

  // Extended ASCII parses as itself, except for EOF
  if (c > 127 && c < 255)
  {
    t._type = Token::Type::SYMBOL;
    t._extendedType = Token::ExtendedType::SYMBOL;
    t._string = c;
    c = eatAndPeekInputChar();
  }

  // Perform appropriate setup for a symbol (including setting up the token
  // string to start with c), eat the input character, and overwrite
  // 'c' with the peeked next input character.
#define THEA_SETUP_SYMBOL(c)                                                \
  {                                                                         \
    t._type = Token::Type::SYMBOL;                                          \
    t._extendedType = Token::ExtendedType::SYMBOL;                          \
    t._string = c;                                                          \
    c = eatAndPeekInputChar();                                              \
  }

  switch (c)
  {
    case '@':                   // Simple symbols -> just themselves.
    case '(':
    case ')':
    case ',':
    case ';':
    case '{':
    case '}':
    case '[':
    case ']':
    case '#':
    case '$':
    case '?':
    case '%':
      THEA_SETUP_SYMBOL(c);
      return t;

    case '-':                   // negative number, -, --, -=, or ->
      THEA_SETUP_SYMBOL(c);

      switch (c)
      {
        case '>':               // ->
        case '-':               // --
        case '=':               // -=
          t._string += c;
          eatInputChar();
          return t;
      }

      if (options.signedNumbers)
      {
        if (isDigit(c) || (c == '.' && isDigit(peekInputChar(1))))
        {
          // Negative number.  'c' is still the first digit, and is
          // the next input char.
          goto numLabel;
        }
        else
        {
          char terminal = peekInputChar(3);

          if (options.simpleFloatSpecials && (c == 'i') && (peekInputChar(1) == 'n') && (peekInputChar(2) == 'f') &&
              ! isAlpha(terminal) && (terminal != '_'))
          {
            // negative infinity
            t._type = Token::Type::NUMBER;
            t._extendedType = Token::ExtendedType::FLOATING_POINT;
            t._string = "-inf";
            eatInputChar(); // i
            eatInputChar(); // n
            eatInputChar(); // f
            return t;
          }
        }
      }

      // plain -
      return t;

    case '+':                   // positive number, +, ++, or +=
      THEA_SETUP_SYMBOL(c);

      switch (c)
      {
        case '+':               // ++
        case '=':               // +=
          t._string += c;
          eatInputChar();
          return t;
      }

      if (options.signedNumbers)
      {
        if (isDigit(c) || (c == '.' && isDigit(peekInputChar(1))))
        {
          // Positive number.  'c' is still the first digit, and is
          // the next input char.
          goto numLabel;
        }
        else
        {
          char terminal = peekInputChar(3);

          if (options.simpleFloatSpecials && (c == 'i') && (peekInputChar(1) == 'n') && (peekInputChar(2) == 'f') &&
              ! isAlpha(terminal) && (terminal != '_'))
          {
            // positive infinity
            t._type = Token::Type::NUMBER;
            t._extendedType = Token::ExtendedType::FLOATING_POINT;
            t._string = "+inf";
            eatInputChar(); // i
            eatInputChar(); // n
            eatInputChar(); // f
            return t;
          }
        }
      }

      return t;

    case ':':                   // : or :: or ::> or ::= or := or :>
      THEA_SETUP_SYMBOL(c);

      if (c == ':')
      {
        t._string += c;
        eatInputChar();

        if (options.proofSymbols)
        {
          c = peekInputChar(0);

          if ((c == '>') || (c == '='))
          {
            t._string += c;
            eatInputChar();
          }
        }
      }
      else if (options.proofSymbols && (c == '=' || c == '>'))
      {
        t._string += c;
        eatInputChar();
      }

      return t;

    case '=':                   // = or == or =>
      THEA_SETUP_SYMBOL(c);

      if (c == '=')
      {
        t._string += c;
        eatInputChar();
        return t;
      }
      else if (options.proofSymbols && (c == '>'))
      {
        t._string += c;
        eatInputChar();
        return t;
      }

      return t;

    case '*':                   // * or *=
    case '/':                   // / or /=
    case '!':                   // ! or !=
    case '~':                   // ~ or ~=
    case '^':                   // ^ or ^=
      THEA_SETUP_SYMBOL(c);

      if (c == '=')
      {
        t._string += c;
        eatInputChar();
        return t;
      }

      return t;

    case '>':                   // >, >>,or >=
    case '<':                   // <<, <<, or <= or <- or <:
    case '|':                   // ||, ||, or |= or |-
    case '&':                   // &, &&, or &=
    {
      int orig_c = c;
      THEA_SETUP_SYMBOL(c);

      if ((c == '=') || (orig_c == c))
      {
        t._string += c;
        eatInputChar();
        return t;
      }
      else if (options.proofSymbols)
      {
        if ((orig_c == '<') && (c == '-'))
        {
          t._string += c;
          eatInputChar();
        }
        else if ((orig_c == '|') && (c == '-'))
        {
          t._string += c;
          eatInputChar();
        }
        else if ((orig_c == '<') && (c == ':'))
        {
          t._string += c;
          c = eatAndPeekInputChar();

          if (c == ':')
          {
            t._string += c;
            eatInputChar();
          }
        }
      }
    }

    return t;

    case '\\':                // backslash or escaped comment char.
      THEA_SETUP_SYMBOL(c);

      if ((options.otherCommentCharacter != '\0'
           && c == options.otherCommentCharacter)
          || (options.otherCommentCharacter2 != '\0'
              && c == options.otherCommentCharacter2))
      {
        // Escaped comment character. Return the raw comment char (no backslash).
        t._string = c;
        eatInputChar();
        return t;
      }

      return t;

    case '.':                   // number, ., .., or ...
      if (isDigit(peekInputChar(1)))
      {
        // We're parsing a float that began without a leading zero
        goto numLabel;
      }

      THEA_SETUP_SYMBOL(c);

      if (c == '.')           // .. or ...
      {
        t._string += c;
        c = eatAndPeekInputChar();

        if (c == '.')       // ...
        {
          t._string += c;
          eatInputChar();
        }

        return t;
      }

      return t;
  } // switch (c)

#undef THEA_SETUP_SYMBOL
numLabel:

  if (isDigit(c) || (c == '.'))
  {
    // A number.  Note-- single dots have been
    // parsed already, so a . indicates a number
    // less than 1 in floating point form.

    // [0-9]*(\.[0-9][f]) or [0-9]+ or 0x[0-9,A-F]+
    if (t._string != "-")
    {
      // If we picked up a leading "-" sign above, keep it,
      // otherwise drop the string parsed thus far
      t._string = "";
    }

    t._type = Token::Type::NUMBER;

    if (c == '.')
    {
      t._extendedType = Token::ExtendedType::FLOATING_POINT;
    }
    else
    {
      t._extendedType = Token::ExtendedType::INTEGER;
    }

    if ((c == '0') && (peekInputChar(1) == 'x'))
    {
      // Hex number
      t._string += "0x";
      // skip the 0x
      eatInputChar();
      eatInputChar();
      c = peekInputChar();

      while (isDigit(c) || ((c >= 'A') && (c <= 'F')) || ((c >= 'a') && (c <= 'f')))
      {
        t._string += c;
        c = eatAndPeekInputChar();
      }
    }
    else
    {
      // Non-hex number

      // Read the part before the decimal.
      while (isDigit(c))
      {
        t._string += c;
        c = eatAndPeekInputChar();
      }

      // True if we are reading a floating-point special type
      bool isSpecial = false;

      // Read the decimal, if one exists
      if (c == '.')
      {
        t._extendedType = Token::ExtendedType::FLOATING_POINT;
        // The '.' character was a decimal point, not the start of a
        // method or range operator
        t._string += c;
        c = eatAndPeekInputChar();

        // Floating point specials (msvc format only)
        if (options.msvcFloatSpecials && (c == '#'))
        {
          isSpecial = true;
          // We are reading a floating point special value
          // of the form -1.#IND00, -1.#INF00, or 1.#INF00
          c = eatAndPeekInputChar();
          char test = c;

          if (! options.caseSensitive)
          {
            test = toupper(c);
          }

          if (test != 'I')
          {
            throw BadMSVCSpecial
            (
              "Incorrect floating-point special (inf or nan) format",
              t.line(), charNumber);
          }

          c = eatAndPeekInputChar();
          test = c;

          if (! options.caseSensitive)
          {
            test = toupper(c);
          }

          if (test != 'N')
          {
            throw BadMSVCSpecial
            (
              "Incorrect floating-point special (inf or nan) format",
              t.line(), charNumber);
          }

          t._string += "#IN";
          c = eatAndPeekInputChar();
          test = c;

          if (! options.caseSensitive)
          {
            test = toupper(c);
          }

          if ((test != 'F') && (test != 'D'))
          {
            throw BadMSVCSpecial
            (
              "Incorrect floating-point special (inf or nan) format",
              t.line(), charNumber);
          }

          t._string += c;

          for (int j = 0; j < 2; ++j)
          {
            c = eatAndPeekInputChar();

            if (c != '0')
            {
              throw BadMSVCSpecial
              (
                "Incorrect floating-point special (inf or nan) format",
                t.line(), charNumber);
            }

            t._string += (char)c;
          }
        }
        else
        {
          // Read the part after the decimal
          while (isDigit((char)c))
          {
            t._string += (char)c;
            c = eatAndPeekInputChar();
          }
        }
      }

      if (! isSpecial && ((c == 'e') || (c == 'E')))
      {
        // Read exponent
        t._extendedType = Token::ExtendedType::FLOATING_POINT;
        t._string += c;
        c = eatAndPeekInputChar();

        if ((c == '-') || (c == '+'))
        {
          t._string += c;
          c = eatAndPeekInputChar();
        }

        while (isDigit(c))
        {
          t._string += c;
          c = eatAndPeekInputChar();
        }
      }

      if (! isSpecial && (t._extendedType == Token::ExtendedType::FLOATING_POINT) && (c == 'f'))
      {
        // Trailing f on a float
        t._string += c;
        c = eatAndPeekInputChar();
      }
    }

    return t;
  }
  else if (isAlpha(c) || (c == '_'))
  {
    // Identifier or keyword
    // [A-Za-z_][A-Za-z_0-9]*
    t._type = Token::Type::SYMBOL;
    t._extendedType = Token::ExtendedType::SYMBOL;
    t._string = "";

    do
    {
      t._string += c;
      c = eatAndPeekInputChar();
    }
    while (isAlpha(c) || isDigit(c) || (c == '_'));

    // See if this symbol is actually a boolean
    if ((options.trueSymbols.size() > 0) || (options.falseSymbols.size() > 0))
    {
      std::string str = t._string;

      if (! options.caseSensitive)
      {
        str = toUpper(str);
      }

      if (options.trueSymbols.find(str) != options.trueSymbols.end())
      {
        t._type = Token::Type::BOOLEAN;
        t._extendedType = Token::ExtendedType::BOOLEAN;
        t._bool = true;
      }
      else if (options.falseSymbols.find(str) != options.falseSymbols.end())
      {
        t._type = Token::Type::BOOLEAN;
        t._extendedType = Token::ExtendedType::BOOLEAN;
        t._bool = false;
      }
    }

    if (options.simpleFloatSpecials && ((t._string == "nan") || (t._string == "inf")))
    {
      t._type = Token::Type::NUMBER;
      t._extendedType = Token::ExtendedType::FLOATING_POINT;
    }

    return t;
  }
  else if (c == '\"')
  {
    // Discard the double-quote.
    eatInputChar();
    // Double quoted string
    parseQuotedString('\"', t);
    return t;
  }
  else if (c == options.singleQuoteCharacter)
  {
    // Discard the single-quote.
    eatInputChar();

    if (options.singleQuotedStrings)
    {
      // Single quoted string
      parseQuotedString(options.singleQuoteCharacter, t);
    }
    else
    {
      t._string = c;
      t._type = Token::Type::SYMBOL;
      t._extendedType = Token::ExtendedType::SYMBOL;
    }

    return t;
  } // end of special case tokens

  if (c == EOF)
  {
    t._type = Token::Type::END;
    t._extendedType = Token::ExtendedType::END;
    t._string = "";
    return t;
  }

  // Some unknown token
  debugAssertM(false,
               format("Unrecognized token type beginning with character '%c' (ASCII %d)",
                      c, (int)c));
  return t;
}

void
TextInputStream::parseQuotedString(char delimiter, Token & t)
{
  t._type = Token::Type::STRING;

  if (delimiter == options.singleQuoteCharacter)
  {
    t._extendedType = Token::ExtendedType::SINGLE_QUOTED;
  }
  else
  {
    t._extendedType = Token::ExtendedType::DOUBLE_QUOTED;
  }

  while (true)
  {
    // We're definitely going to consume the next input char, so we get
    // it right now.  This makes the condition handling below a bit easier.
    int c = eatInputChar();

    if (c == EOF)
    {
      // Type::END inside a quoted string.  (We finish the string.)
      break;
    }

    if (options.escapeSequencesInStrings && (c == '\\'))
    {
      // An escaped character.  We're definitely going to consume it,
      // so we get it (and consume it) now.
      c = eatInputChar();

      switch (c)
      {
        case 'r':
          t._string += '\r';
          break;

        case 'n':
          t._string += '\n';
          break;

        case 't':
          t._string += '\t';
          break;

        case '0':
          t._string += '\0';
          break;

        case '\\':
        case '\"':
          t._string += (char)c;
          break;

        default:
          if (c ==  options.singleQuoteCharacter)
          {
            t._string += (char)c;
            break;
          }

          if (((c == options.otherCommentCharacter) &&
               (options.otherCommentCharacter != '\0')) ||
              ((c == options.otherCommentCharacter2) &&
               (options.otherCommentCharacter2 != '\0')))
          {
            t._string += c;
          }

          // otherwise, some illegal escape sequence; skip it.
          break;
      } // switch
    }
    else if (c == delimiter)
    {
      // End of the string.  Already consumed the character.
      break;
    }
    else
    {
      // All other chars, go on to the string.  Already consumed the
      // character.
      t._string += (char)c;
    }
  }
}

bool
TextInputStream::readBoolean()
{
  Token t(read());

  if (t._type == Token::Type::BOOLEAN)
  {
    return t.boolean();
  }

  // Push initial token back, and throw an error.  We intentionally
  // indicate that the wrong type is the type of the initial token.
  // Logically, the number started there.
  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::BOOLEAN, t._type);
}

double
TextInputStream::readNumber()
{
  Token t(read());

  if (t._type == Token::Type::NUMBER)
  {
    return t.number();
  }

  // Even if signedNumbers is disabled, readNumber attempts to
  // read a signed number, so we handle that case here.
  if (! options.signedNumbers
      && (t._type == Token::Type::SYMBOL)
      && ((t._string == "-")
          || (t._string == "+")))
  {
    Token t2(read());

    if ((t2._type == Token::Type::NUMBER)
        && (t2._character == t._character + 1))
    {
      if (t._string == "-")
      {
        return -t2.number();
      }
      else
      {
        return t2.number();
      }
    }

    // push back the second token.
    push(t2);
  }

  // Push initial token back, and throw an error.  We intentionally
  // indicate that the wrong type is the type of the initial token.
  // Logically, the number started there.
  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::NUMBER, t._type);
}

Token
TextInputStream::readStringToken()
{
  Token t(read());

  if (t._type == Token::Type::STRING)  // fast path
  {
    return t;
  }

  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::STRING, t._type);
}

std::string
TextInputStream::readString()
{
  return readStringToken()._string;
}

void
TextInputStream::readString(std::string const & s)
{
  Token t(readStringToken());

  if (t._string == s)  // fast path
  {
    return;
  }

  push(t);
  throw WrongString(getName(), t.line(), t.character(),
                    s, t._string);
}

Token
TextInputStream::readCommentToken()
{
  Token t(read());

  if (t._type == Token::Type::COMMENT)  // fast path
  {
    return t;
  }

  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::COMMENT, t._type);
}

std::string
TextInputStream::readComment()
{
  return readCommentToken()._string;
}

void
TextInputStream::readComment(std::string const & s)
{
  Token t(readCommentToken());

  if (t._string == s)  // fast path
  {
    return;
  }

  push(t);
  throw WrongString(getName(), t.line(), t.character(),
                    s, t._string);
}

Token
TextInputStream::readNewlineToken()
{
  Token t(read());

  if (t._type == Token::Type::NEWLINE)  // fast path
  {
    return t;
  }

  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::NEWLINE, t._type);
}

std::string
TextInputStream::readNewline()
{
  return readNewlineToken()._string;
}

void
TextInputStream::readNewline(std::string const & s)
{
  Token t(readNewlineToken());

  if (t._string == s)  // fast path
  {
    return;
  }

  push(t);
  throw WrongString(getName(), t.line(), t.character(),
                    s, t._string);
}

Token
TextInputStream::readSymbolToken()
{
  Token t(read());

  if (t._type == Token::Type::SYMBOL)  // fast path
  {
    return t;
  }

  push(t);
  throw WrongTokenType(getName(), t.line(), t.character(),
                       Token::Type::SYMBOL, t._type);
}

std::string
TextInputStream::readSymbol()
{
  return readSymbolToken()._string;
}

void
TextInputStream::readSymbol(std::string const & symbol)
{
  Token t(readSymbolToken());

  if (t._string == symbol)  // fast path
  {
    return;
  }

  push(t);
  throw WrongSymbol(getName(), t.line(), t.character(),
                    symbol, t._string);
}

TextInputStream::TextInputStream(std::string const & path_, Settings const & opt)
: NamedObject(FilePath::objectName(path_)),
  path(FileSystem::resolve(path_)),
  options(opt)
{
  errorSourceName = getName();

  init();
  std::string input = FileSystem::readWholeFile(path);

  size_t n = input.size();
  buffer.resize((size_t)n);
  std::memcpy(&buffer[0], input.c_str(), n);
}

TextInputStream::TextInputStream(FS fs, std::string const & str, Settings const & opt)
: options(opt)
{
  (void)fs;
  init();

  if (str.length() < 14)
  {
    setName(std::string("\"") + str + "\"");
  }
  else
  {
    setName(std::string("\"") + str.substr(0, 10) + "...\"");
  }

  errorSourceName = getName();

  buffer.resize((size_t)str.length());  // we don't bother copying trailing NUL
  std::memcpy(&buffer[0], str.c_str(), str.length());
}

///////////////////////////////////////////////////////////////////////////////////

TextInputStream::TokenException::TokenException(
  std::string const & src,
  int                 ln,
  int                 ch,
  std::string const & m)
: ParseError(src, ln, ch, format("%s(%d): %s", src.c_str(), ln, m.c_str()))
{
}

///////////////////////////////////////////////////////////////////////////////////

static char const *
tokenTypeToString(Token::Type t)
{
  switch (t)
  {
    case Token::Type::SYMBOL:
      return "SYMBOL";

    case Token::Type::STRING:
      return "STRING";

    case Token::Type::NUMBER:
      return "NUMBER";

    case Token::Type::END:
      return "END";

    case Token::Type::NEWLINE:
      return "NEWLINE";

    default:
      debugAssertM(false, "Fell through switch");
      return "?";
  }
}

TextInputStream::WrongTokenType::WrongTokenType(
  std::string const & src,
  int                 ln,
  int                 ch,
  Token::Type         e,
  Token::Type         a)
: TokenException(src, ln, ch, format("Expected token of type %s, found type %s", tokenTypeToString(e), tokenTypeToString(a))),
  expected(e),
  actual(a)
{
}

TextInputStream::BadMSVCSpecial::BadMSVCSpecial(
  std::string const & src,
  int                 ln,
  int                 ch)
: TokenException(src, ln, ch, "Bad MSVC special character")
{
}

TextInputStream::WrongSymbol::WrongSymbol(
  std::string const & src,
  int                 ln,
  int                 ch,
  std::string const & e,
  std::string const & a)
: TokenException(src, ln, ch, format("Expected symbol '%s', found symbol '%s'", e.c_str(), a.c_str())),
  expected(e),
  actual(a)
{
}

TextInputStream::WrongString::WrongString(
  std::string const & src,
  int                 ln,
  int                 ch,
  std::string const & e,
  std::string const & a)
: TokenException(src, ln, ch, format("Expected string '%s', found string '%s'", e.c_str(), a.c_str())),
  expected(e),
  actual(a)
{
}

} // namespace Thea

#ifdef _MSC_VER
#   pragma warning (pop)
#endif
