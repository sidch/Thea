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

 @file TextInputStream.h

 Simple text lexer/tokenizer.

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @cite Based on a lexer written by Aaron Orenstein.

 @created 2002-11-27
 @edited  2010-07-03

 Copyright 2000-2010, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_TextInputStream_hpp__
#define __Thea_TextInputStream_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "NamedObject.hpp"
#include "Noncopyable.hpp"
#include "UnorderedSet.hpp"
#include <queue>
#include <stdexcept>

namespace Thea {

/** A token read from a TextInputStream. */
class THEA_API Token
{
  public:
    /** Detailed categorization of a token (enum class). */
    struct THEA_API ExtendedType
    {
      /** Supported values. */
      enum Value
      {
        DOUBLE_QUOTED,   ///< Double-quoted string.
        SINGLE_QUOTED,   ///< Single-quoted string.
        SYMBOL,          ///< Symbol (sequence of characters not fitting the other types).
        FLOATING_POINT,  ///< Floating-point number.
        INTEGER,         ///< Integer.
        BOOLEAN,         ///< Boolean value (true/false).
        LINE_COMMENT,    ///< Single-line comment.
        BLOCK_COMMENT,   ///< Multi-line (block) comment.
        NEWLINE,         ///< Newline character.
        END              ///< End-of-stream.
      };

      THEA_ENUM_CLASS_BODY(ExtendedType)
    };

    /** Basic categorization of a token. */
    struct THEA_API Type
    {
      /** Supported values. */
      enum Value
      {
        STRING   =  ExtendedType::DOUBLE_QUOTED,   ///< Quoted string.
        SYMBOL   =  ExtendedType::SYMBOL,          ///< Symbol (sequence of characters not fitting the other types).
        NUMBER   =  ExtendedType::FLOATING_POINT,  ///< Number.
        BOOLEAN  =  ExtendedType::BOOLEAN,         ///< Boolean value (true/false).
        COMMENT  =  ExtendedType::LINE_COMMENT,    ///< Single or multi-line comment.
        NEWLINE  =  ExtendedType::NEWLINE,         ///< Newline character.
        END      =  ExtendedType::END              ///< End-of-stream.
      };

      THEA_ENUM_CLASS_BODY(Type)
    };

  private:
    friend class TextInputStream;

    /** Holds the actual text value, which might be any type. If a number, it will be parsed at runtime. */
    std::string             _string;

    bool                    _bool;
    int                     _line;
    int                     _character;
    uint64                  _bytePosition;
    Type                    _type;
    ExtendedType            _extendedType;

    /** Construct a token. */
    Token(Type t, ExtendedType e, std::string const & s, int L, int c, uint64 byte)
    : _string(s), _bool(false), _line(L), _character(c), _bytePosition(byte), _type(t), _extendedType(e) {}

    Token(Type t, ExtendedType e, std::string const & s, bool b, int L, int c, uint64 byte)
    : _string(s), _bool(b), _line(L), _character(c), _bytePosition(byte), _type(t), _extendedType(e) {}

  public:
    /** Default constructor. */
    Token()
    : _bool(false),
      _line(0),
      _character(0),
      _bytePosition(0),
      _type(Type::END),
      _extendedType(ExtendedType::END)
    {}

    /** Get the type of the token. */
    Type type() const
    {
      return _type;
    }

    /** Get the extended type of the token. */
    ExtendedType extendedType() const
    {
      return _extendedType;
    }

    /**
     * The value of a single- or double-quoted string (not including the quotes), the name of a symbol, or the exact textual
     * representation of a number as parsed from the input.
     */
    std::string const & string() const
    {
      return _string;
    }

    /** Return the value of this token as a boolean. Undefined if the token is not, in fact, a boolean. */
    bool boolean() const
    {
      return _bool;
    }

    /** Starting line of the input from which this token was parsed. Starts at 1. */
    int line() const
    {
      return _line;
    }

    /** Starting character position in the input line from which this token was parsed. Starts at 1. */
    int character() const
    {
      return _character;
    }

    /** Number of bytes from the beginning of the buffer that this token was parsed from. Begins at 0 */
    uint64 bytePosition() const
    {
      return _bytePosition;
    }

    /** Return the numeric value for a number token. Undefined if the token is not, in fact, a number. */
    double number() const;

}; // class Token

/** Thrown by TextInputStream and other parsers on unexpected input. */
class THEA_API ParseError : public std::runtime_error
{
  public:
    enum { UNKNOWN = -1 };

    /** Path to the file being parsed or a short description of the source string. Empty means unknown. */
    std::string     src;

    /** For a binary file, the location of the parse error. -1 if unknown. */
    int64           byte;

    /**
     * For a text file, the line number is the line number of start of token which caused the exception. 1 is the first line of
     * the file. -1 means unknown. Note that you can use TextInputStream::Settings::startingLineNumberOffset to shift the
     * effective line number that is reported by that class.
    */
    int             line;

    /**
     * Character number (in the line) of the start of the token which caused the exception. 1 is the first character in the
     * line. May be -1 if unknown.
     */
    int             character;

    /** Default constructor. */
    ParseError() : std::runtime_error(""), byte(UNKNOWN), line(UNKNOWN), character(UNKNOWN) {}

    /** Initialize from a file path, line number, character number and error message. */
    ParseError(std::string const & src_, int l, int c, std::string const & m)
    : std::runtime_error(m), src(src_), byte(UNKNOWN), line(l), character(c) {}

    /** Initialize from a file path, byte position and error message. */
    ParseError(std::string const & src_, int64 b, std::string const & m)
    : std::runtime_error(m), src(src_), byte(b), line(UNKNOWN), character(UNKNOWN) {}

    /** Destructor. */
    ~ParseError() throw() {}

}; // class ParseError

/**
 * A simple style tokenizer for reading text files. TextInputStream handles a superset of C++, Java, Matlab, and Bash code text
 * including single line comments, block comments, quoted strings with escape sequences, and operators. TextInputStream
 * recognizes several categories of tokens, which are separated by white space, quotation marks, or the end of a
 * recognized operator:
 *
 * - Token::ExtendedType::SINGLE_QUOTED: string of characters surrounded by single quotes, e.g. 'x', '\\0', 'foo'.
 * - Token::ExtendedType::DOUBLE_QUOTED: string of characters surrounded by double quotes, e.g. "x", "abc\txyz", "b o b".
 * - Token::ExtendedType::SYMBOL: legal C++ operators, keywords, and identifiers. e.g. >=, Foo, _X, class, {
 * - Token::ExtendedType::INTEGER: numbers without decimal places or exponential notation. e.g. 10, 0x17F, 32, 0, -155
 * - Token::ExtendedType::FLOATING_POINT: numbers with decimal places or exponential notation. e.g. 1e3, -1.2, .4, 0.5
 * - Token::ExtendedType::BOOLEAN: special symbols like "true" and "false"; the exact details can be configured in
 *   TextInputStream::Settings
 * - Token::ExtendedType::LINE_COMMENT: (disabled by default); generated for line comments as specified by
 *   TextInputStream::Settings
 * - Token::ExtendedType::BLOCK_COMMENT: (disabled by default); generated for C-style block comments as specified by
 *   TextInputStream::Settings
 * - Token::ExtendedType::NEWLINE: (disabled by default); generated for any of "\\r", "\\n" or "\\r\\n"
 *
 * The special ".." and "..." tokens are always recognized in addition to normal C++ operators. Additional tokens can be made
 * available by changing the Settings.
 *
 * Negative numbers are handled specially because of the ambiguity between unary minus and negative numbers -- see the note for
 * TextInputStream::read().
 *
 * Inside quoted strings escape sequences are converted. Thus the string token for ["a\\nb"] is 'a', followed by a newline,
 * followed by 'b'. Outside of quoted strings, escape sequences are not converted, so the token sequence for [a\\nb] is symbol
 * 'a', symbol '\\', symbol 'nb' (this matches what a C++ parser would do). The exception is that a specified
 * TextInputStream::Settings::otherCommentCharacter preceeded by a backslash is assumed to be an escaped comment character and
 * is returned as a symbol token instead of being parsed as a comment (this is what a LaTeX or VRML parser would do).
 *
 * Assumes that the file is not modified once opened.
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 *
 * <b>Examples</b>
 *
 * <pre>
 * TextInputStream ti(TextInputStream::FROM_STRING, "name = 'Max', height = 6");
 *
 * Token t;
 *
 * t = ti.read();
 * assert(t.type == Token::Type::SYMBOL);
 * assert(t.sval == "name");
 *
 * t = ti.read();
 * assert(t.type == Token::Type::SYMBOL);
 * assert(t.sval == "=");
 *
 * std::string name = ti.read().sval;
 * ti.read();
 * </pre>
 *
 * <pre>
 * TextInputStream ti(TextInputStream::FROM_STRING, "name = 'Max', height = 6");
 * ti.readSymbols("name", "=");
 * std::string name = ti.readString();
 * ti.readSymbols(",", "height", "=");
 * double height = ti.readNumber();
 * </pre>
 */
class THEA_API TextInputStream : public virtual NamedObject, private Noncopyable
{
  public:
    /** Extract a number from a string. Includes MSVC specials parsing */
    static double parseNumber(std::string const & _string);

    /**
     * Extract a boolean value from a string.
     *
     * @return True if toLower(_string) == "true", else false.
     */
    static bool parseBoolean(std::string const & _string);

    /** Tokenizer configuration options. */
    struct THEA_API Settings
    {
      /**
       * If true, C-style slash-star marks a multi-line comment. See generateCommentTokens for rules on how this is applied.
       * Default is true.
       */
      bool                             cppBlockComments;

      /**
       * If true, // begins a single line comment. See generateCommentTokens for rules on how this is applied. Default is
       * true.
       */
      bool                             cppLineComments;

      /**
       * If true, otherCommentCharacter and otherCommentCharacter2 are used to begin single line comments in the same way
       * cppLineComments is. See generateCommentTokens for rules on how this is applied. Default is true.
       */
      bool                             otherLineComments;

      /**
       * If true, \\r, \\n, \\t, \\0, \\\\ and other escape sequences inside strings are converted to the equivalent C++
       * escaped character. If false, backslashes are treated literally. It is convenient to set to false if reading Windows
       * paths, for example, like c:\\foo\\bar. Default is true.
       */
      bool                             escapeSequencesInStrings;

      /**
       * If not '\\0', specifies a character that begins single line comments ('#' and '%' are popular choices). This is
       * independent of the cppLineComments flag. If the character appears in text with a backslash in front of it, it is
       * considered escaped and is not treated as a comment character. Default is '\\0'.
       */
      char                             otherCommentCharacter;

      /** Another (optional) 1-comment character. Useful for files that support multiple comment syntaxes. Default is '\\0'. */
      char                             otherCommentCharacter2;

      /**
       * If true, comments enabled by cppBlockComments, cppLineComments and otherLineComments will generate their respective
       * tokens. If false, the same settings will enable parsing and ignoring comments. Default is false.
       */
      bool                             generateCommentTokens;

      /**
       * If true, newlines will generate tokens. If false, newlines will be discarded as whitespace when parsed outside of other
       * tokens. Default is false.
       */
      bool                             generateNewlineTokens;

      /** If true, "-1" parses as the number -1 instead of the symbol "-" followed by the number 1. Default is true.*/
      bool                             signedNumbers;

      /**
       * If true, strings can be marked with single quotes (such as 'aaa'). If false, the quote character is parsed as a symbol.
       * Default is true. Backquote (`) is always parsed as a symbol.
       */
      bool                             singleQuotedStrings;

      /**
       * The character to use as a single quote. Defaults to "'" (backquote), occasionally useful to set to "`" (forward quote)
       * or to "," (comma) for reading CSV files.
       */
      char                             singleQuoteCharacter;

      /**
       * Added to the line number reported by peekLineNumber and in exceptions. Useful for concatenating files that are parsed
       * separately. Default is zero.
       */
      int                              startingLineNumberOffset;

      /**
       * Parse -1.#IND00 as the floating point number returned by Thea::Math::nan(), -1.#INF00 as -Thea::Math::inf(), and
       * 1.#INF00 as Thea::Math::inf(). Note that the C99 standard specifies that a variety of formats like "nan" are to be
       * used; these are supported by Thea::TextInputStream::Settings::simpleFloatSpecials.
       *
       * An alternative to specifying msvcFloatSpecials is to read numbers as:
       * <pre>
       *   Token x = t.read();
       *   Token y = t.peek();
       *   if ((x.string() == "-1.") &&
       *        (y.string() == "#INF00") &&
       *        (y.character() == x.character() + 3) &&
       *        (y.line() == x.line())
       *   {
       *        t.read();
       *        return nan();
       *   }
       *   // ... similar cases for inf
       * </pre>
       *
       * If the single-comment character was #, the floating point special format overrides the comment and will be parsed
       * instead.
       *
       * If signedNumbers is false msvcFloatSpecials will not be parsed.
       *
       * Default is true.
       */
      bool                             msvcFloatSpecials;

      /** Parses "+inf', "-inf", "inf", "nan" as floats instead of symbols. Defaults to true. */
      bool                             simpleFloatSpecials;

      /**
       * Parse the following set of useful proof symbols:
       *
       * =>
       * ::>
       * <::
       * :>
       * <:
       * |-
       * ::=
       * :=
       * <-
       *
       * Default is false.
       */
      bool                             proofSymbols;

      /** When parsing booleans and msvcFloatSpecials, is case significant? Default is true. */
      bool                             caseSensitive;

      /**
       * All symbols that will become the 'true' boolean token. See also caseSensitive. Clear this value to disable parsing of
       * true booleans. Default is {true}.
       */
      TheaUnorderedSet<std::string>    trueSymbols;

      /** See trueSymbols. Default is {false}. */
      TheaUnorderedSet<std::string>    falseSymbols;

      /** Constucts the default settings. */
      Settings();

      /** Get the default settings. */
      static Settings const & defaults() { static Settings const s; return s; }

    }; // struct Settings

  private:
    /** Stack of buffered tokens. */
    std::deque<Token>       stack;

    /** Characters to be tokenized. */
    TheaArray<char>         buffer;

    /** Offset of current character (the next character to consumed) in input buffer. */
    array_size_t            currentCharOffset;

    /**
     * Line number of next character to be consumed from the input buffer. (1 indicates first line of input.) Note that this is
     * the line number of the <em>next</em> character to be consumed from the input, not the line number of the <em>last</em>
     * character consumed!
     */
    int                     lineNumber;

    /**
     * Character number (within the line) of the next character to be consumed from the input buffer. (1 indicates first
     * character of the line). Note that this is the character number of the <em>next</em> character to be consumed from the
     * input, not the character number of the <em>last</em> character consumed!
     */
    int                     charNumber;

    /** Full path to file. */
    std::string             path;

    /** Full path to file or short description of string. */
    std::string             errorSourceName;

    /** Configuration options. */
    Settings                options;

    /** Initialize the parser. */
    void init();

    /**
     * Consumes the next character from the input buffer, and returns that character. Updates lineNumber and charNumber to
     * reflect the location of the next character in the input buffer.
     *
     * @note You shouldn't be using the return value of this function in most cases. In general, you should peekInputChar() to
     * get the next character, determine what to do with it, then consume it with this function (or with eatAndPeekInputChar).
     * Given that usage, in most instances you already know what this function would return!
     */
    int eatInputChar();

    /**
     * Returns the next character from the input buffer, without consuming any characters. Can also be used to look deeper into
     * the input buffer. Does not modify lineNumber or charNumber
     *
     * \param distance Index of the character in the input buffer to peek at, relative to the next character. Default is 0, for
     *   the next character in the input buffer.
     */
    int peekInputChar(int distance = 0);

    /**
     * Helper function to consume the next character in the input buffer and peek at the one following (without consuming it).
     */
    int eatAndPeekInputChar()
    {
      eatInputChar();
      return peekInputChar(0);
    }

    /** Read the next token, returning an Type::END token if no more input is available. */
    Token nextToken();

    /**
     * Helper for nextToken. Appends characters to t._string until the end delimiter is reached. When called, the next character
     * in the input buffer should be first the first character after the opening delimiter character.
     */
    void parseQuotedString(char delimiter, Token & t);

  public:
    THEA_DEF_POINTER_TYPES(TextInputStream, shared_ptr, weak_ptr)

    /** Thrown when a token cannot be read. */
    class THEA_API TokenException : public ParseError
    {
      public:
        /** Destructor. */
        ~TokenException() throw() {}

      protected:
        /** Constructor. */
        TokenException(
          std::string const & src,
          int                 ln,
          int                 ch,
          std::string const & m);
    };

    /** Thrown while parsing a number of the form 1.\#IN?00, if ? was not 'D' or 'F'. */
    class THEA_API BadMSVCSpecial : public TokenException
    {
      public:
        /** Constructor. */
        BadMSVCSpecial(
          std::string const & src,
          int                 ln,
          int                 ch);

        /** Destructor. */
        ~BadMSVCSpecial() throw() {}
    };

    /** Thrown by the read methods if a token is not of the expected type. */
    class THEA_API WrongTokenType : public TokenException
    {
      public:
        Token::Type expected;
        Token::Type actual;

        /** Constructor. */
        WrongTokenType(
          std::string const & src,
          int                 ln,
          int                 ch,
          Token::Type         e,
          Token::Type         a);

        /** Destructor. */
        ~WrongTokenType() throw() {}
    };

    /** Thrown by the read methods if a symbol string does not match the expected string. */
    class THEA_API WrongSymbol : public TokenException
    {
      public:
        std::string expected;
        std::string actual;

        /** Constructor. */
        WrongSymbol(
          std::string const & src,
          int                 ln,
          int                 ch,
          std::string const & e,
          std::string const & a);

        /** Destructor. */
        ~WrongSymbol() throw() {}
    };

    /** String read from input did not match expected string. */
    class THEA_API WrongString : public TokenException
    {
      public:
        std::string expected;
        std::string actual;

        /** Constructor. */
        WrongString(
          std::string const & src,
          int                 ln,
          int                 ch,
          std::string const & e,
          std::string const & a);

        /** Destructor. */
        ~WrongString() throw() {}
    };

    /** A flag indicting the source of a stream. */
    enum FS { FROM_STRING };

    /** Open a file for reading formatted text input. */
    explicit TextInputStream(std::string const & path_, Settings const & settings = Settings::defaults());

    /** Creates input directly from a string. The first argument must be
        TextInputStream::FROM_STRING.
    */
    TextInputStream(FS fs, std::string const & str, Settings const & settings = Settings::defaults());

    /**
     * Get the path to the file from which this input is drawn, or the first few characters of the string if created from a
     * string.
     */
    std::string getPath() const
    {
      return path;
    }

    /** Returns true while there are tokens remaining. */
    bool hasMore();

    /**
     * Read the next token (which will be the Type::END token if !hasMore()). Signed numbers can be handled in one of two
     * modes. If the option TextInputStream::Settings::signedNumbers is true, a '+' or '-' immediately before a number is
     * prepended onto that number and if there is intervening whitespace, it is read as a separate symbol. If
     * TextInputStream::Settings::signedNumbers is false, read() does not distinguish between a plus or minus symbol next to a
     * number and a positive/negative number itself. For example, "x - 1" and "x -1"  will be parsed the same way by read(). In
     * both cases, readNumber() will contract a leading "-" or "+" onto  a number.
     */
    Token read();

    /** Calls read() until the result is not a newline or comment. */
    Token readSignificant();

    /**
     * Read one token (or possibly two) as a number. If the first token in the input is a number, it is returned directly. If
     * TextInputStream::Settings::signedNumbers is false and the input stream contains a '+' or '-' symbol token immediately
     * followed by a number token, both tokens will be consumed and a single token will be returned by this method.
     *
     * WrongTokenType will be thrown if one of the input conditions described above is not satisfied. When an exception is
     * thrown, no tokens are consumed.
     */
    double readNumber();

    /** Read a boolean. If the next input token is not a boolean, throws WrongTokenType. */
    bool readBoolean();

    /**
     * Read a string token and return it. Use this method (rather than readString) if you want the token's location as well as
     * its value.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a string. When an exception is thrown, no
     * tokens are consumed.
     */
    Token readStringToken();

    /**
     * Like readStringToken, but returns the token's string. Use this method (rather than readStringToken) if you want the
     * token's value but don't really care about its location in the input. Use of readStringToken is encouraged for better
     * error reporting.
     */
    std::string readString();

    /**
     * Read a specific string token. If the next token in the input is a string matching s, it will be consumed. Use this method
     * if you want to match a specific string from the input. In that case, typically error reporting related to the token is
     * only going to occur because of a mismatch, so no location information is needed by the caller.
     *
     * WrongTokenType will be thrown if the next token in the input stream
     * is not a string. WrongString will be thrown if the next token in the
     * input stream is a string but does not match the \a s parameter. When
     * an exception is thrown, no tokens are consumed.
     *
     * @see readString(), readStringToken(), readLine()
     */
    void readString(std::string const & s);

    /**
     * Read from the beginning of the next token until the following newline and return the result as a string, ignoring all
     * parsing in between. The newline is not returned in the string, and the following token read will be a newline or end of
     * file token (if they are enabled for parsing).
     */
    std::string readLine();

    /**
     * Read a comment token and return it. Use this method (rather than readComment) if you want the token's location as well as
     * its value.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a comment. When an exception is thrown, no
     * tokens are consumed.
     */
    Token readCommentToken();

    /**
     * Like readCommentToken(), but returns the token's string.
     *
     * Use this method (rather than readCommentToken) if you want the token's value but don't really care about its location in
     * the input. Use of readCommentToken is encouraged for better error reporting.
    */
    std::string readComment();

    /**
     * Read a specific comment token. If the next token in the input is a comment matching s, it will be consumed. Use this
     * method if you want to match a specific comment from the input. In that case, typically error reporting related to the
     * token is only going to occur because of a mismatch, so no location information is needed by the caller.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a comment. WrongString will be thrown if the
     * next token in the input stream is a comment but does not match the s parameter. When an exception is thrown, no tokens
     * are consumed.
     */
    void readComment(std::string const & s);

    /**
     * Read a newline token and return it. Use this method (rather than readNewline) if you want the token's location as well as
     * its value. WrongTokenType will be thrown if the next token in the input stream is not a newline. When an exception is
     * thrown, no tokens are consumed.
     */
    Token readNewlineToken();

    /**
     * Like readNewlineToken(), but returns the token's string. Use this method (rather than readNewlineToken) if you want the
     * token's value but don't really care about its location in the input. Use of readNewlineToken() is encouraged for better
     * error reporting.
     */
    std::string readNewline();

    /**
     * Read a specific newline token. If the next token in the input is a newline matching \a s, it will be consumed. Use this
     * method if you want to match a specific newline from the input. In that case, typically error reporting related to the
     * token is only going to occur because of a mismatch, so no location information is needed by the caller.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a newline. WrongString will be thrown if the
     * next token in the input stream is a newline but does not match the \a s parameter. When an exception is thrown, no
     * tokens are consumed.
     */
    void readNewline(std::string const & s);

    /**
     * Read a symbol token and return it. Use this method (rather than readSymbol) if you want the token's location as well as
     * its value.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a symbol. When an exception is thrown, no
     * tokens are consumed.
     */
    Token readSymbolToken();

    /**
     * Like readSymbolToken(), but returns the token's string. Use this method (rather than readSymbolToken) if you want the
     * token's value but don't really care about its location in the input. Use of readSymbolToken() is encouraged for better
     * error reporting.
     */
    std::string readSymbol();

    /**
     * Read a specific symbol token. If the next token in the input is a %symbol matching \a symbol, it will be consumed. Use
     * this method if you want to match a specific %symbol from the input. In that case, typically error reporting related to
     * the token is only going to occur because of a mismatch, so no location information is needed by the caller.
     *
     * WrongTokenType will be thrown if the next token in the input stream is not a %symbol. WrongSymbol will be thrown if the
     * next token in the input stream is a %symbol but does not match the \a symbol parameter. When an exception is thrown, no
     * tokens are consumed.
     */
    void readSymbol(std::string const & symbol);

    /** Read a series of two specific symbols. See readSymbol(). */
    void readSymbols(std::string const & s1, std::string const & s2)
    {
      readSymbol(s1);
      readSymbol(s2);
    }

    /** Read a series of three specific symbols. See readSymbol(). */
    void readSymbols(
      std::string const & s1,
      std::string const & s2,
      std::string const & s3)
    {
      readSymbol(s1);
      readSymbol(s2);
      readSymbol(s3);
    }

    /** Read a series of four specific symbols. See readSymbol(). */
    void readSymbols(
      std::string const & s1,
      std::string const & s2,
      std::string const & s3,
      std::string const & s4)
    {
      readSymbol(s1);
      readSymbol(s2);
      readSymbol(s3);
      readSymbol(s4);
    }

    /** Get a copy of the next token in the input stream, but don't remove it from the input stream. */
    Token peek();

    /**
     * Get the line number for the <em>next</em> token.
     *
     * @see peek(), peekCharacterNumber()
     */
    int peekLineNumber();

    /**
     * Get the character number (relative to the line) for the <em>next</em> token in the input stream.
     *
     * @see peek(), peekLineNumber()
     */
    int peekCharacterNumber();

    /**
     * Take a previously read token and push it back at the front of the input stream.
     *
     * Can be used in the case where more than one token of read-ahead is needed (i.e. when peek doesn't suffice).
     */
    void push(Token const & t);

}; // class TextInputStream

} // namespace Thea

#endif
