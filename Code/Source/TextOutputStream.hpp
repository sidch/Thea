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

 @file TextOutputStream.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu
 @created 2004-06-21
 @edited  2006-10-24

 Copyright 2000-2007, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_TextOutputStream_hpp__
#define __Thea_TextOutputStream_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "NamedObject.hpp"
#include "Noncopyable.hpp"

namespace Thea {

/**
 * Convenient formatting of ASCII text written to a file. The core writeString(), writeNumber(), and writeSymbol() methods map
 * to TextInputStream's methods. writeNumber() and writeSymbol() each print an additional trailing space that is used to
 * separate adjacent tokens.
 *
 * TextOutputStream::printf allows arbitrary text to be conveniently dumped en-masse.
 *
 * When a word-wrap line break occurs, all whitespace between words is replaced with a single newline (the newline may be two
 * characters -- see TextOutputStream::Options::NewlineStyle).  Word wrapping occurs against the number of columns specified by
 * Options::numColumns, <em>minus</em> the current indent level.
 *
 * Indenting adds the specified number of spaces immediately after a newline. If a newline was followed by spaces in the
 * original string, these are added to the indent spaces.  Indenting <b>will</b> indent blank lines and will leave indents after
 * the last newline of a file (if the indent level is non-zero at the end).
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API TextOutputStream : public virtual NamedObject, private Noncopyable
{
  public:
    /**
     * Word-wrap settings (enum class).
     *
     * - NONE               Word wrapping is disabled.
     * - WITHOUT_BREAKING   Word-wrap, but don't break continuous lines that are longer than numColumns (default).
     * - ALWAYS             Wrap even if it means breaking a continuous line or a quoted string.
     *
     * Word wrapping is only allowed at whitespaces ('\\n', '\\r', '\\t', ' '); it will not occur after commas, punctuation,
     * minus signs, or any other characters
     */
    struct THEA_API WordWrap
    {
      /** Supported values. */
      enum Value
      {
        NONE,
        WITHOUT_BREAKING,
        ALWAYS
      };

      THEA_ENUM_CLASS_BODY(WordWrap)
    };

    /** Newline style (enum class). */
    struct THEA_API NewlineStyle
    {
      enum Value
      {
        WINDOWS,
        UNIX
      };

      THEA_ENUM_CLASS_BODY(NewlineStyle);
    };

    /** Output configuration optionss. */
    class THEA_API Settings
    {
      public:
        /** Word wrap mode. Defaults to WordWrap::WITHOUT_BREAKING. */
        WordWrap            wordWrap;

        /** Is word-wrapping allowed to insert newlines inside double quotes? Default: false. */
        bool                allowWordWrapInsideDoubleQuotes;

        /** Number of columns for word wrapping. Default: 8. */
        int                 numColumns;

        /** Number of spaces in each indent. Default: 4. */
        int                 spacesPerIndent;

        /**
         * Style of newline used by word wrapping and by (optionsal) conversion. Default: Windows: NewlineStyle::WINDOWS, Linux,
         * OS X: NewlineStyle::UNIX.
         */
        NewlineStyle        newlineStyle;

        /** If true, all newlines are converted to newlineStyle regardless of how they start out. Default: true. */
        bool                convertNewlines;

        /** The symbol used to represent true. Used by writeBoolean */
        std::string         trueSymbol;

        /** The symbol used to represent false. Used by writeBoolean */
        std::string         falseSymbol;

        /** Constructs the default settings. */
        Settings();

        /** Get the default settings. */
        static Settings const & defaults() { static Settings const s; return s; }

    }; // class Settings

  private:
    /**
     * Used by indentAndAppend to tell when we are writing the first character of a new line. For push/popIndent to work
     * correctly, we cannot indent immediately after writing a newline. Instead we must indent on writing the first character
     * <em>after</em> that newline.
     */
    bool                    startingNewLine;

    /** Number of characters at the end of the buffer since the last newline. */
    int                     currentColumn;

    /** True if we have seen an opening " and no closing ". */
    bool                    inDQuote;

    /** Path to the current file, or "<memory>" if writing to a memory stream. */
    std::string             path;

    /** Buffered data to write. */
    TheaArray<char>         data;

    /** Number of indents to prepend before each line.  Always set using setIndentLevel. */
    int                     indentLevel;

    /** Actual number of spaces to indent. */
    int                     indentSpaces;

    /** The newline character(s). */
    std::string             newline;

    /** Have any errors been encountered? */
    bool                    m_ok;

    /** Has the stream changed since the last commit? */
    bool                    changed;

    /** Configuration options. */
    Settings                options;

    /** Set the configuration options. */
    void setOptions(Settings const & _opt);

    /** Converts to the desired newlines. Called from vprintf(). */
    void convertNewlines(std::string const & in, std::string & out);

    /** Called from vprintf(). */
    void wordWrapIndentAppend(std::string const & str);

    /** Set the indentation level. */
    void setIndentLevel(int i);

    /** Appends the character to data, indenting whenever a newline is encountered. Called from wordWrapIndentAppend(). */
    void indentAppend(char c);

    /** Commit data to disk, optionally forcing a write even if the buffer is empty. */
    bool _commit(bool flush, bool force);

  public:
    THEA_DEF_POINTER_TYPES(TextOutputStream, shared_ptr, weak_ptr)

    /**
     * Construct a stream that writes to a file. Use "<memory>" as the path if you're going to commit to memory -- this has the
     * same effect as the default constructor. The file and its parent directories will be created if they do not exist, and the
     * file initialized to be blank. If the file cannot be constructed, ok() will return false.
     */
    explicit TextOutputStream(std::string const & path_, Settings const & settings_ = Settings::defaults());

    /** Construct a text output stream in memory that can later be committed to a string instead of a file. */
    explicit TextOutputStream(Settings const & settings_ = Settings::defaults());

    /** Destructor. Automatically calls commit(). */
    ~TextOutputStream();

    /** True if no errors have been encountered.*/
    bool ok() const;

    /** Get the path to the current file being written ("<memory>" for memory streams). */
    std::string getPath() const
    {
      return path;
    }

    /**
     * Write the buffered data to disk. Parent directories are created as needed if they do not exist.
     *
     * @param flush If true (default) the file is ready for reading when the method returns, otherwise the method returns
     *   immediately and writes the file in the background.
     *
     * @return True if the commit succeeded, else false.
     */
    bool commit(bool flush = true);

    /** Commit the buffered data to a string. */
    void commitToString(std::string & str);

    /** Commit the buffered data to a string and return it. */
    std::string commitToString();

    /**
     * Increase indent level by one.
     *
     * @see popIndent()
     */
    void pushIndent();

    /**
     * Decrease indent level by one.
     *
     * @see pushIndent()
     */
    void popIndent();

    /**
     * Write a quoted string. Special characters in the string (e.g. \\, \\t, \\n) are escaped so that TextInputStream will
     * produce the identical string on reading.
     */
    void writeString(std::string const & str);

    /** Write a boolean value, followed by a space. Uses the values of trueSymbol and falseSymbol in the settings. */
    void writeBoolean(bool b);

    /**
     * Write a double precision number, followed by a space. A 64-bit double precision number has enough bits to store a 32-bit
     * integer exactly.
     */
    void writeNumber(double n);

    /** Write a newline character(s). */
    void writeNewline();

    /** Write several newline character(s). */
    void writeNewlines(int numLines);

    /**
     * Write a symbol that can be read by TextInputStream, followed by a space. The symbol is written without quotes. Symbols
     * are required to begin with a letter or underscore and contain only letters, underscores, and numbers or be a C++ symbol
     * (e.g. "{", "(", "++" etc.) so that they may be properly parsed by TextInputStream::readSymbol().
     */
    void writeSymbol(std::string const & str);

    /**
     * Convenience function for writing multiple symbols in a row, separated by spaces, e.g.\ writeSymbols("name", "="). Empty
     * symbols are not written.
     */
    void writeSymbols(
      std::string const & a,
      std::string const & b = "",
      std::string const & c = "",
      std::string const & d = "",
      std::string const & e = "",
      std::string const & f = "");

// Indices shifted by one for member functions
#ifdef __GNUC__
#  define THEA_TEXTOUTPUTSTREAM_CHECK_PRINTF_ARGS   __attribute__((__format__(__printf__, 2, 3)))
#  define THEA_TEXTOUTPUTSTREAM_CHECK_VPRINTF_ARGS  __attribute__((__format__(__printf__, 2, 0)))
#else
#  define THEA_TEXTOUTPUTSTREAM_CHECK_PRINTF_ARGS
#  define THEA_TEXTOUTPUTSTREAM_CHECK_VPRINTF_ARGS
#endif

    /**
     * Use C-style %printf syntax to write to the stream. Follows normal %printf conventions. Note that the output will be
     * reformatted for word-wrapping and newlines.
     */
    void __cdecl printf(char const * fmt, ...) THEA_TEXTOUTPUTSTREAM_CHECK_PRINTF_ARGS;

    /**
     * Use C-style %printf syntax to write to the stream, where the format string is a std::string. Follows normal %printf
     * conventions. Note that the output will be reformatted for word-wrapping and newlines.
     *
     * @note Can't pass \a fmt by reference because that confuses va_start.
     */
    void __cdecl printf(std::string const fmt, ...);

    /**
     * Use C-style %vprintf syntax to write to the stream. Follows normal %vprintf conventions. Note that the output will be
     * reformatted for word-wrapping and newlines.
     */
    void __cdecl vprintf(char const * fmt, va_list argPtr) THEA_TEXTOUTPUTSTREAM_CHECK_VPRINTF_ARGS;

#undef THEA_TEXTOUTPUTSTREAM_CHECK_PRINTF_ARGS
#undef THEA_TEXTOUTPUTSTREAM_CHECK_VPRINTF_ARGS

}; // class TextOutputStream

} // namespace Thea

#endif
