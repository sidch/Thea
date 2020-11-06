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
#include <cstdarg>

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
class THEA_API TextOutputStream : public NamedObject, private Noncopyable
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
    Array<char>             data;

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
    THEA_DECL_SMART_POINTERS(TextOutputStream)

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

    /**
     * Use C-style %printf syntax to write to the stream. Follows normal %printf conventions. Note that the output will be
     * reformatted for word-wrapping and newlines.
     */
    void THEA_CDECL printf(char const * fmt, ...) THEA_CHECK_MEMBER_PRINTF_ARGS;

    /**
     * Use C-style %vprintf syntax to write to the stream. Follows normal %vprintf conventions. Note that the output will be
     * reformatted for word-wrapping and newlines.
     */
    void THEA_CDECL vprintf(char const * fmt, va_list arg_list) THEA_CHECK_MEMBER_VPRINTF_ARGS;

}; // class TextOutputStream

} // namespace Thea

#endif
