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

#ifndef __Thea_Log_hpp__
#define __Thea_Log_hpp__

#include "Platform.hpp"
#include "BasicStringAlg.hpp"
#include "Noncopyable.hpp"
#include "Spinlock.hpp"
#include <iostream>
#include <string>

namespace Thea {

namespace LogInternal {

extern Spinlock lock;

// Extract the filename from a full path.
THEA_API std::string stripPathFromFilename(std::string const & full_path);

// Get the current date and time as a string (not threadsafe).
THEA_API std::string currentDateTimeToString();

} // namespace LogInternal

/**
 * A temporary object that locks an output stream on construction and unlocks it (after optionally writing a newline) on
 * destruction. All objects piped to the object in a single line are written atomically. Useful e.g. for writing log messages to
 * the same file/console from multiple threads.
 *
 * Example:
 * \code
 *   LockedOutputStream<>(std::cout).getStream() << "This " << " line " << " will " << " be " << " written " << " atomically";
 * \endcode
 *
 * @note Currently, all objects of this class share the same lock, regardless of the wrapped stream.
 */
template < typename StreamT = std::basic_ostream<char> >
class /* THEA_API */ LockedOutputStream : public Noncopyable
{
  public:
    typedef StreamT StreamType;  ///< The type of wrapped output stream.

    /** Constructor. */
    LockedOutputStream(StreamT & stream_, bool append_newline_ = false)
    : stream(&stream_), append_newline(append_newline_)
    {
      setOutputLock(true);
    }

    /** Constructor. Writes a prefix to the stream after acquiring the output lock. */
    LockedOutputStream(StreamT & stream_, std::string const & prefix, bool append_newline_ = false)
    : stream(&stream_), append_newline(append_newline_)
    {
      setOutputLock(true);
      getStream() << prefix;
    }

    /** Destructor. Writes a newline to the output stream if the object was constructed with <tt>append_newline = true</tt>. */
    ~LockedOutputStream()
    {
      if (append_newline)
        getStream() << std::endl;

      setOutputLock(false);
    }

    /** Get the locked stream. */
    StreamT & getStream() { return *stream; }

  private:
    StreamT * stream;  ///< The wrapped output stream.
    bool append_newline;  ///< If true, a newline is written to the output stream when this object is destroyed.

    /** Lock/unlock the common lock */
    void setOutputLock(bool value)
    {
      if (value)
        LogInternal::lock.lock();
      else
        LogInternal::lock.unlock();
    }

}; // class LockedOutputStream

} // namespace Thea

// Fully qualify references in #defines so they can be used in client programs in non-Thea namespaces without namespace errors.

#define THEA_LOG_STANDARD_PREFIX Thea::format("[%s] %s:%ld: ", \
                                              Thea::LogInternal::currentDateTimeToString().c_str(), \
                                              Thea::LogInternal::stripPathFromFilename(__FILE__).c_str(), \
                                              (long)__LINE__)

/**
 * Synchronized console output stream, with no line prefix. Outputs a newline at the end of every sequence of stream operations
 * (i.e. after every stack such as <code>THEA_CONSOLE << a << b << c;</code>).
 */
#define THEA_CONSOLE Thea::LockedOutputStream<>(std::cout, true).getStream()

/**
 * Synchronized logging stream, with a prefix indicating the time, source file and line number. Outputs a newline at the end of
 * every sequence of stream operations (i.e. after every stack such as <code>THEA_LOG << a << b << c;</code>).
 */
#define THEA_LOG Thea::LockedOutputStream<>(std::cout, THEA_LOG_STANDARD_PREFIX, true).getStream()

#ifdef THEA_DEBUG_BUILD
/**
 * Synchronized stream for debug messages, with a prefix indicating the time, source file and line number. Outputs a newline at
 * the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_DEBUG << a << b << c;</code>).
 *
 * Deactivated in release mode.
 */
#  define THEA_DEBUG Thea::LockedOutputStream<>(std::cout, THEA_LOG_STANDARD_PREFIX, true).getStream()
#else
/**
 * Synchronized stream for debug messages, with a prefix indicating the time, source file and line number. Outputs a newline at
 * the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_DEBUG << a << b << c;</code>).
 *
 * Deactivated in release mode.
 */
#  define THEA_DEBUG while (false) std::cout
#endif

/**
 * Synchronized stream for error messages, with a prefix indicating the time, source file and line number. Outputs a newline at
 * the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_ERROR << a << b << c;</code>).
 */
#define THEA_ERROR Thea::LockedOutputStream<>(std::cerr, THEA_LOG_STANDARD_PREFIX + "ERROR: ", true).getStream()

/**
 * Synchronized stream for warning messages, with a prefix indicating the time, source file and line number. Outputs a newline
 * at the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_WARNING << a << b << c;</code>).
 */
#define THEA_WARNING Thea::LockedOutputStream<>(std::cerr, THEA_LOG_STANDARD_PREFIX + "WARNING: ", true).getStream()

#endif
