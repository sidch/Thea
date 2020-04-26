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
// First version: 2009
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
                                              (intx)__LINE__)

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
