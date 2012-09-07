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

#ifndef Thea_Log_hpp
#define Thea_Log_hpp

#include "G3D/G3D.h"
#include "Platform.hpp"
#include <iostream>
#include <string>

namespace Thea {
namespace LogInternal {

// Extract the filename from a full path.
THEA_API std::string stripPathFromFilename(std::string const & fullPath);

// Get the current date and time as a string (not threadsafe).
THEA_API std::string currentDateTimeToString();

// A temporary object that locks an output stream on construction and unlocks it (after writing a newline) on destruction. All
// objects of this class share the same lock, regardless of the wrapped stream.
class THEA_API LockedOutputStream
{
  public:
    // Constructor.
    LockedOutputStream(std::ostream & stream_);

    // Constructor. Writes a prefix string to the stream after acquiring the output lock.
    LockedOutputStream(std::ostream & stream_, std::string const & prefix);

    // Destructor.
    ~LockedOutputStream();

    // Get the locked stream.
    std::ostream & getStream() { return *stream; }

  private:
    std::ostream * stream;

    // Lock/unlock the common lock
    void setOutputLock(bool value);

}; // class LockedOutputStream

} // namespace LogInternal
} // namespace Thea

// Fully qualify references in #defines so they can be used in client programs in non-Thea namespaces without namespace errors.

#define THEA_LOG_STANDARD_PREFIX G3D::format("[%s] %s:%ld: ", \
                                             Thea::LogInternal::currentDateTimeToString().c_str(), \
                                             Thea::LogInternal::stripPathFromFilename(__FILE__).c_str(), \
                                             (long)__LINE__)

/**
 * Synchronized console output stream, with no line prefix. Outputs a newline at the end of every sequence of stream operations
 * (i.e. after every stack such as <code>THEA_CONSOLE << a << b << c;</code>).
 */
#define THEA_CONSOLE Thea::LogInternal::LockedOutputStream(std::cout).getStream()

/**
 * Synchronized logging stream, with a prefix indicating the time, source file and line number. Outputs a newline at the end of
 * every sequence of stream operations (i.e. after every stack such as <code>THEA_LOG << a << b << c;</code>).
 */
#define THEA_LOG Thea::LogInternal::LockedOutputStream(std::cout, THEA_LOG_STANDARD_PREFIX).getStream()

#ifdef THEA_DEBUG_BUILD
/**
 * Synchronized stream for debug messages, with a prefix indicating the time, source file and line number. Outputs a newline at
 * the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_DEBUG << a << b << c;</code>).
 *
 * Deactivated in release mode.
 */
#  define THEA_DEBUG Thea::LogInternal::LockedOutputStream(std::cout, THEA_LOG_STANDARD_PREFIX).getStream()
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
#define THEA_ERROR Thea::LogInternal::LockedOutputStream(std::cerr, THEA_LOG_STANDARD_PREFIX + "ERROR: ").getStream()

/**
 * Synchronized stream for warning messages, with a prefix indicating the time, source file and line number. Outputs a newline
 * at the end of every sequence of stream operations (i.e. after every stack such as <code>THEA_WARNING << a << b << c;</code>).
 */
#define THEA_WARNING Thea::LogInternal::LockedOutputStream(std::cerr, THEA_LOG_STANDARD_PREFIX + "WARNING: ").getStream()

#endif
