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

#ifndef __Thea_Error_hpp__
#define __Thea_Error_hpp__

#include "Platform.hpp"
#include "BasicStringAlg.hpp"
#include "Log.hpp"
#include <stdexcept>
#include <string>

namespace Thea {

/**
 * An error class.
 */
typedef std::runtime_error Error;

/**
 * A fatal error class. Doesn't derive from std::exception so that it is not caught by one of the usual handlers.
 */
class THEA_API FatalError
{
  public:
    /** Constructor. */
    FatalError(std::string const & message_) : message(message_) {}

    /** Destructor. */
    ~FatalError() {}

    /** Assignment operator. */
    FatalError & operator=(FatalError const & src) { message = src.message; return *this; }

    /** The error message. */
    char const * what() const { return message.c_str(); }

  private:
    std::string message;

}; // class FatalError

/**
 * Generate a standard series of catch blocks that print a message on the specified stream and then perform an action. A fatal
 * error always generates a standard message on the error stream and the exception is rethrown, regardless of the specified
 * stream, message and action.
 *
 * Example:
 * <pre>
 *   THEA_STANDARD_CATCH_BLOCKS(continue;, WARNING, "A non-fatal error occurred in iteration %d", i)
 * </pre>
 */
#define THEA_STANDARD_CATCH_BLOCKS(action__, stream__, message__, ...) \
    catch (FatalError & e__) \
    { \
      THEA_ERROR << format( "A fatal error occurred (%s)", e__.what() ); \
      throw; \
    } \
    catch (Error & e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__.what() ); \
      action__ \
    } \
    catch (std::exception & e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__.what() ); \
      action__ \
    } \
    catch (std::string & e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__.c_str() ); \
      action__ \
    } \
    catch (const char * e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__ ); \
      action__ \
    } \
    catch (...) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (Unknown error)")).c_str(), __VA_ARGS__ ); \
      action__ \
    }

} // namespace Thea

#endif
