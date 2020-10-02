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

#ifndef __Thea_Error_hpp__
#define __Thea_Error_hpp__

#include "Platform.hpp"
#include "BasicStringAlg.hpp"
#include "Log.hpp"
#include <stdexcept>
#include <string>

namespace Thea {

/** An error class. */
typedef std::runtime_error Error;

/** A fatal error class. Doesn't derive from std::exception so that it is not caught by one of the usual handlers. */
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
    catch (FatalError const & e__) \
    { \
      THEA_ERROR << format( "A fatal error occurred (%s)", e__.what() ); \
      throw; \
    } \
    catch (std::exception const & e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__.what() ); \
      action__ \
    } \
    catch (std::string const & e__) \
    { \
      THEA_ ## stream__ << format( (message__ + std::string(" (%s)")).c_str(), __VA_ARGS__, e__.c_str() ); \
      action__ \
    } \
    catch (char const * e__) \
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
