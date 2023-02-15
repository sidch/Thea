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

#ifndef __Thea_System_hpp__
#define __Thea_System_hpp__

#include "Common.hpp"
#include <chrono>

namespace Thea {

/** Low-level system information and profiling functions. */
class THEA_API System
{
  public:
    /** Get the hardware concurrency (approximate number of thread contexts). */
    static intx concurrency();

    /** Get the machine endianness. */
    static Endianness endianness() { return Endianness::machine(); }

    /** Pause the current thread for a given number of milliseconds. */
    static void sleep(intx ms);

    /** The time in seconds since Jan 1 1970 midnight. Adjusted for local timezone and daylight savings time. */
    static double time()
    {
      typedef std::chrono::duration<double> DSecs;

      return std::chrono::duration_cast<DSecs>(std::chrono::system_clock::now().time_since_epoch()).count();
    }

  private:
    /** Constructor. */
    System();

    /** Singleton instance. */
    static System & instance() { static System s; return s; }

}; // class System

} // namespace Thea

#endif
