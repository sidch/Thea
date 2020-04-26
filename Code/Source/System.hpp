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

#ifndef THEA_WINDOWS
#  include <sys/time.h>
#endif

namespace Thea {

/** Low-level system information and profiling functions. */
class THEA_API System
{
  public:
    /** Get the hardware concurrency (approximate number of thread contexts). */
    static intx concurrency();

    /** Get the machine endianness. */
    static Endianness endianness()
    {
      return Endianness::machine();
    }

    /** Pause the current thread for a given number of milliseconds. */
    static void sleep(intx ms);

    /**
     * The actual time (measured in seconds since Jan 1 1970 midnight). Adjusted for local timezone and daylight savings time.
     * This is as accurate and fast as getCycleCount().
     */
    static double time();

    /**
     * Begin a timing operation. To count the number of cycles a given operation takes:
     *
     * <pre>
     *   uint64 count;
     *   System::beginCycleCount(count);
     *   ...
     *   System::endCycleCount(count);
     *   // count now contains the cycle count for the intervening operation.
     * </pre>
     *
     * @see endCycleCount()
     */
    static void beginCycleCount(uint64 & cycle_count);

    /**
     * End a timing operation.
     *
     * @see beginCycleCount()
     */
    static void endCycleCount(uint64 & cycle_count);

  private:
    /** Constructor. */
    System();

    /** Called to initialize timing operations when the object is constructed. */
    void initTime();

    /** Get the current cycle count. */
    static uint64 getCycleCount();

    /** Singleton instance. */
    static System & instance()
    {
      static System s; return s;
    }

#ifdef THEA_WINDOWS
    int64          m_start;
    int64          m_counterFrequency;
#else
    struct timeval m_start;
#endif

#ifdef THEA_MAC
    int32          m_OSXCPUSpeed;
    double         m_secondsPerNS;
#endif

    double         m_realWorldGetTickTime0;

}; // class System

} // namespace Thea

#endif
