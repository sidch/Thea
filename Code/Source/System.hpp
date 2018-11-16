//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Stanford University
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
    static long concurrency();

    /** Get the machine endianness. */
    static Endianness endianness()
    {
      return Endianness::machine();
    }

    /** Pause the current thread for a given number of milliseconds. */
    static void sleep(long ms);

    /**
     * The actual time (measured in seconds since Jan 1 1970 midnight). Adjusted for local timezone and daylight savings time.
     * This is as accurate and fast as getCycleCount().
     */
    static double time();

    /**
     * Begin a timing operation. To count the number of cycles a given operation takes:
     *
     * <pre>
     *   unsigned long count;
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
