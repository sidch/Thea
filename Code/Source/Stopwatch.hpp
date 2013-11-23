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

 @file Stopwatch.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2005-10-05
 @edited  2009-05-10

 Copyright 2000-2009, Morgan McGuire.
 All rights reserved.
*/

#ifndef __Thea_Stopwatch_hpp__
#define __Thea_Stopwatch_hpp__

#include "Common.hpp"

namespace Thea {

/**
 * Accurately measure durations and framerates.
 *
 * Example 1: For profiling code in the context of a rendering loop:
 * <pre>
 *   sw.tick();
 *   ...timed code...
 *   sw.tock();
 *
 *   printf("%f\n", sw.smoothFPS());
 * </pre>
 *
 * Example 2: For profiling pieces of a sequence:
 * <pre>
 *   Stopwatch sw;
 *   slowOperation();
 *   sw.after("slowOperation");
 *   kdTree.balance();
 *   sw.after("Balance tree");
 * </pre>
 *
 * Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API Stopwatch
{
  private:
    std::string             m_name;
    double                  startTime;
    std::string             prevMark;
    double                  prevTime;

    /** True between tick and tock */
    bool                    inBetween;

    /** The initial cycle count. */
    uint64                  cycleStart;

    /** The time at which tick was called. */
    double                  timeStart;

    /** The time at which the previous tock was called, -1 if never. */
    double                  lastTockTime;

    double                  lastDuration;
    int64                   lastCycleCount;

    /** Frames per second. */
    double                  m_fps;

    /** Weighted fps */
    double                  emwaFPS;
    double                  m_smoothFPS;

    /** Weighted duration */
    double                  emwaDuration;

    /** The overhead for calling into the class. */
    int64                   cycleOverhead;

    /** Called from the constructor. */
    void computeOverhead();

  public:
    /** Constructor. */
    Stopwatch(std::string const & name = "Stopwatch");

    /** Amount of time between the most recent tick() and tock() calls. 0 if tick() has never been called. */
    double elapsedTime() const
    {
      return lastDuration;
    }

    /**
     * Time-smoothed value that is stable to the nearest 1%. This is useful if you are displaying elapsed time in real time and
     * want a stable number.
     */
    double smoothElapsedTime() const
    {
      return emwaDuration;
    }

    /**
     * The elapsed cycle time between tick() and tock(). An attempt is made to factor out all tick/tock overhead, so that
     * back-to-back calls should return zero. Unreliable on non-x86 platforms.
     */
    uint64 elapsedCycles() const
    {
      return lastCycleCount;
    }

    /**
     * Get the number of times tick() was called per wall-clock second. If tick() is called once every frame, this measures
     * frames-per-second.
     */
    double fps() const
    {
      return m_fps;
    }

    /**
     * Time-smoothed value of fps that is stable to the nearest integer for fps > 10 and to the first decimal place for
     * fps <= 10. This is useful if you are displaying the frame rate in real-time and want a stable (readable) number.
     */
    double smoothFPS() const
    {
      return m_smoothFPS;
    }

    /** Begin a timing operation. */
    void tick();

    /** End a timing operation. */
    void tock();

    /** Reset the start time used by after() and the smoothed average measures.*/
    void reset();

    /**
     * Call after an operation has completed, with the name of the operation, to print a debug message listing the time since
     * the previous after() call.
     */
    void after(std::string const & s = "");

}; // class Stopwatch

} // namespace Thea

#endif
