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
