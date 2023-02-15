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
 ORIGINAL HEADER ([SC] Note: As of Feb 2023, this class has been completely
 rewritten to use C++17 standard library functions and remove unused
 functionality. Apart from some function signatures, it has no further
 connection with G3D.)

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
#include <chrono>

namespace Thea {

/** Measure duration between calls to tick() and tock(). */
class THEA_API Stopwatch
{
  private:
    using Clock = std::chrono::steady_clock;

    Clock::time_point time_start;
    double last_duration;
    bool in_between;

  public:
    /** Constructor. */
    Stopwatch() : last_duration(0), in_between(false) {}

    /** Begin a timing operation. */
    void tick()
    {
      if (in_between) { throw Error("Stopwatch: tick() called twice in a row"); }

      in_between = true;
      time_start = Clock::now();
    }

    /** End a timing operation. */
    void tock()
    {
      std::chrono::duration<double> secs = Clock::now() - time_start;
      last_duration = secs.count();

      if (!in_between) { throw Error("Stopwatch: tock() called without matching tick()"); }
      in_between = false;
    }

    /** Duration in seconds between the last tick() and tock() calls. Returns 0 if tick() has never been called. */
    double elapsedTime() const { return last_duration; }

}; // class Stopwatch

} // namespace Thea

#endif
