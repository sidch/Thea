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

 @file Stopwatch.cpp

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2005-10-05
 @edited  2009-03-14

 Copyright 2000-2009, Morgan McGuire.
 All rights reserved.
*/

#include "Stopwatch.hpp"
#include "System.hpp"
#include "Math.hpp"

namespace Thea {

Stopwatch::Stopwatch(const std::string & name)
: m_name(name),
  inBetween(false),
  lastTockTime(-1),
  lastDuration(0),
  lastCycleCount(0),
  m_fps(0),
  emwaFPS(0),
  m_smoothFPS(0),
  emwaDuration(0)
{
  computeOverhead();
  reset();
}

void
Stopwatch::computeOverhead()
{
  cycleOverhead = 0;
  tick();
  tock();
  cycleOverhead = elapsedCycles();
}

void
Stopwatch::tick()
{
  // This is 'alwaysAssertM' instead of 'debugAssertM' since people rarely profile in debug mode.
  alwaysAssertM(!inBetween, "Stopwatch: tick() called twice in a row");

  inBetween = true;

  // We read RDTSC twice here, but it is more abstract to implement this way and at least we're reading the cycle count last.
  timeStart = System::time();
  System::beginCycleCount(cycleStart);
}

void
Stopwatch::tock()
{
  System::endCycleCount(cycleStart);
  double now = System::time();
  lastDuration = now - timeStart;

  if (std::fabs(emwaDuration - lastDuration) > std::max(emwaDuration, lastDuration) * 0.50)
  {
    // Off by more than 50%
    emwaDuration = lastDuration;
  }
  else
  {
    emwaDuration = lastDuration * 0.05 + emwaDuration * 0.95;
  }

  lastCycleCount = cycleStart - cycleOverhead;

  if (lastCycleCount < 0)
  {
    lastCycleCount = 0;
  }

  if (lastTockTime != -1.0)
  {
    m_fps = 1.0 / (now - lastTockTime);

    double const BLEND = 0.01;
    emwaFPS = m_fps * BLEND + emwaFPS * (1.0 - BLEND);

    double const MAX_DISCREPANCY_PERCENTAGE = 0.25;
    if (std::fabs(emwaFPS - m_fps) > std::max(emwaFPS, m_fps) * MAX_DISCREPANCY_PERCENTAGE)
    {
      // The difference between emwa and m_fps is way off, so
      // update emwa directly.
      emwaFPS = m_fps * 0.20 + emwaFPS * 0.80;
    }

    // Update m_smoothFPS only when the value varies significantly.
    // We round so as to not mislead the user as to the accuracy of
    // the number.
    if (m_smoothFPS == 0)
    {
      m_smoothFPS = m_fps;
    }
    else if (emwaFPS <= 20)
    {
      if (std::fabs(m_smoothFPS - emwaFPS) > 0.75)
      {
        // Small number and display is off by more than 0.75; round to the nearest 0.1
        m_smoothFPS = floor(emwaFPS * 10.0 + 0.5) / 10.0;
      }
    }
    else if (std::fabs(m_smoothFPS - emwaFPS) > 1.25)
    {
      // Large number and display is off by more than 1.25; round to the nearest 1.0
      m_smoothFPS = std::floor(emwaFPS + 0.5);
    }
  }

  lastTockTime = now;

  alwaysAssertM(inBetween, "Stopwatch: tock() called without matching tick().");
  inBetween = false;
}

void
Stopwatch::reset()
{
  prevTime = startTime = System::time();
  prevMark = "start";
}

void
Stopwatch::after(std::string const & s)
{
  double now = System::time();
  THEA_CONSOLE << m_name << ": " << s << " - " << now - prevTime << "s since " << prevMark
               << " (" << now - startTime << "s since start)",
  prevTime = now;
  prevMark = s;
}

} // namespace Thea
