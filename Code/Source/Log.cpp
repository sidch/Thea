//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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

#include "Log.hpp"
#include <ctime>
#include <iomanip>
#include <sstream>

namespace Thea {
namespace LogInternal {

std::string
stripPathFromFilename(std::string const & fullPath)
{
  std::size_t lastSlash;
  lastSlash = fullPath.find_last_of("/\\");
  return (lastSlash >= fullPath.length()) ? fullPath : fullPath.substr(lastSlash + 1);
}

std::string
currentDateTimeToString()
{
  std::time_t rawTime;
  std::tm const * timeInfo;

  std::time(&rawTime);  // get raw system time
  timeInfo = std::localtime(&rawTime);  // break it down into components w.r.t. the current locale
  std::ostringstream os;
  os << 1900 + timeInfo->tm_year << '-'
     << std::setw(2) << std::right << std::setfill('0') << timeInfo->tm_mon << '-'
     << std::setw(2) << std::right << std::setfill('0') << timeInfo->tm_mday << ' '
     << std::setw(2) << std::right << std::setfill('0') << timeInfo->tm_hour << ':'
     << std::setw(2) << std::right << std::setfill('0') << timeInfo->tm_min  << ':'
     << std::setw(2) << std::right << std::setfill('0') << timeInfo->tm_sec;
  return os.str();
}

LockedOutputStream::LockedOutputStream(std::ostream & stream_)
: stream(&stream_)
{
  setOutputLock(true);
}

LockedOutputStream::LockedOutputStream(std::ostream & stream_, std::string const & prefix)
: stream(&stream_)
{
  setOutputLock(true);
  getStream() << prefix;
}

LockedOutputStream::~LockedOutputStream()
{
  getStream() << std::endl;
  setOutputLock(false);
}

void
LockedOutputStream::setOutputLock(bool value)
{
  static G3D::GMutex mutex;

  if (value)
    mutex.lock();
  else
    mutex.unlock();
}

} // namespace LogInternal
} // namespace Thea
