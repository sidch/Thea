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

#include "Log.hpp"
#include "Common.hpp"
#include <ctime>
#include <iomanip>
#include <sstream>

namespace Thea {
namespace LogInternal {

Spinlock lock;

std::string
stripPathFromFilename(std::string const & full_path)
{
  std::size_t last_slash;
  last_slash = full_path.find_last_of("/\\");
  return (last_slash >= full_path.length()) ? full_path : full_path.substr(last_slash + 1);
}

std::string
currentDateTimeToString()
{
  std::time_t raw_time;
  std::tm const * time_info;

  std::time(&raw_time);  // get raw system time
  time_info = std::localtime(&raw_time);  // break it down into components w.r.t. the current locale
  std::ostringstream os;
  os << 1900 + time_info->tm_year << '-'
     << std::setw(2) << std::right << std::setfill('0') << 1 + time_info->tm_mon << '-'
     << std::setw(2) << std::right << std::setfill('0') << time_info->tm_mday << ' '
     << std::setw(2) << std::right << std::setfill('0') << time_info->tm_hour << ':'
     << std::setw(2) << std::right << std::setfill('0') << time_info->tm_min  << ':'
     << std::setw(2) << std::right << std::setfill('0') << time_info->tm_sec;
  return os.str();
}

} // namespace LogInternal
} // namespace Thea
