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

#include "System.hpp"
#include <thread>

#ifdef THEA_WINDOWS
#  include <windows.h>
#else
#  include <unistd.h>
#endif

namespace Thea {

System::System()
{}

intx
System::concurrency()
{
  intx cc = (intx)std::thread::hardware_concurrency();
  if (cc <= 0)
    return 1;  // operate in single-threaded mode as fallback
  else
    return cc;
}

void
System::sleep(intx ms)
{
#ifdef THEA_WINDOWS
  Sleep(static_cast<DWORD>(ms));
#else
  usleep(static_cast<useconds_t>(ms * 1000));
#endif
}

} // namespace Thea
