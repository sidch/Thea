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

#ifndef __Thea_Spinlock_hpp__
#define __Thea_Spinlock_hpp__

#include "AtomicInt32.hpp"

namespace Thea {

/**
 * A simple lock that can be shared by multiple threads for synchronization. Busy-waits for a lock, and unlocks immediately.
 * Implemented with a single atomic integer.
 */
class THEA_API Spinlock
{
  public:
    /** Constructor. The object is initially unlocked. */
    Spinlock() : atomic(0) {}

    /**
     * Acquire the lock by atomically setting the value of an integer to non-zero, busy-waiting until the lock is obtained.
     *
     * @see unlock()
     */
    void lock()
    {
      while (atomic.compareAndSet(0, 1) != 0) {}
    }

    /**
     * Release the lock by atomically setting the value of an integer to zero. This operation is immediate and does not wait.
     *
     * @see unlock()
     */
    void unlock()
    {
      atomic = 0;
    }

  private:
    AtomicInt32 atomic;

}; // class Spinlock

} // namespace Thea

#endif
