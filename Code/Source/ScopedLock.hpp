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

#ifndef __Thea_ScopedLock_hpp__
#define __Thea_ScopedLock_hpp__

namespace Thea {

/** A wrapper for a shared lock object that acquires the lock upon construction and releases it upon destruction. */
template <typename LockT>
class ScopedLock
{
  public:
    typedef LockT Lock;  ///< The wrapped lock class.

    /** Constructor, which acquires the lock. The wrapped lock object must persist as long as this wrapper does. */
    ScopedLock(LockT * lock_) : lock(lock_)
    {
      if (lock)
        lock->lock();
    }

    /** Destructor, which releases the lock. */
    ~ScopedLock()
    {
      if (lock)
        lock->unlock();
    }

  private:
    Lock * lock;

}; // class ScopedLock

} // namespace Thea

#endif
