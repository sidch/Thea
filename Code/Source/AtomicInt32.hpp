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

 @file AtomicInt32.h

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2005-09-01
 @edited  2006-06-21
 */

#ifndef __Thea_AtomicInt32_hpp__
#define __Thea_AtomicInt32_hpp__

#include "Platform.hpp"
#include "NumericType.hpp"

#if defined(THEA_MAC)
#  include <libkern/OSAtomic.h>
#elif defined(THEA_WINDOWS)
#  include <windows.h>
#endif

namespace Thea {

/**
 * An integer that may safely be used on different threads without external locking. On Win32, Linux, FreeBSD, and Mac OS X this
 * is implemented without locks. Derived from the G3D library: http://g3d.sourceforge.net
 */
class THEA_API AtomicInt32
{
  private:
#if defined(THEA_WINDOWS)
    typedef LONG ImplT;
    volatile ImplT m_value;
#elif defined(THEA_MAC)
    typedef int32_t ImplT;
    ImplT m_value;
#else
    typedef int32 ImplT;
    volatile ImplT m_value;
#endif

  public:
    /** Initial value is undefined. */
    AtomicInt32() {}

    /** Atomic initialization */
    explicit AtomicInt32(int32 x)
    {
      m_value = (ImplT)x;
    }

    /** Atomic copy constructor. */
    AtomicInt32(AtomicInt32 const & x)
    {
      m_value = x.m_value;
    }

    /** Atomic set */
    AtomicInt32 const & operator=(int32 x)
    {
      m_value = (ImplT)x;
      return *this;
    }

    /** Atomic set */
    void operator=(AtomicInt32 const & x)
    {
      m_value = x.m_value;
    }

    /** Get the current value */
    int32 value() const
    {
      return m_value;
    }

    /** Increment by \a x and return the old value, before the addition. */
    int32 add(int32 const x)
    {
#if defined(THEA_WINDOWS)
      return InterlockedExchangeAdd(&m_value, x);
#elif defined(THEA_LINUX) || defined(THEA_FREEBSD)
      int32 old;
      asm volatile ("lock; xaddl %0,%1"
                    : "=r"(old), "=m"(m_value) /* outputs */
                    : "0"(x), "m"(m_value)   /* inputs */
                    : "memory", "cc");
      return old;
#elif defined(THEA_MAC)
      int32 old = m_value;
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
      OSAtomicAdd32(x, &m_value);
#  pragma clang diagnostic pop
      return old;
#endif
    }

    /** Decrement by \a x and return the old value, before the subtraction. */
    int32 sub(int32 const x)
    {
      return add(-x);
    }

    /** Increment by 1. */
    void increment()
    {
#if defined(THEA_WINDOWS)
      // Note: returns the newly incremented value
      InterlockedIncrement(&m_value);
#elif defined(THEA_LINUX) || defined(THEA_FREEBSD)
      add(1);
#elif defined(THEA_MAC)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
      // Note: returns the newly incremented value
      OSAtomicIncrement32(&m_value);
#  pragma clang diagnostic pop
#endif
    }

    /**
     * Decrement by 1.
     *
     * @return Zero if the result is zero after decrement, non-zero otherwise.
     */
    int32 decrement()
    {
#if defined(THEA_WINDOWS)
      // Note: returns the newly decremented value
      return InterlockedDecrement(&m_value);
#elif defined(THEA_LINUX)  || defined(THEA_FREEBSD)
      unsigned char nz;
      asm volatile ("lock; decl %1;\n\t"
                    "setnz %%al"
                    : "=a" (nz)
                    : "m" (m_value)
                    : "memory", "cc");
      return nz;
#elif defined(THEA_MAC)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
      // Note: returns the newly decremented value
      return OSAtomicDecrement32(&m_value);
#  pragma clang diagnostic pop
#endif
    }

    /**
     * Atomic test-and-set: if <code>*this == comperand</code> then <code>*this := exchange</code>, else do nothing. In both
     * cases, returns the old value of <code>*this</code>.
     *
     * Performs an atomic comparison of the integer with the \a comperand value. If this is equal to the comperand, the
     * \a exchange value is stored in this. Otherwise, no operation is performed.
     *
     * @note Under VC6 the sign bit may be lost.
     */
    int32 compareAndSet(int32 comperand, int32 exchange)
    {
#if defined(THEA_WINDOWS)
      return InterlockedCompareExchange(&m_value, exchange, comperand);
#elif defined(THEA_LINUX) || defined(THEA_FREEBSD) || defined(THEA_MAC)
      // Based on Apache Portable Runtime
      // http://koders.com/c/fid3B6631EE94542CDBAA03E822CA780CBA1B024822.aspx
      int32 ret;
      asm volatile ("lock; cmpxchgl %1, %2"
                    : "=a" (ret)
                    : "r" (exchange), "m" (m_value), "0"(comperand)
                    : "memory", "cc");
      return ret;
      // Note that OSAtomicCompareAndSwap32 does not return a useful value for us
      // so it can't satisfy the cmpxchgl contract.
#endif
    }

}; // class AtomicInt32

} // namespace Thea

#endif
