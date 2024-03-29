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
// First version: 2019
//
//============================================================================

//============================================================================
//
// ORIGINAL LICENSE
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
//============================================================================

#ifndef __Thea_ThreadGroup_hpp__
#define __Thea_ThreadGroup_hpp__

#include "Common.hpp"
#include "List.hpp"
#include <algorithm>
#include <mutex>
#include <shared_mutex>
#include <thread>

namespace Thea {

/**
 * Manages a group of threads and allows waiting for all of them to finish. Exists because <code>std::thread_group</code> did
 * not make it to C++11. Lifted straight from the Boost source code (1.67.0).
 */
class THEA_API ThreadGroup
{
  private:
    /** Disable copy construction. */
    ThreadGroup(ThreadGroup const &);

    /** Disable assignment. */
    ThreadGroup & operator=(ThreadGroup const &);

  public:
    /** Default constructor. */
    ThreadGroup() {}

    /** Destructor. */
    ~ThreadGroup()
    {
      for (auto it = threads.begin(), end = threads.end(); it != end; ++it)
        delete *it;
    }

    /** Check if the group contains the current thread. */
    bool containsThisThread() const
    {
      std::thread::id id = std::this_thread::get_id();
      std::shared_lock<std::shared_mutex> guard(m);

      for (auto it = threads.begin(), end = threads.end(); it != end; ++it)
      {
        if ((*it)->get_id() == id)
          return true;
      }

      return false;
    }

    /** Check if the group contains a particular thread. */
    bool containsThread(std::thread * thrd) const
    {
      if (thrd)
      {
        std::thread::id id = thrd->get_id();
        std::shared_lock<std::shared_mutex> guard(m);

        for (auto it = threads.begin(), end = threads.end(); it != end; ++it)
        {
          if ((*it)->get_id() == id)
            return true;
        }
      }

      return false;
    }

    /** Create a thread wrapping a functor and add it to the group. */
    template <typename F> std::thread * createThread(F thread_func)
    {
      std::lock_guard<std::shared_mutex> guard(m);
      std::unique_ptr<std::thread> new_thread(new std::thread(thread_func));
      threads.push_back(new_thread.get());
      return new_thread.release();
    }

    /** Add an existing thread to the group. */
    void addThread(std::thread * thrd)
    {
      if (thrd)
      {
        alwaysAssertM(!containsThread(thrd), "ThreadGroup: Thread is already in group");

        std::lock_guard<std::shared_mutex> guard(m);
        threads.push_back(thrd);
      }
    }

    /** Remove a thread from the group. */
    void removeThread(std::thread * thrd)
    {
      std::lock_guard<std::shared_mutex> guard(m);

      auto it = std::find(threads.begin(), threads.end(), thrd);
      if (it != threads.end())
        threads.erase(it);
    }

    /** Wait for all threads to finish. */
    void joinAll()
    {
      alwaysAssertM(!containsThisThread(), "ThreadGroup: Can't join the executing thread");
      std::shared_lock<std::shared_mutex> guard(m);

      for (auto it = threads.begin(), end = threads.end(); it != end; ++it)
      {
        if ((*it)->joinable())
          (*it)->join();
      }
    }

    size_t numThreads() const
    {
      std::shared_lock<std::shared_mutex> guard(m);
      return threads.size();
    }

  private:
    List<std::thread *> threads;  // Holds the set of threads.
    mutable std::shared_mutex m;  // Mutex allowing multiple-reader/single-writer access to the thread list.

}; // class ThreadGroup

} // namespace Thea

#endif
