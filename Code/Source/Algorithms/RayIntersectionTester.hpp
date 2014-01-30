//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_RayIntersectionTester_hpp__
#define __Thea_Algorithms_RayIntersectionTester_hpp__

#include "../Common.hpp"
#include "../RayIntersectableN.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Helper class for RayIntersectionTester. Specializations of this class actually test for intersection. This is required
 * because C++ does unexpected things with specialized and overloaded function templates (see
 * http://www.gotw.ca/publications/mill17.htm).
 *
 * The default implementation works for objects implementing the RayIntersectable3 interface.
 *
 * @see RayIntersectionTester
 */
template <typename T, long N, typename ScalarT>
struct /* THEA_API */ RayIntersectionTesterImpl
{
  /**
   * Get the time taken for a ray to intersect an object.
   *
   * @param ray The ray to test for intersection with the object.
   * @param obj The object to test for intersection with the ray.
   * @param max_time Maximum allowable hit time, ignored if negative.
   */
  static ScalarT rayIntersectionTime(RayN<N, ScalarT> const & ray, T const & obj, ScalarT const & max_time)
  { return obj.rayIntersectionTime(ray, max_time); }

  /**
   * Get a description of the intersection point of a ray with an object.
   *
   * @param ray The ray to test for intersection with the object.
   * @param obj The object to test for intersection with the ray.
   * @param max_time Maximum allowable hit time, ignored if negative.
   */
  static RayIntersectionN<N, ScalarT> rayIntersection(RayN<N, ScalarT> const & ray, T const & obj, ScalarT const & max_time)
  { return obj.rayIntersection(ray, max_time); }

}; // struct RayIntersectionTesterImpl

/**
 * Ray intersection queries on objects.
 *
 * To add support for intersection queries between custom types, add specializations (full or partial) of the helper class
 * RayIntersectionTesterImpl. This class is required because C++ does unexpected things with specialized and overloaded function
 * templates (see http://www.gotw.ca/publications/mill17.htm). Do <b>not</b> try specializing
 * RayIntersectionTester::rayIntersectionTime() or RayIntersectionTester::rayIntersection().
 */
class THEA_API RayIntersectionTester
{
  public:
    /**
     * Get the time taken for a ray to intersect an object, via the helper class RayIntersectionTesterImpl. Add specializations
     * of RayIntersectionTesterImpl as required.
     *
     * @param ray The ray to test for intersection with the object.
     * @param obj The object to test for intersection with the ray.
     * @param max_time Maximum allowable hit time, ignored if negative.
     */
    template <long N, typename ScalarT, typename T>
    static ScalarT rayIntersectionTime(RayN<N, ScalarT> const & ray, T const & obj, ScalarT max_time = -1)
    {
      return RayIntersectionTesterImpl<T, N, ScalarT>::rayIntersectionTime(ray, obj, max_time);
    }

    /**
     * Get a description of the intersection point of a ray with an object, via the helper class RayIntersectionTesterImpl. Add
     * specializations of RayIntersectionTesterImpl as required.
     *
     * @param ray The ray to test for intersection with the object.
     * @param obj The object to test for intersection with the ray.
     * @param max_time Maximum allowable hit time, ignored if negative.
     */
    template <long N, typename ScalarT, typename T>
    static RayIntersectionN<N, ScalarT> rayIntersection(RayN<N, ScalarT> const & ray, T const & obj, ScalarT max_time = -1)
    {
      return RayIntersectionTesterImpl<T, N, ScalarT>::rayIntersection(ray, obj, max_time);
    }

}; // class RayIntersectionTester

// Support for pointer types
template <typename T, long N, typename ScalarT>
struct /* THEA_API */ RayIntersectionTesterImpl<T *, N, ScalarT>
{
  static ScalarT rayIntersectionTime(RayN<N, ScalarT> const & ray, T const * obj, ScalarT const & max_time)
  { return RayIntersectionTester::rayIntersectionTime<N, ScalarT>(ray, *obj, max_time); }

  static RayIntersectionN<N, ScalarT> rayIntersection(RayN<N, ScalarT> const & ray, T const * obj, ScalarT const & max_time)
  { return RayIntersectionTester::rayIntersection<N, ScalarT>(ray, *obj, max_time); }
};

} // namespace Algorithms
} // namespace Thea

#endif
