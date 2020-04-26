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
template <typename T, int N, typename ScalarT>
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
    template <int N, typename ScalarT, typename T>
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
    template <int N, typename ScalarT, typename T>
    static RayIntersectionN<N, ScalarT> rayIntersection(RayN<N, ScalarT> const & ray, T const & obj, ScalarT max_time = -1)
    {
      return RayIntersectionTesterImpl<T, N, ScalarT>::rayIntersection(ray, obj, max_time);
    }

}; // class RayIntersectionTester

// Support for pointer types
template <typename T, int N, typename ScalarT>
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
