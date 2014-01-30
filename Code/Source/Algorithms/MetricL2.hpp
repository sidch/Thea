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

#ifndef __Thea_Algorithms_MetricL2_hpp__
#define __Thea_Algorithms_MetricL2_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../Triangle3.hpp"
#include "../VectorN.hpp"
#include "PointTraitsN.hpp"
#include "TransformedObject.hpp"
#include <boost/type_traits/is_pointer.hpp>
#include <boost/utility/enable_if.hpp>
#include <cmath>

namespace Thea {
namespace Algorithms {

/**
 * Helper class for MetricL2. Specializations of this class actually compute the metric. This is required because C++ does
 * unexpected things with specialized and overloaded function templates (see http://www.gotw.ca/publications/mill17.htm).
 *
 * @see MetricL2
 *
 * @todo Complete closestPoints() specializations for all supported types.
 */
template <typename A, typename B, long N, typename T, typename Enable = void>
struct /* THEA_API */ MetricL2Impl
{
  public:
    /** Measure the Euclidean (L2) distance between two objects. */
    static T distance(A const & a, B const & b)
    {
      return static_cast<T>(a.distance(b));
    }

    /**
     * Measure the square of the Euclidean (L2) distance between two objects, which is an efficiently computable (avoids the
     * square root) monotonic approximation to the true distance.
     */
    static T monotoneApproxDistance(A const & a, B const & b)
    {
      return static_cast<T>(a.squaredDistance(b));
    }

    /**
     * Find the closest pair of points between two objects.
     *
     * @return A monotonic approximation (the square of the Euclidean distance in this case) to the shortest distance between
     * the objects, i.e. the value of monotoneApproxDistance(\a a, \a b).
     */
    static T closestPoints(A const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
    {
      return static_cast<T>(a.closestPoints(b, cpa, cpb));
    }

}; // class MetricL2Impl

/**
 * Distances and closest pairs of points between various types of objects according to the Euclidean (L2) metric. When distances
 * will only be compared to one another to see which is larger, a monotonic function of the true distance works just as well and
 * may be more efficient to compute. In this case, the square of the Euclidean distance avoids a costly square root operation.
 * For this reason, it is better to use the monotoneApproxDistance() function instead of distance() wherever possible. The
 * functions computeMonotoneApprox() and invertMonotoneApprox() switch between the true and approximate distances.
 *
 * To add support for distances between custom types, add specializations (full or partial) of the helper class MetricL2Impl.
 * This class is required because C++ does unexpected things with specialized and overloaded function templates (see
 * http://www.gotw.ca/publications/mill17.htm). Note that to commutatively support distances between distinct types A and B, you
 * must specialize both MetricL2Impl<A, B> and MetricL2Impl<B, A>. Do <b>not</b> try specializing MetricL2::distance<N, T>() and
 * similar functions.
 *
 * MetricL2 defines the standard interface for a metric -- its interface must be supported by all other metric classes.
 *
 * @see MetricL2Impl
 */
class THEA_API MetricL2
{
  public:
    /** Compute a fixed monotone approximation (here, square) of a distance. */
    template <typename T> static T computeMonotoneApprox(T const & d) { return d * d; }

    /** Invert the fixed monotone approximation of a distance to get the true distance (here, by taking the square root). */
    template <typename T> static T invertMonotoneApprox(T const & d) { return std::sqrt(d); }

    /**
     * Measure the Euclidean (L2) distance between two objects via the helper class MetricL2Impl. Add specializations of
     * MetricL2Impl as required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     */
    template <long N, typename T, typename A, typename B> static T distance(A const & a, B const & b)
    { return MetricL2Impl<A, B, N, T>::distance(a, b); }

    /**
     * Measure the square of the Euclidean (L2) distance between two objects, which is an efficiently computable (avoids the
     * square root) monotonic approximation to the true distance. Add specializations of the helper class MetricL2Impl as
     * required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     */
    template <long N, typename T, typename A, typename B> static T monotoneApproxDistance(A const & a, B const & b)
    { return MetricL2Impl<A, B, N, T>::monotoneApproxDistance(a, b); }

    /**
     * Find the closest pair of points between two objects, via the helper class MetricL2Impl. Add specializations of
     * MetricL2Impl as required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     *
     * @return A monotonic approximation (the square of the Euclidean distance in this case) to the shortest distance between
     * the objects, i.e. the value of monotoneApproxDistance(\a a, \a b).
     */
    template <long N, typename T, typename A, typename B>
    static T closestPoints(A const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
    { return MetricL2Impl<A, B, N, T>::closestPoints(a, b, cpa, cpb); }

}; // class MetricL2

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename A, typename B, long N, typename T>
struct /* THEA_API */ MetricL2Impl<A *, B *, N, T>
{
  public:
    static T distance(A const * a, B const * b) { return MetricL2::distance<N, T>(*a, *b); }
    static T monotoneApproxDistance(A const * a, B const * b) { return MetricL2::monotoneApproxDistance<N, T>(*a, *b); }
    static T closestPoints(A const * a, B const * b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
    { return MetricL2::closestPoints<N, T>(*a, *b, cpa, cpb); }
};

template <typename A, typename B, long N, typename T>
struct /* THEA_API */ MetricL2Impl<A, B *, N, T>
{
  public:
    static T distance(A const & a, B const * b) { return MetricL2::distance<N, T>(a, *b); }
    static T monotoneApproxDistance(A const & a, B const * b) { return MetricL2::monotoneApproxDistance<N, T>(a, *b); }
    static T closestPoints(A const & a, B const * b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
    { return MetricL2::closestPoints<N, T>(a, *b, cpa, cpb); }
};

template <typename A, typename B, long N, typename T>
struct /* THEA_API */ MetricL2Impl<A *, B, N, T>
{
  public:
    static T distance(A const * a, B const & b) { return MetricL2::distance<N, T>(*a, b); }
    static T monotoneApproxDistance(A const * a, B const & b) { return MetricL2::monotoneApproxDistance<N, T>(*a, b); }
    static T closestPoints(A const * a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
    { return MetricL2::closestPoints<N, T>(*a, b, cpa, cpb); }
};

//=============================================================================================================================
// Default specializations
//=============================================================================================================================

namespace MetricL2Internal {

template <typename T> struct TransformedObjectCheck { static bool const value = false; };

template <typename ObjT, typename TransT> struct TransformedObjectCheck< TransformedObject<ObjT, TransT> >
{
  static bool const value = true;
};

} // namespace MetricL2Internal

// Both objects are points
template <typename A, typename B, long N, typename T>
struct MetricL2Impl<A, B, N, T, typename boost::enable_if_c< IsPointN<A, N>::value && IsPointN<B, N>::value >::type>
{
  static T distance(A const & a, B const & b)
  { return (PointTraitsN<A, N, T>::getPosition(a) - PointTraitsN<B, N, T>::getPosition(b)).length(); }

  static T monotoneApproxDistance(A const & a, B const & b)
  { return (PointTraitsN<A, N, T>::getPosition(a) - PointTraitsN<B, N, T>::getPosition(b)).squaredLength(); }

  static T closestPoints(A const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  {
    cpa = PointTraitsN<A, N, T>::getPosition(a);
    cpb = PointTraitsN<B, N, T>::getPosition(b);
    return (cpa - cpb).squaredLength();
  }
};

// Only the second object is a point
template <typename A, typename B, long N, typename T>
struct MetricL2Impl<A, B, N, T, typename boost::enable_if_c< !IsPointN<A, N>::value
                                                          && !MetricL2Internal::TransformedObjectCheck<A>::value
                                                          && !boost::is_pointer<A>::value
                                                          && IsPointN<B, N>::value >::type>
{
  static T distance(A const & a, B const & b) { return a.distance(PointTraitsN<B, N, T>::getPosition(b)); }

  static T monotoneApproxDistance(A const & a, B const & b)
  {
    VectorN<N, T> p = PointTraitsN<B, N, T>::getPosition(b);
    return (a.closestPoint(p) - p).squaredLength();
  }

  static T closestPoints(A const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  {
    cpb = PointTraitsN<B, N, T>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredLength();
  }
};

// Only the first object is a point
template <typename A, typename B, long N, typename T>
struct MetricL2Impl<A, B, N, T, typename boost::enable_if_c< IsPointN<A, N>::value
                                                          && !IsPointN<B, N>::value
                                                          && !MetricL2Internal::TransformedObjectCheck<B>::value
                                                          && !boost::is_pointer<B>::value >::type>
{
  static T distance(A const & a, B const & b)
  { return MetricL2Impl<B, A, N, T>::distance(b, a); }

  static T monotoneApproxDistance(A const & a, B const & b)
  { return MetricL2Impl<B, A, N, T>::monotoneApproxDistance(b, a); }

  static T closestPoints(A const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  { return MetricL2Impl<B, A, N, T>::closestPoints(b, a, cpb, cpa); }
};

//=============================================================================================================================
// Specializations where the distance() member function is not defined in both classes
//=============================================================================================================================

template <long N, typename T>
struct MetricL2Impl<AxisAlignedBoxN<N, T>, BallN<N, T>, N, T>
{
  static T distance(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b) { return b.squaredDistance(a); }
  static T closestPoints(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  { return b.closestPoints(a, cpb, cpa); }
};

template <typename VertexTripleT, typename T>
struct MetricL2Impl<BallN<3, T>, Triangle3<VertexTripleT>, 3, T>
{
  static T distance(BallN<3, T> const & a, Triangle3<VertexTripleT> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(BallN<3, T> const & a, Triangle3<VertexTripleT> const & b) { return b.squaredDistance(a); }
  static T closestPoints(BallN<3, T> const & a, Triangle3<VertexTripleT> const & b, VectorN<3, T> & cpa, VectorN<3, T> & cpb)
  { return b.closestPoints(a, cpb, cpa); }
};

} // namespace Algorithms
} // namespace Thea

#include "MetricL2_Transform.hpp"

#endif
