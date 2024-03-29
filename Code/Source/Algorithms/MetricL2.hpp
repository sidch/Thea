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

#ifndef __Thea_Algorithms_MetricL2_hpp__
#define __Thea_Algorithms_MetricL2_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../LineN.hpp"
#include "../LineSegmentN.hpp"
#include "../MatVec.hpp"
#include "../RayN.hpp"
#include "../TriangleN.hpp"
#include "PointTraitsN.hpp"
#include "TransformedObject.hpp"
#include <cmath>
#include <type_traits>

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
template <typename A, typename B, int N, typename T, typename Enable = void>
struct /* THEA_API */ MetricL2Impl
{
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
  static T closestPoints(A const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    return static_cast<T>(a.squaredDistance(b, &cpa, &cpb));
  }

}; // struct MetricL2Impl

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
    template <int N, typename T, typename A, typename B> static T distance(A const & a, B const & b)
    { return MetricL2Impl<A, B, N, T>::distance(a, b); }

    /**
     * Measure the square of the Euclidean (L2) distance between two objects, which is an efficiently computable (avoids the
     * square root) monotonic approximation to the true distance. Add specializations of the helper class MetricL2Impl as
     * required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     */
    template <int N, typename T, typename A, typename B> static T monotoneApproxDistance(A const & a, B const & b)
    { return MetricL2Impl<A, B, N, T>::monotoneApproxDistance(a, b); }

    /**
     * Find the closest pair of points between two objects, via the helper class MetricL2Impl. Add specializations of
     * MetricL2Impl as required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     *
     * @return A monotonic approximation (the square of the Euclidean distance in this case) to the shortest distance between
     * the objects, i.e. the value of monotoneApproxDistance(\a a, \a b).
     */
    template <int N, typename T, typename A, typename B>
    static T closestPoints(A const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
    { return MetricL2Impl<A, B, N, T>::closestPoints(a, b, cpa, cpb); }

}; // class MetricL2

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ MetricL2Impl<A *, B *, N, T>
{
  static T distance(A const * a, B const * b) { return MetricL2::distance<N, T>(*a, *b); }
  static T monotoneApproxDistance(A const * a, B const * b) { return MetricL2::monotoneApproxDistance<N, T>(*a, *b); }
  static T closestPoints(A const * a, B const * b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return MetricL2::closestPoints<N, T>(*a, *b, cpa, cpb); }
};

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ MetricL2Impl<A, B *, N, T>
{
  static T distance(A const & a, B const * b) { return MetricL2::distance<N, T>(a, *b); }
  static T monotoneApproxDistance(A const & a, B const * b) { return MetricL2::monotoneApproxDistance<N, T>(a, *b); }
  static T closestPoints(A const & a, B const * b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return MetricL2::closestPoints<N, T>(a, *b, cpa, cpb); }
};

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ MetricL2Impl<A *, B, N, T>
{
  static T distance(A const * a, B const & b) { return MetricL2::distance<N, T>(*a, b); }
  static T monotoneApproxDistance(A const * a, B const & b) { return MetricL2::monotoneApproxDistance<N, T>(*a, b); }
  static T closestPoints(A const * a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
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
template <typename A, typename B, int N, typename T>
struct MetricL2Impl<A, B, N, T, typename std::enable_if< IsNonReferencedPointN<A, N>::value
                                                      && IsNonReferencedPointN<B, N>::value >::type>
{
  static T distance(A const & a, B const & b)
  { return (PointTraitsN<A, N, T>::getPosition(a) - PointTraitsN<B, N, T>::getPosition(b)).norm(); }

  static T monotoneApproxDistance(A const & a, B const & b)
  { return (PointTraitsN<A, N, T>::getPosition(a) - PointTraitsN<B, N, T>::getPosition(b)).squaredNorm(); }

  static T closestPoints(A const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    cpa = PointTraitsN<A, N, T>::getPosition(a);
    cpb = PointTraitsN<B, N, T>::getPosition(b);
    return (cpa - cpb).squaredNorm();
  }
};

// Only the second object is a point
template <typename A, typename B, int N, typename T>
struct MetricL2Impl<A, B, N, T, typename std::enable_if< !IsPointN<A, N>::value
                                                      && !MetricL2Internal::TransformedObjectCheck<A>::value
                                                      && !std::is_pointer<A>::value
                                                      && IsNonReferencedPointN<B, N>::value >::type>
{
  static T distance(A const & a, B const & b) { return a.distance(PointTraitsN<B, N, T>::getPosition(b)); }

  static T monotoneApproxDistance(A const & a, B const & b)
  {
    Vector<N, T> p = PointTraitsN<B, N, T>::getPosition(b);
    return (a.closestPoint(p) - p).squaredNorm();
  }

  static T closestPoints(A const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    cpb = PointTraitsN<B, N, T>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredNorm();
  }
};

// Only the first object is a point
template <typename A, typename B, int N, typename T>
struct MetricL2Impl<A, B, N, T, typename std::enable_if< IsNonReferencedPointN<A, N>::value
                                                      && !IsPointN<B, N>::value
                                                      && !MetricL2Internal::TransformedObjectCheck<B>::value
                                                      && !std::is_pointer<B>::value >::type>
{
  static T distance(A const & a, B const & b)
  { return MetricL2Impl<B, A, N, T>::distance(b, a); }

  static T monotoneApproxDistance(A const & a, B const & b)
  { return MetricL2Impl<B, A, N, T>::monotoneApproxDistance(b, a); }

  static T closestPoints(A const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return MetricL2Impl<B, A, N, T>::closestPoints(b, a, cpb, cpa); }
};

//=============================================================================================================================
// Specializations where the distance()/squaredDistance() member functions are not defined in both classes
//=============================================================================================================================

template <int N, typename T>
struct MetricL2Impl<LineN<N, T>, AxisAlignedBoxN<N, T>, N, T>
{
  static T distance(LineN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(LineN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return b.squaredDistance(a); }
  static T closestPoints(LineN<N, T> const & a, AxisAlignedBoxN<N, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename T>
struct MetricL2Impl<LineSegmentN<N, T>, AxisAlignedBoxN<N, T>, N, T>
{
  static T distance(LineSegmentN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(LineSegmentN<N, T> const & a, AxisAlignedBoxN<N, T> const & b)
  { return b.squaredDistance(a); }
  static T closestPoints(LineSegmentN<N, T> const & a, AxisAlignedBoxN<N, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename T>
struct MetricL2Impl<RayN<N, T>, AxisAlignedBoxN<N, T>, N, T>
{
  static T distance(RayN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(RayN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return b.squaredDistance(a); }
  static T closestPoints(RayN<N, T> const & a, AxisAlignedBoxN<N, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename T>
struct MetricL2Impl<AxisAlignedBoxN<N, T>, BallN<N, T>, N, T>
{
  static T distance(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b) { return b.squaredDistance(a); }
  static T closestPoints(AxisAlignedBoxN<N, T> const & a, BallN<N, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename VertexTripleT, typename T>
struct MetricL2Impl<LineN<N, T>, TriangleN<N, VertexTripleT, T>, N, T>
{
  static T distance(LineN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(LineN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b)
  { return b.squaredDistance(a); }
  static T closestPoints(LineN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename VertexTripleT, typename T>
struct MetricL2Impl<LineSegmentN<N, T>, TriangleN<N, VertexTripleT, T>, N, T>
{
  static T distance(LineSegmentN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(LineSegmentN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b)
  { return b.squaredDistance(a); }
  static T closestPoints(LineSegmentN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b,
                         Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename VertexTripleT, typename T>
struct MetricL2Impl<RayN<N, T>, TriangleN<N, VertexTripleT, T>, N, T>
{
  static T distance(RayN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(RayN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b)
  { return b.squaredDistance(a); }
  static T closestPoints(RayN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

template <int N, typename VertexTripleT, typename T>
struct MetricL2Impl<BallN<N, T>, TriangleN<N, VertexTripleT, T>, N, T>
{
  static T distance(BallN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b) { return b.distance(a); }
  static T monotoneApproxDistance(BallN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b)
  { return b.squaredDistance(a); }
  static T closestPoints(BallN<N, T> const & a, TriangleN<N, VertexTripleT, T> const & b,
                         Vector<N, T> & cpa, Vector<N, T> & cpb)
  { return b.squaredDistance(a, &cpb, &cpa); }
};

} // namespace Algorithms
} // namespace Thea

#include "MetricL2_Transform.hpp"

#endif
