//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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
#include "../AxisAlignedBox3.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../Ball3.hpp"
#include "../Line3.hpp"
#include "../LineSegment3.hpp"
#include "../Ray3.hpp"
#include "../Triangle3.hpp"
#include "../VectorN.hpp"
#include "PointTraitsN.hpp"
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
template <typename A, typename B, typename Enable = void>
struct /* THEA_API */ MetricL2Impl
{
  private:
    typedef char UnspecifiedPointT;

  public:
    /** Measure the Euclidean (L2) distance between two objects. */
    static double distance(A const & a, B const & b);

    /**
     * Measure the square of the Euclidean (L2) distance between two objects, which is an efficiently computable (avoids the
     * square root) monotonic approximation to the true distance.
     */
    static double monotoneApproxDistance(A const & a, B const & b);

    /**
     * Find the closest pair of points between two objects.
     *
     * @return A monotonic approximation (the square of the Euclidean distance in this case) to the shortest distance between
     * the objects, i.e. the value of monotoneApproxDistance(\a a, \a b).
     */
    static double closestPoints(A const & a, B const & b, UnspecifiedPointT & cpa, UnspecifiedPointT & cpb);

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
 * must specialize both MetricL2Impl<A, B> and MetricL2Impl<B, A>. Do <b>not</b> try specializing MetricL2::distance() and
 * similar functions.
 *
 * MetricL2 defines the standard interface for a metric -- its interface must be supported by all other metric classes.
 *
 * @see MetricL2Impl
 */
class /* THEA_API */ MetricL2  // static members defined in header, so can't use dllimport
{
  public:
    /** Compute a fixed monotone approximation (here, square) of a distance. */
    static double computeMonotoneApprox(double d) { return d * d; }

    /** Invert the fixed monotone approximation of a distance to get the true distance (here, by taking the square root). */
    static double invertMonotoneApprox(double d) { return std::sqrt(d); }

    /**
     * Measure the Euclidean (L2) distance between two objects via the helper class MetricL2Impl. Add specializations of
     * MetricL2Impl as required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     */
    template <typename A, typename B> static double distance(A const & a, B const & b)
    { return MetricL2Impl<A, B>::distance(a, b); }

    /**
     * Measure the square of the Euclidean (L2) distance between two objects, which is an efficiently computable (avoids the
     * square root) monotonic approximation to the true distance. Add specializations of the helper class MetricL2Impl as
     * required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     */
    template <typename A, typename B> static double monotoneApproxDistance(A const & a, B const & b)
    { return MetricL2Impl<A, B>::monotoneApproxDistance(a, b); }

    /**
     * Find the closest pair of points between two objects, via the helper class MetricL2Impl. Add specializations of
     * MetricL2Impl s required for specific types of objects. Note that if types A and B differ, you must explicitly specialize
     * both MetricL2Impl<A, B> and MetricL2Impl<B, A>.
     *
     * @return A monotonic approximation (the square of the Euclidean distance in this case) to the shortest distance between
     * the objects, i.e. the value of monotoneApproxDistance(\a a, \a b).
     */
    template <typename A, typename B, typename PointT>
    static double closestPoints(A const & a, B const & b, PointT & cpa, PointT & cpb)
    { return MetricL2Impl<A, B>::closestPoints(a, b, cpa, cpb); }

}; // class MetricL2

// Support for pointer types
template <typename A, typename B>
struct /* THEA_API */ MetricL2Impl<A *, B *>
{
  public:
    static double distance(A const * a, B const * b) { return MetricL2::distance(*a, *b); }
    static double monotoneApproxDistance(A const * a, B const * b) { return MetricL2::monotoneApproxDistance(*a, *b); }
    template <typename PointT> static double closestPoints(A const * a, B const * b, PointT & cpa, PointT & cpb)
    { return MetricL2::closestPoints(*a, *b, cpa, cpb); }
};

template <typename A, typename B>
struct /* THEA_API */ MetricL2Impl<A, B *>
{
  public:
    static double distance(A const & a, B const * b) { return MetricL2::distance(a, *b); }
    static double monotoneApproxDistance(A const & a, B const * b) { return MetricL2::monotoneApproxDistance(a, *b); }
    template <typename PointT> static double closestPoints(A const & a, B const * b, PointT & cpa, PointT & cpb)
    { return MetricL2::closestPoints(a, *b, cpa, cpb); }
};

template <typename A, typename B>
struct /* THEA_API */ MetricL2Impl<A *, B>
{
  public:
    static double distance(A const * a, B const & b) { return MetricL2::distance(*a, b); }
    static double monotoneApproxDistance(A const * a, B const & b) { return MetricL2::monotoneApproxDistance(*a, b); }
    template <typename PointT> static double closestPoints(A const * a, B const & b, PointT & cpa, PointT & cpb)
    { return MetricL2::closestPoints(*a, b, cpa, cpb); }
};

// Default specializations
template <>
struct /* THEA_API */ MetricL2Impl<float, float>
{
  static double distance(float a, float b) { return std::fabs(a - b); }
  static double monotoneApproxDistance(float a, float b) { float x = a - b; return x * x; }

  static double closestPoints(float a, float b, float & cpa, float & cpb)
  { cpa = a; cpb = b; return monotoneApproxDistance(a, b); }
};

template <>
struct /* THEA_API */ MetricL2Impl<double, double>
{
  static double distance(double a, double b) { return std::fabs(a - b); }
  static double monotoneApproxDistance(double a, double b) { double x = a - b; return x * x; }

  static double closestPoints(double a, double b, double & cpa, double & cpb)
  { cpa = a; cpb = b; return monotoneApproxDistance(a, b); }
};

template <typename A, typename B>
struct /* THEA_API */ MetricL2Impl<A, B, typename boost::enable_if_c< IsPointN<A, 2>::value && IsPointN<B, 2>::value >::type>
{
  static double distance(A const & a, B const & b)
  { return (PointTraitsN<A, 2>::getPosition(a) - PointTraitsN<B, 2>::getPosition(b)).length(); }

  static double monotoneApproxDistance(A const & a, B const & b)
  { return (PointTraitsN<A, 2>::getPosition(a) - PointTraitsN<B, 2>::getPosition(b)).squaredLength(); }

  static double closestPoints(A const & a, B const & b, Vector2 & cpa, Vector2 & cpb)
  {
    cpa = PointTraitsN<A, 2>::getPosition(a);
    cpb = PointTraitsN<B, 2>::getPosition(b);
    return MetricL2::monotoneApproxDistance(cpa, cpb);
  }
};

template <typename A, typename B>
struct /* THEA_API */ MetricL2Impl<A, B, typename boost::enable_if_c< IsPointN<A, 3>::value && IsPointN<B, 3>::value >::type>
{
  static double distance(A const & a, B const & b)
  { return (PointTraitsN<A, 3>::getPosition(a) - PointTraitsN<B, 3>::getPosition(b)).length(); }

  static double monotoneApproxDistance(A const & a, B const & b)
  { return (PointTraitsN<A, 3>::getPosition(a) - PointTraitsN<B, 3>::getPosition(b)).squaredLength(); }

  static double closestPoints(A const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    cpa = PointTraitsN<A, 3>::getPosition(a);
    cpb = PointTraitsN<B, 3>::getPosition(b);
    return MetricL2::monotoneApproxDistance(cpa, cpb);
  }
};

template <long N, typename T>
struct /* THEA_API */ MetricL2Impl< VectorN<N, T>, VectorN<N, T>, typename boost::enable_if_c< N != 2 && N != 3 >::type >
{
  static double distance(VectorN<N, T> const & a, VectorN<N, T> const & b) { return (a - b).length(); }
  static double monotoneApproxDistance(VectorN<N, T> const & a, VectorN<N, T> const & b) { return (a - b).squaredLength(); }

  static double closestPoints(VectorN<N, T> const & a, VectorN<N, T> const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  { cpa = a; cpb = b; return monotoneApproxDistance(a, b); }
};

template <typename B>
struct /* THEA_API */ MetricL2Impl<Line3, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  static double distance(Line3 const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(Line3 const & a, B const & b)
  {
    Vector3 p = PointTraitsN<B, 3>::getPosition(b);
    return (a.closestPoint(p) - p).squaredLength();
  }

  static double closestPoints(Line3 const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    cpb = PointTraitsN<B, 3>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredLength();
  }
};

template <typename A>
struct /* THEA_API */ MetricL2Impl<A, Line3, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  static double distance(A const & a, Line3 const & b)
  { return MetricL2Impl<Line3, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, Line3 const & b)
  { return MetricL2Impl<Line3, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, Line3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Line3, A>::closestPoints(b, a, cpb, cpa); }
};

template <typename B>
struct /* THEA_API */ MetricL2Impl<Ray3, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  static double distance(Ray3 const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(Ray3 const & a, B const & b)
  { return a.squaredDistance(PointTraitsN<B, 3>::getPosition(b)); }

  static double closestPoints(Ray3 const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    cpb = PointTraitsN<B, 3>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredLength();
  }
};

template <typename A>
struct /* THEA_API */ MetricL2Impl<A, Ray3, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  static double distance(A const & a, Ray3 const & b)
  { return MetricL2Impl<Ray3, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, Ray3 const & b)
  { return MetricL2Impl<Ray3, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, Ray3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Ray3, A>::closestPoints(b, a, cpb, cpa); }
};

template <typename B>
struct /* THEA_API */ MetricL2Impl<LineSegment3, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  static double distance(LineSegment3 const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(LineSegment3 const & a, B const & b)
  {
    Vector3 p = PointTraitsN<B, 3>::getPosition(b);
    return (a.closestPoint(p) - p).squaredLength();
  }

  static double closestPoints(LineSegment3 const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    cpb = PointTraitsN<B, 3>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredLength();
  }
};

template <typename A>
struct /* THEA_API */ MetricL2Impl<A, LineSegment3, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  static double distance(A const & a, LineSegment3 const & b)
  { return MetricL2Impl<LineSegment3, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, LineSegment3 const & b)
  { return MetricL2Impl<LineSegment3, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, LineSegment3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<LineSegment3, A>::closestPoints(b, a, cpb, cpa); }
};

template <typename B>
struct /* THEA_API */ MetricL2Impl<AxisAlignedBox3, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  static double distance(AxisAlignedBox3 const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(AxisAlignedBox3 const & a, B const & b)
  { return a.squaredDistance(PointTraitsN<B, 3>::getPosition(b)); }

  static double closestPoints(AxisAlignedBox3 const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <typename A>
struct /* THEA_API */ MetricL2Impl<A, AxisAlignedBox3, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  static double distance(A const & a, AxisAlignedBox3 const & b)
  { return MetricL2Impl<AxisAlignedBox3, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, AxisAlignedBox3 const & b)
  { return MetricL2Impl<AxisAlignedBox3, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, AxisAlignedBox3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<AxisAlignedBox3, A>::closestPoints(b, a, cpb, cpa); }
};

template <>
struct /* THEA_API */ MetricL2Impl<AxisAlignedBox3, AxisAlignedBox3>
{
  static double distance(AxisAlignedBox3 const & a, AxisAlignedBox3 const & b) { return a.distance(b); }
  static double monotoneApproxDistance(AxisAlignedBox3 const & a, AxisAlignedBox3 const & b) { return a.squaredDistance(b); }

  static double closestPoints(AxisAlignedBox3 const & a, AxisAlignedBox3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <long N, typename T>
struct /* THEA_API */ MetricL2Impl< AxisAlignedBoxN<N, T>, VectorN<N, T>,
                                    typename boost::enable_if_c< N != 2 && N != 3 >::type >
{
  static double distance(AxisAlignedBoxN<N, T> const & a, VectorN<N, T> const & b) { return a.distance(b); }

  static double monotoneApproxDistance(AxisAlignedBoxN<N, T> const & a, VectorN<N, T> const & b)
  { return a.squaredDistance(b); }

  static double closestPoints(AxisAlignedBoxN<N, T> const & a, VectorN<N, T> const & b, VectorN<N, T> & cpa,
                              VectorN<N, T> & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <long N, typename T>
struct /* THEA_API */ MetricL2Impl< VectorN<N, T>, AxisAlignedBoxN<N, T>,
                                    typename boost::enable_if_c< N != 2 && N != 3 >::type >
{
  static double distance(VectorN<N, T> const & a, AxisAlignedBoxN<N, T> const & b)
  { return MetricL2Impl< AxisAlignedBoxN<N, T>, VectorN<N, T> >::distance(b, a); }

  static double monotoneApproxDistance(VectorN<N, T> const & a, AxisAlignedBoxN<N, T> const & b)
  { return MetricL2Impl< AxisAlignedBoxN<N, T>, VectorN<N, T> >::monotoneApproxDistance(b, a); }

  static double closestPoints(VectorN<N, T> const & a, AxisAlignedBoxN<N, T> const & b, VectorN<N, T> & cpa,
                              VectorN<N, T> & cpb)
  { return MetricL2Impl< AxisAlignedBoxN<N, T>, VectorN<N, T> >::closestPoints(b, a, cpb, cpa); }
};

template <long N, typename T>
struct /* THEA_API */ MetricL2Impl< AxisAlignedBoxN<N, T>, AxisAlignedBoxN<N, T> >
{
  static double distance(AxisAlignedBoxN<N, T> const & a, AxisAlignedBoxN<N, T> const & b) { return a.distance(b); }

  static double monotoneApproxDistance(AxisAlignedBoxN<N, T> const & a, AxisAlignedBoxN<N, T> const & b)
  { return a.squaredDistance(b); }

  static double closestPoints(AxisAlignedBoxN<N, T> const & a, AxisAlignedBoxN<N, T> const & b, VectorN<N, T> & cpa,
                              VectorN<N, T> & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <typename B>
struct /* THEA_API */ MetricL2Impl<Ball3, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  static double distance(Ball3 const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(Ball3 const & a, B const & b)
  { return a.squaredDistance(PointTraitsN<B, 3>::getPosition(b)); }

  static double closestPoints(Ball3 const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <typename A>
struct /* THEA_API */ MetricL2Impl<A, Ball3, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  static double distance(A const & a, Ball3 const & b)
  { return MetricL2Impl<Ball3, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, Ball3 const & b)
  { return MetricL2Impl<Ball3, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, Ball3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Ball3, A>::closestPoints(b, a, cpb, cpa); }
};

template <>
struct /* THEA_API */ MetricL2Impl<Ball3, Ball3>
{
  static double distance(Ball3 const & a, Ball3 const & b) { return a.distance(b); }
  static double monotoneApproxDistance(Ball3 const & a, Ball3 const & b) { return a.squaredDistance(b); }

  static double closestPoints(Ball3 const & a, Ball3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <>
struct /* THEA_API */ MetricL2Impl<Ball3, AxisAlignedBox3>
{
  static double distance(Ball3 const & a, AxisAlignedBox3 const & b) { return a.distance(b); }
  static double monotoneApproxDistance(Ball3 const & a, AxisAlignedBox3 const & b) { return a.squaredDistance(b); }

  static double closestPoints(Ball3 const & a, AxisAlignedBox3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return monotoneApproxDistance(a, b); /* TODO */ }
};

template <>
struct /* THEA_API */ MetricL2Impl<AxisAlignedBox3, Ball3>
{
  static double distance(AxisAlignedBox3 const & a, Ball3 const & b)
  { return MetricL2Impl<Ball3, AxisAlignedBox3>::distance(b, a); }

  static double monotoneApproxDistance(AxisAlignedBox3 const & a, Ball3 const & b)
  { return MetricL2Impl<Ball3, AxisAlignedBox3>::monotoneApproxDistance(b, a); }

  static double closestPoints(AxisAlignedBox3 const & a, Ball3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Ball3, AxisAlignedBox3>::closestPoints(b, a, cpb, cpa); }
};

template <typename VertexTripleT, typename B>
struct /* THEA_API */ MetricL2Impl<Triangle3<VertexTripleT>, B, typename boost::enable_if< IsPointN<B, 3> >::type>
{
  typedef Triangle3<VertexTripleT> Triangle;

  static double distance(Triangle const & a, B const & b) { return a.distance(PointTraitsN<B, 3>::getPosition(b)); }

  static double monotoneApproxDistance(Triangle const & a, B const & b)
  { return a.squaredDistance(PointTraitsN<B, 3>::getPosition(b)); }

  static double closestPoints(Triangle const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    cpb = PointTraitsN<B, 3>::getPosition(b);
    cpa = a.closestPoint(cpb);
    return (cpa - cpb).squaredLength();
  }
};

template <typename A, typename VertexTripleT>
struct /* THEA_API */ MetricL2Impl<A, Triangle3<VertexTripleT>, typename boost::enable_if< IsPointN<A, 3> >::type>
{
  typedef Triangle3<VertexTripleT> Triangle;

  static double distance(A const & a, Triangle const & b)
  { return MetricL2Impl<Triangle, A>::distance(b, a); }

  static double monotoneApproxDistance(A const & a, Triangle const & b)
  { return MetricL2Impl<Triangle, A>::monotoneApproxDistance(b, a); }

  static double closestPoints(A const & a, Triangle const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Triangle, A>::closestPoints(b, a, cpb, cpa); }
};

template <typename VertexTripleT1, typename VertexTripleT2>
struct /* THEA_API */ MetricL2Impl< Triangle3<VertexTripleT1>, Triangle3<VertexTripleT2> >
{
  typedef Triangle3<VertexTripleT1> Triangle3_1;
  typedef Triangle3<VertexTripleT2> Triangle3_2;

  static double distance(Triangle3_1 const & a, Triangle3_2 const & b) { return a.distance(b); }
  static double monotoneApproxDistance(Triangle3_1 const & a, Triangle3_2 const & b) { return a.squaredDistance(b); }

  static double closestPoints(Triangle3_1 const & a, Triangle3_2 const & b, Vector3 & cpa, Vector3 & cpb)
  { return a.closestPoints(b, cpa, cpb); }
};

template <typename VertexTripleT>
struct /* THEA_API */ MetricL2Impl< Triangle3<VertexTripleT>, Ball3 >
{
  typedef Triangle3<VertexTripleT> Triangle;

  static double distance(Triangle const & a, Ball3 const & b) { return a.distance(b); }
  static double monotoneApproxDistance(Triangle const & a, Ball3 const & b) { return a.squaredDistance(b); }

  static double closestPoints(Triangle const & a, Ball3 const & b, Vector3 & cpa, Vector3 & cpb)
  { return a.closestPoints(b, cpa, cpb); }
};

template <typename VertexTripleT>
struct /* THEA_API */ MetricL2Impl< Ball3, Triangle3<VertexTripleT> >
{
  typedef Triangle3<VertexTripleT> Triangle;

  static double distance(Ball3 const & a, Triangle const & b)
  { return MetricL2Impl<Triangle, Ball3>::distance(b, a); }

  static double monotoneApproxDistance(Ball3 const & a, Triangle const & b)
  { return MetricL2Impl<Triangle, Ball3>::monotoneApproxDistance(b, a); }

  static double closestPoints(Ball3 const & a, Triangle const & b, Vector3 & cpa, Vector3 & cpb)
  { return MetricL2Impl<Triangle, Ball3>::closestPoints(b, a, cpb, cpa); }
};

} // namespace Algorithms
} // namespace Thea

#include "MetricL2_Transform.hpp"

#endif
