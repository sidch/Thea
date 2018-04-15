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

#ifndef __Thea_Algorithms_BoundedTraitsN_hpp__
#define __Thea_Algorithms_BoundedTraitsN_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../BoxN.hpp"
#include "../LineSegmentN.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_pointer.hpp>

namespace Thea {
namespace Algorithms {

//=============================================================================================================================
// Flag that should be set to true for bounded objects
//=============================================================================================================================

/**
 * Has boolean member <code>value = true</code> if <code>T</code> is a geometrically bounded object in n-D space, else
 * <code>value = false</code>. Unless you specialize the class to set the value to true, it is false by default.
 *
 * @see BoundedTraitsN
 */
template <typename T, long N, typename Enable = void>
class IsBoundedN
{
  public:
    static bool const value = false;
};

// Partial specialization for const and pointer types
template <typename T, long N> class IsBoundedN<T const, N> { public: static bool const value = IsBoundedN<T, N>::value; };
template <typename T, long N> class IsBoundedN<T *, N>     { public: static bool const value = IsBoundedN<T, N>::value; };

// Points are bounded
template <typename T, long N>
class IsBoundedN< T, N, typename boost::enable_if< IsRawPointN<T, N> >::type >
{ public: static bool const value = true; };

// ... as are many simple geometric objects
template <long N, typename S> class IsBoundedN< AxisAlignedBoxN<N, S>, N >         { public: static bool const value = true; };
template <long N, typename S> class IsBoundedN< BallN<N, S>, N >                   { public: static bool const value = true; };
template <long N, typename S> class IsBoundedN< BoxN<N, S>, N >                    { public: static bool const value = true; };
template <long N, typename S> class IsBoundedN< LineSegmentN<N, S>, N >            { public: static bool const value = true; };
template <typename VertexTripleT> class IsBoundedN< Triangle3<VertexTripleT>, 3 >  { public: static bool const value = true; };

/** Same as IsBoundedN (no need to specialize it separately), except false for const or pointer types. */
template <typename T, long N>
class /* THEA_API */ IsRawBoundedN
{
  public:
    static bool const value = IsBoundedN<T, N>::value
                           && !boost::is_const<T>::value
                           && !boost::is_pointer<T>::value;
};

/** Same as IsBoundedN (no need to specialize it separately), except false for pointer types. */
template <typename T, long N>
class /* THEA_API */ IsNonReferencedBoundedN
{
  public:
    static bool const value = IsBoundedN<T, N>::value
                           && !boost::is_pointer<T>::value;
};

//=============================================================================================================================
// Traits class to get bounding volumes of bounded objects
//=============================================================================================================================

/**
 * Traits class for a bounded object in N-space. Default implementation assumes class T has a T::getBounds() function that
 * returns an N-dimensional axis-aligned bounding box. Different specializations exist for specific classes.
 */
template <typename T, long N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ BoundedTraitsN
{
  public:
    /** Get a bounding box for an object. */
    static void getBounds(T const & t, AxisAlignedBoxN<N, ScalarT> & bounds) { bounds = t.getBounds(); }

    /** Get a bounding sphere for an object. */
    static void getBounds(T const & t, BallN<N, ScalarT> & bounds)
    {
      // Can we make this tighter in general?
      BoundedTraitsN<AxisAlignedBoxN<N, ScalarT>, N, ScalarT>::getBounds(t.getBounds(), bounds);
    }

    /** Get the center of the object. */
    static VectorN<N, ScalarT> getCenter(T const & t) { return t.getBounds().getCenter(); }

    /** Get the maximum position of the object along a particular coordinate axis. */
    static ScalarT getHigh(T const & t, long coord) { return t.getBounds().getHigh()[coord]; }

    /** Get the minimum position of the object along a particular coordinate axis. */
    static ScalarT getLow(T const & t, long coord) { return t.getBounds().getLow()[coord]; }

}; // class BoundedTraitsN

// Specialization for const types
template <typename T, long N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T const, N, ScalarT>
{
  template <typename RangeT> static void getBounds(T const & t, RangeT & bounds)
  { BoundedTraitsN<T, N, ScalarT>::getBounds(t, bounds); }

  static VectorN<N, ScalarT> getCenter(T const & t)  { return BoundedTraitsN<T, N, ScalarT>::getCenter(t);      }
  static ScalarT getHigh(T const & t, long coord)    { return BoundedTraitsN<T, N, ScalarT>::getHigh(t, coord); }
  static ScalarT getLow(T const & t, long coord)     { return BoundedTraitsN<T, N, ScalarT>::getLow(t, coord);  }
};

// Specialization for pointer types
template <typename T, long N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T *, N, ScalarT>
{
  template <typename RangeT> static void getBounds(T const * t, RangeT & bounds)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    BoundedTraitsN<T, N, ScalarT>::getBounds(*t, bounds);
  }

  static VectorN<N, ScalarT> getCenter(T const * t)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get center of null object");
    return BoundedTraitsN<T, N, ScalarT>::getCenter(*t);
  }

  static ScalarT getHigh(T const * t, long coord)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    return BoundedTraitsN<T, N, ScalarT>::getHigh(*t, coord);
  }

  static ScalarT getLow(T const * t, long coord)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    return BoundedTraitsN<T, N, ScalarT>::getLow(*t, coord);
  }
};

// Specialization for points
template <typename T, long N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T, N, ScalarT, typename boost::enable_if< IsRawPointN<T, N> >::type>
{
  static void getBounds(T const & t, AxisAlignedBoxN<N, ScalarT> & bounds)
  {
    VectorN<N, ScalarT> pos = PointTraitsN<T, N, ScalarT>::getPosition(t);
    bounds.set(pos, pos);
  }

  static void getBounds(T const & t, BallN<N, ScalarT> & bounds)
  {
    bounds = BallN<N, ScalarT>(PointTraitsN<T, N, ScalarT>::getPosition(t), 0);
  }

  static VectorN<N, ScalarT> getCenter(T const & t) { return PointTraitsN<T, N, ScalarT>::getPosition(t); }
  static ScalarT getHigh(T const & t, long coord) { return PointTraitsN<T, N, ScalarT>::getPosition(t)[coord]; }
  static ScalarT getLow(T const & t, long coord)  { return PointTraitsN<T, N, ScalarT>::getPosition(t)[coord];  }
};

// Specialization for AxisAlignedBoxN
template <long N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN< AxisAlignedBoxN<N, ScalarT>, N, ScalarT >
{
  static void getBounds(AxisAlignedBoxN<N, ScalarT> const & t, AxisAlignedBoxN<N, ScalarT> & bounds) { bounds = t; }

  static void getBounds(AxisAlignedBoxN<N, ScalarT> const & t, BallN<N, ScalarT> & bounds)
  { bounds = BallN<N, ScalarT>(t.getCenter(), 0.5f * t.getExtent().length()); }

  static VectorN<N, ScalarT> getCenter(AxisAlignedBoxN<N, ScalarT> const & t) { return t.getCenter(); }
  static ScalarT getHigh(AxisAlignedBoxN<N, ScalarT> const & t, long coord) { return t.getHigh()[coord]; }
  static ScalarT getLow(AxisAlignedBoxN<N, ScalarT> const & t, long coord)  { return t.getLow()[coord];  }
};

// Specialization for BallN
template <long N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN< BallN<N, ScalarT>, N, ScalarT >
{
  static void getBounds(BallN<N, ScalarT> const & t, AxisAlignedBoxN<N, ScalarT> & bounds) { bounds = t.getBounds(); }
  static void getBounds(BallN<N, ScalarT> const & t, BallN<N, ScalarT> & bounds)           { bounds = t; }

  static VectorN<N, ScalarT> getCenter(BallN<N, ScalarT> const & t) { return t.getCenter(); }
  static ScalarT getHigh(BallN<N, ScalarT> const & t, long coord) { return t.getCenter()[coord] + t.getRadius(); }
  static ScalarT getLow(BallN<N, ScalarT> const & t, long coord)  { return t.getCenter()[coord] - t.getRadius(); }
};

// Specialization for Triangle3
template <typename VertexTripleT, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN< Triangle3<VertexTripleT>, 3, ScalarT >
{
  typedef Triangle3<VertexTripleT> Triangle;

  static void getBounds(Triangle const & t, AxisAlignedBox3 & bounds) { bounds = t.getBounds(); }

  // TODO: Make this tighter
  static void getBounds(Triangle const & t, Ball3 & bounds)
  { BoundedTraitsN<AxisAlignedBox3, 3, ScalarT>::getBounds(t.getBounds(), bounds); }

  static VectorN<3, ScalarT> getCenter(Triangle const & t) { return (t.getVertex(0) + t.getVertex(1) + t.getVertex(2)) / 3; }

  static ScalarT getHigh(Triangle const & t, long coord)
  { return (ScalarT)std::max(std::max(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }

  static ScalarT getLow(Triangle const & t, long coord)
  { return (ScalarT)std::min(std::min(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }
};

} // namespace Algorithms
} // namespace Thea

#include "BoundedTraitsN_Transform.hpp"

#endif
