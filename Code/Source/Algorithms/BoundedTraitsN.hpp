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

#ifndef __Thea_Algorithms_BoundedTraitsN_hpp__
#define __Thea_Algorithms_BoundedTraitsN_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../BoxN.hpp"
#include "../LineSegmentN.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include <type_traits>

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
template <typename T, int N, typename Enable = void>
class IsBoundedN
{
  public:
    static bool const value = false;
};

// Partial specialization for const and pointer types
template <typename T, int N> class IsBoundedN<T const, N> { public: static bool const value = IsBoundedN<T, N>::value; };
template <typename T, int N> class IsBoundedN<T *, N>     { public: static bool const value = IsBoundedN<T, N>::value; };

// Points are bounded
template <typename T, int N>
class IsBoundedN< T, N, typename std::enable_if< IsRawPointN<T, N>::value >::type >
{ public: static bool const value = true; };

// ... as are many simple geometric objects
template <int N, typename S> class IsBoundedN< AxisAlignedBoxN<N, S>, N >         { public: static bool const value = true; };
template <int N, typename S> class IsBoundedN< BallN<N, S>, N >                   { public: static bool const value = true; };
template <int N, typename S> class IsBoundedN< BoxN<N, S>, N >                    { public: static bool const value = true; };
template <int N, typename S> class IsBoundedN< LineSegmentN<N, S>, N >            { public: static bool const value = true; };
template <typename VertexTripleT> class IsBoundedN< Triangle3<VertexTripleT>, 3 >  { public: static bool const value = true; };

/** Same as IsBoundedN (no need to specialize it separately), except false for const or pointer types. */
template <typename T, int N>
class /* THEA_API */ IsRawBoundedN
{
  public:
    static bool const value = IsBoundedN<T, N>::value
                           && !std::is_const<T>::value
                           && !std::is_pointer<T>::value;
};

/** Same as IsBoundedN (no need to specialize it separately), except false for pointer types. */
template <typename T, int N>
class /* THEA_API */ IsNonReferencedBoundedN
{
  public:
    static bool const value = IsBoundedN<T, N>::value
                           && !std::is_pointer<T>::value;
};

//=============================================================================================================================
// Traits class to get bounding volumes of bounded objects
//=============================================================================================================================

/**
 * Traits class for a bounded object in N-space. Default implementation assumes class T has a T::getBounds() function that
 * returns an N-dimensional axis-aligned bounding box. Different specializations exist for specific classes.
 */
template <typename T, int N, typename ScalarT = Real, typename Enable = void>
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
    static Vector<N, ScalarT> getCenter(T const & t) { return t.getBounds().getCenter(); }

    /** Get the maximum position of the object along a particular coordinate axis. */
    static ScalarT getHigh(T const & t, intx coord) { return t.getBounds().getHigh()[coord]; }

    /** Get the minimum position of the object along a particular coordinate axis. */
    static ScalarT getLow(T const & t, intx coord) { return t.getBounds().getLow()[coord]; }

}; // class BoundedTraitsN

// Specialization for const types
template <typename T, int N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T const, N, ScalarT>
{
  template <typename RangeT> static void getBounds(T const & t, RangeT & bounds)
  { BoundedTraitsN<T, N, ScalarT>::getBounds(t, bounds); }

  static Vector<N, ScalarT> getCenter(T const & t)  { return BoundedTraitsN<T, N, ScalarT>::getCenter(t);      }
  static ScalarT getHigh(T const & t, intx coord)    { return BoundedTraitsN<T, N, ScalarT>::getHigh(t, coord); }
  static ScalarT getLow(T const & t, intx coord)     { return BoundedTraitsN<T, N, ScalarT>::getLow(t, coord);  }
};

// Specialization for pointer types
template <typename T, int N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T *, N, ScalarT>
{
  template <typename RangeT> static void getBounds(T const * t, RangeT & bounds)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    BoundedTraitsN<T, N, ScalarT>::getBounds(*t, bounds);
  }

  static Vector<N, ScalarT> getCenter(T const * t)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get center of null object");
    return BoundedTraitsN<T, N, ScalarT>::getCenter(*t);
  }

  static ScalarT getHigh(T const * t, intx coord)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    return BoundedTraitsN<T, N, ScalarT>::getHigh(*t, coord);
  }

  static ScalarT getLow(T const * t, intx coord)
  {
    debugAssertM(t, "BoundedTraitsN: Can't get bounds of null object");
    return BoundedTraitsN<T, N, ScalarT>::getLow(*t, coord);
  }
};

// Specialization for points
template <typename T, int N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN<T, N, ScalarT, typename std::enable_if< IsRawPointN<T, N>::value >::type>
{
  static void getBounds(T const & t, AxisAlignedBoxN<N, ScalarT> & bounds)
  {
    Vector<N, ScalarT> pos = PointTraitsN<T, N, ScalarT>::getPosition(t);
    bounds.set(pos, pos);
  }

  static void getBounds(T const & t, BallN<N, ScalarT> & bounds)
  {
    bounds = BallN<N, ScalarT>(PointTraitsN<T, N, ScalarT>::getPosition(t), 0);
  }

  static Vector<N, ScalarT> getCenter(T const & t) { return PointTraitsN<T, N, ScalarT>::getPosition(t); }
  static ScalarT getHigh(T const & t, intx coord) { return PointTraitsN<T, N, ScalarT>::getPosition(t)[coord]; }
  static ScalarT getLow(T const & t, intx coord)  { return PointTraitsN<T, N, ScalarT>::getPosition(t)[coord];  }
};

// Specialization for AxisAlignedBoxN
template <int N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN< AxisAlignedBoxN<N, ScalarT>, N, ScalarT >
{
  static void getBounds(AxisAlignedBoxN<N, ScalarT> const & t, AxisAlignedBoxN<N, ScalarT> & bounds) { bounds = t; }

  static void getBounds(AxisAlignedBoxN<N, ScalarT> const & t, BallN<N, ScalarT> & bounds)
  { bounds = BallN<N, ScalarT>(t.getCenter(), 0.5f * t.getExtent().norm()); }

  static Vector<N, ScalarT> getCenter(AxisAlignedBoxN<N, ScalarT> const & t) { return t.getCenter(); }
  static ScalarT getHigh(AxisAlignedBoxN<N, ScalarT> const & t, intx coord) { return t.getHigh()[coord]; }
  static ScalarT getLow(AxisAlignedBoxN<N, ScalarT> const & t, intx coord)  { return t.getLow()[coord];  }
};

// Specialization for BallN
template <int N, typename ScalarT>
struct /* THEA_API */ BoundedTraitsN< BallN<N, ScalarT>, N, ScalarT >
{
  static void getBounds(BallN<N, ScalarT> const & t, AxisAlignedBoxN<N, ScalarT> & bounds) { bounds = t.getBounds(); }
  static void getBounds(BallN<N, ScalarT> const & t, BallN<N, ScalarT> & bounds)           { bounds = t; }

  static Vector<N, ScalarT> getCenter(BallN<N, ScalarT> const & t) { return t.getCenter(); }
  static ScalarT getHigh(BallN<N, ScalarT> const & t, intx coord) { return t.getCenter()[coord] + t.getRadius(); }
  static ScalarT getLow(BallN<N, ScalarT> const & t, intx coord)  { return t.getCenter()[coord] - t.getRadius(); }
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

  static Vector<3, ScalarT> getCenter(Triangle const & t) { return (t.getVertex(0) + t.getVertex(1) + t.getVertex(2)) / 3; }

  static ScalarT getHigh(Triangle const & t, intx coord)
  { return (ScalarT)std::max(std::max(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }

  static ScalarT getLow(Triangle const & t, intx coord)
  { return (ScalarT)std::min(std::min(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }
};

} // namespace Algorithms
} // namespace Thea

#include "BoundedTraitsN_Transform.hpp"

#endif
