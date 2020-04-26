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

#ifndef __Thea_Algorithms_IntersectionTester_hpp__
#define __Thea_Algorithms_IntersectionTester_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../BoxN.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include "TransformedObject.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

/**
 * Helper class for IntersectionTester. Specializations of this class actually test for intersection. This is required because
 * C++ does unexpected things with specialized and overloaded function templates (see
 * http://www.gotw.ca/publications/mill17.htm).
 *
 * @see IntersectionTester
 */
template <typename A, typename B, int N, typename T, typename Enable = void>
struct /* THEA_API */ IntersectionTesterImpl
{
  /** Check if two objects intersect. */
  static bool intersects(A const & a, B const & b) { return a.intersects(b); }

}; // class IntersectionTesterImpl

/**
 * %Intersection queries on objects.
 *
 * To add support for intersection queries between custom types, add specializations (full or partial) of the helper class
 * IntersectionTesterImpl. This class is required because C++ does unexpected things with specialized and overloaded function
 * templates (see http://www.gotw.ca/publications/mill17.htm). Note that to commutatively support queries between distinct types
 * A and B, you must specialize both IntersectionTesterImpl<A, B> and IntersectionTesterImpl<B, A>. Do <b>not</b> try
 * specializing IntersectionTester::intersects<N, T>().
 */
class THEA_API IntersectionTester
{
  public:
    /**
     * Check if two objects intersect, via the helper class IntersectionTesterImpl. Add specializations of
     * IntersectionTesterImpl as required for specific types of objects. Note that to commutatively support intersection queries
     * on distinct types A and B, you must explicitly specialize both IntersectionTesterImpl<A, B> and
     * IntersectionTesterImpl<B, A>.
     */
    template <int N, typename T, typename A, typename B> static bool intersects(A const & a, B const & b)
    { return IntersectionTesterImpl<A, B, N, T>::intersects(a, b); }

}; // class IntersectionTester

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ IntersectionTesterImpl<A *, B *, N, T>
{
  static bool intersects(A const * a, B const * b) { return IntersectionTester::intersects<N, T>(*a, *b); }
};

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ IntersectionTesterImpl<A, B *, N, T>
{
  static bool intersects(A const & a, B const * b) { return IntersectionTester::intersects<N, T>(a, *b); }
};

template <typename A, typename B, int N, typename T>
struct /* THEA_API */ IntersectionTesterImpl<A *, B, N, T>
{
  static bool intersects(A const * a, B const & b) { return IntersectionTester::intersects<N, T>(*a, b); }
};

//=============================================================================================================================
// Default specializations
//=============================================================================================================================

namespace IntersectionTesterInternal {

template <typename T> struct TransformedObjectCheck { static bool const value = false; };

template <typename ObjT, typename TransT> struct TransformedObjectCheck< TransformedObject<ObjT, TransT> >
{
  static bool const value = true;
};

} // namespace IntersectionTesterInternal

// Only the second object is a point
template <typename A, typename B, int N, typename T>
struct IntersectionTesterImpl<A, B, N, T,
                              typename std::enable_if< !IsPointN<A, N>::value
                                                    && !IntersectionTesterInternal::TransformedObjectCheck<A>::value
                                                    && !std::is_pointer<A>::value
                                                    && IsNonReferencedPointN<B, N>::value >::type>
{
  static bool intersects(A const & a, B const & b) { return a.intersects(PointTraitsN<B, N, T>::getPosition(b)); }
};

// Only the first object is a point
template <typename A, typename B, int N, typename T>
struct IntersectionTesterImpl<A, B, N, T,
                              typename std::enable_if< IsNonReferencedPointN<A, N>::value
                                                    && !IsPointN<B, N>::value
                                                    && !IntersectionTesterInternal::TransformedObjectCheck<B>::value
                                                    && !std::is_pointer<B>::value >::type>
{
  static bool intersects(A const & a, B const & b) { return b.intersects(PointTraitsN<A, N, T>::getPosition(a)); }
};

//=============================================================================================================================
// Specializations where the intersects() member function is not defined in both classes
//=============================================================================================================================

template <typename VertexTripleT, typename T>
struct IntersectionTesterImpl<AxisAlignedBoxN<3, T>, Triangle3<VertexTripleT>, 3, T>
{
  static bool intersects(AxisAlignedBoxN<3, T> const & a, Triangle3<VertexTripleT> const & b) { return b.intersects(a); }
};

template <typename VertexTripleT, typename T>
struct IntersectionTesterImpl<BallN<3, T>, Triangle3<VertexTripleT>, 3, T>
{
  static bool intersects(BallN<3, T> const & a, Triangle3<VertexTripleT> const & b) { return b.intersects(a); }
};

template <typename VertexTripleT, typename T>
struct IntersectionTesterImpl<BoxN<3, T>, Triangle3<VertexTripleT>, 3, T>
{
  static bool intersects(BoxN<3, T> const & a, Triangle3<VertexTripleT> const & b) { return b.intersects(a); }
};

template <int N, typename T>
struct IntersectionTesterImpl<AxisAlignedBoxN<N, T>, BoxN<N, T>, N, T>
{
  static bool intersects(AxisAlignedBoxN<N, T> const & a, BoxN<N, T> const & b) { return b.intersects(a); }
};

} // namespace Algorithms
} // namespace Thea

#include "IntersectionTester_Transform.hpp"

#endif
