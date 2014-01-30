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

#ifndef __Thea_Algorithms_IntersectionTester_hpp__
#define __Thea_Algorithms_IntersectionTester_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BallN.hpp"
#include "../BoxN.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include "TransformedObject.hpp"
#include <boost/type_traits/is_pointer.hpp>
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

/**
 * Helper class for IntersectionTester. Specializations of this class actually test for intersection. This is required because
 * C++ does unexpected things with specialized and overloaded function templates (see
 * http://www.gotw.ca/publications/mill17.htm).
 *
 * @see IntersectionTester
 */
template <typename A, typename B, long N, typename T, typename Enable = void>
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
    template <long N, typename T, typename A, typename B> static bool intersects(A const & a, B const & b)
    { return IntersectionTesterImpl<A, B, N, T>::intersects(a, b); }

}; // class IntersectionTester

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename A, typename B, long N, typename T>
struct /* THEA_API */ IntersectionTesterImpl<A *, B *, N, T>
{
  static bool intersects(A const * a, B const * b) { return IntersectionTester::intersects<N, T>(*a, *b); }
};

template <typename A, typename B, long N, typename T>
struct /* THEA_API */ IntersectionTesterImpl<A, B *, N, T>
{
  static bool intersects(A const & a, B const * b) { return IntersectionTester::intersects<N, T>(a, *b); }
};

template <typename A, typename B, long N, typename T>
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
template <typename A, typename B, long N, typename T>
struct IntersectionTesterImpl<A, B, N, T,
                              typename boost::enable_if_c< !IsPointN<A, N>::value
                                                        && !IntersectionTesterInternal::TransformedObjectCheck<A>::value
                                                        && !boost::is_pointer<A>::value
                                                        && IsPointN<B, N>::value >::type>
{
  static bool intersects(A const & a, B const & b) { return a.intersects(PointTraitsN<B, N, T>::getPosition(b)); }
};

// Only the first object is a point
template <typename A, typename B, long N, typename T>
struct IntersectionTesterImpl<A, B, N, T,
                              typename boost::enable_if_c< IsPointN<A, N>::value
                                                        && !IsPointN<B, N>::value
                                                        && !IntersectionTesterInternal::TransformedObjectCheck<B>::value
                                                        && !boost::is_pointer<B>::value >::type>
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

template <long N, typename T>
struct IntersectionTesterImpl<AxisAlignedBoxN<N, T>, BoxN<N, T>, N, T>
{
  static bool intersects(AxisAlignedBoxN<N, T> const & a, BoxN<N, T> const & b) { return b.intersects(a); }
};

} // namespace Algorithms
} // namespace Thea

#include "IntersectionTester_Transform.hpp"

#endif
