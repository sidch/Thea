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

#ifndef __Thea_Algorithms_BoundedObjectTraits3_hpp__
#define __Thea_Algorithms_BoundedObjectTraits3_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Ball3.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

/** Traits class for a bounded object in 3-space. */
template <typename T, typename Enable = void>
class /* THEA_API */ BoundedObjectTraits3
{
  public:
    /** Get a bounding range for an object. */
    template <typename RangeT> static void getBounds(T const & t, RangeT & bounds) { bounds = t.getBounds(); }

    /** Get the center of the object. */
    static Vector3 getCenter(T const & t) { return t.getBounds().getCenter(); }

    /** Get the maximum position of the object along a particular coordinate axis. */
    static Real getHigh(T const & t, int coord) { return t.getBounds().getHigh()[coord]; }

    /** Get the minimum position of the object along a particular coordinate axis. */
    static Real getLow(T const & t, int coord) { return t.getBounds().getLow()[coord]; }

}; // class BoundedObjectTraits3

// Specialization for pointer types
template <typename T>
struct /* THEA_API */ BoundedObjectTraits3<T *>
{
  template <typename RangeT> static void getBounds(T const * t, RangeT & bounds)
  {
    debugAssertM(t, "BoundedObjectTraits3: Can't get bounds of null object");
    BoundedObjectTraits3<T>::getBounds(*t, bounds);
  }

  static Vector3 getCenter(T const * t)
  {
    debugAssertM(t, "BoundedObjectTraits3: Can't get center of null object");
    return BoundedObjectTraits3<T>::getCenter(*t);
  }

  static Real getHigh(T const * t, int coord)
  {
    debugAssertM(t, "BoundedObjectTraits3: Can't get bounds of null object");
    return BoundedObjectTraits3<T>::getHigh(*t, coord);
  }

  static Real getLow(T const * t, int coord)
  {
    debugAssertM(t, "BoundedObjectTraits3: Can't get bounds of null object");
    return BoundedObjectTraits3<T>::getLow(*t, coord);
  }
};

// Specialization for 3D points
template <typename T>
struct /* THEA_API */ BoundedObjectTraits3<T, typename boost::enable_if< IsPointN<T, 3> >::type>
{
  static void getBounds(T const & t, AxisAlignedBox3 & bounds)
  {
    Vector3 pos = PointTraitsN<T, 3>::getPosition(t);
    bounds.set(pos, pos);
  }

  static void getBounds(T const & t, Ball3 & bounds) { bounds = Ball3(PointTraitsN<T, 3>::getPosition(t), 0); }
  static Vector3 getCenter(T const & t) { return PointTraitsN<T, 3>::getPosition(t); }
  static Real getHigh(T const & t, int coord) { return PointTraitsN<T, 3>::getPosition(t)[coord]; }
  static Real getLow(T const & t, int coord)  { return PointTraitsN<T, 3>::getPosition(t)[coord];  }
};

// Specialization for AxisAlignedBox3
template <>
struct /* THEA_API */ BoundedObjectTraits3<AxisAlignedBox3>
{
  static void getBounds(AxisAlignedBox3 const & t, AxisAlignedBox3 & bounds) { bounds = t; }

  static void getBounds(AxisAlignedBox3 const & t, Ball3 & bounds)
  { bounds = Ball3(t.getCenter(), 0.5f * t.getExtent().length()); }

  static Vector3 getCenter(AxisAlignedBox3 const & t) { return t.getCenter(); }
  static Real getHigh(AxisAlignedBox3 const & t, int coord) { return t.getHigh()[coord]; }
  static Real getLow(AxisAlignedBox3 const & t, int coord)  { return t.getLow()[coord];  }
};

// Specialization for Ball3
template <>
struct /* THEA_API */ BoundedObjectTraits3<Ball3>
{
  static void getBounds(Ball3 const & t, AxisAlignedBox3 & bounds) { bounds = t.getBounds(); }
  static void getBounds(Ball3 const & t, Ball3 & bounds)           { bounds = t; }

  static Vector3 getCenter(Ball3 const & t) { return t.getCenter(); }
  static Real getHigh(Ball3 const & t, int coord) { return t.getCenter()[coord] + t.getRadius(); }
  static Real getLow(Ball3 const & t, int coord)  { return t.getCenter()[coord] - t.getRadius(); }
};

// Specialization for Triangle3
template <typename VertexTripleT>
struct /* THEA_API */ BoundedObjectTraits3< Triangle3<VertexTripleT> >
{
  typedef Triangle3<VertexTripleT> Triangle;

  static void getBounds(Triangle const & t, AxisAlignedBox3 & bounds) { bounds = t.getBounds(); }

  // TODO: Make this tighter
  static void getBounds(Triangle const & t, Ball3 & bounds)
  { BoundedObjectTraits3<AxisAlignedBox3>::getBounds(t.getBounds(), bounds); }

  static Vector3 getCenter(Triangle const & t) { return (t.getVertex(0) + t.getVertex(1) + t.getVertex(2)) / 3.0f; }

  static Real getHigh(Triangle const & t, int coord)
  { return std::max(std::max(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }

  static Real getLow(Triangle const & t, int coord)
  { return std::min(std::min(t.getVertex(0)[coord], t.getVertex(1)[coord]), t.getVertex(2)[coord]); }
};

} // namespace Algorithms
} // namespace Thea

#include "BoundedObjectTraits3_Transform.hpp"

#endif
