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

#ifndef __Thea_Algorithms_Transformer_hpp__
#define __Thea_Algorithms_Transformer_hpp__

#include "../Common.hpp"
#include "../RigidTransform3.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Box3.hpp"
#include "../Ball3.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

/**
 * Helper class for Transformer. Specializations of this class actually transform objects. This is required because C++ does
 * unexpected things with specialized and overloaded function templates (see http://www.gotw.ca/publications/mill17.htm).
 *
 * Every specialization of this class must suitably define the <code>Result</code> subtype, describing the result of the
 * transformation.
 *
 * @see Transformer
 */
template <typename ObjT, typename TransT, typename Enable = void>
struct /* THEA_API */ TransformerImpl
{
  typedef char Result;
  static Result transform(ObjT const & obj, TransT const & tr);

}; // struct TransformerImpl

/**
 * Apply transformations to objects.
 *
 * To add support for custom object and transformation types, add specializations (full or partial) of the helper class
 * TransformerImpl. This class is required because C++ does unexpected things with specialized and overloaded function templates
 * (see http://www.gotw.ca/publications/mill17.htm). Do <b>not</b> try specializing Transformer::transform().
 */
class THEA_API Transformer
{
  public:
    /** Apply a transformation to an object. */
    template <typename ObjT, typename TransT>
    static typename TransformerImpl<ObjT, TransT>::Result transform(ObjT const & obj, TransT const & tr)
    {
      return TransformerImpl<ObjT, TransT>::transform(obj, tr);
    }

}; // class Transformer

// Support for pointer types
template <typename ObjT, typename TransT>
struct /* THEA_API */ TransformerImpl<ObjT *, TransT *>
{
  typedef typename TransformerImpl<ObjT, TransT>::Result Result;
  static Result transform(ObjT const * obj, TransT const * tr) { return Transformer::transform(*obj, *tr); }
};

template <typename ObjT, typename TransT>
struct /* THEA_API */ TransformerImpl<ObjT, TransT *>
{
  typedef typename TransformerImpl<ObjT, TransT>::Result Result;
  static Result transform(ObjT const & obj, TransT const * tr) { return Transformer::transform(obj, *tr); }
};

template <typename ObjT, typename TransT>
struct /* THEA_API */ TransformerImpl<ObjT *, TransT>
{
  typedef typename TransformerImpl<ObjT, TransT>::Result Result;
  static Result transform(ObjT const * obj, TransT const & tr) { return Transformer::transform(*obj, tr); }
};

// Default specializations
template <typename ObjT, typename TransT>
struct TransformerImpl<ObjT, TransT, typename boost::enable_if< IsPointN<ObjT, 2> >::type>
{
  typedef Vector2 Result;
  static Result transform(ObjT const & obj, TransT const & tr) { return tr * PointTraitsN<ObjT, 2>::getPosition(obj); }
};

template <typename ObjT, typename TransT>
struct TransformerImpl<ObjT, TransT, typename boost::enable_if< IsPointN<ObjT, 3> >::type>
{
  typedef Vector3 Result;
  static Result transform(ObjT const & obj, TransT const & tr) { return tr * PointTraitsN<ObjT, 3>::getPosition(obj); }
};

template <typename TransT>
struct TransformerImpl<Vector4, TransT>
{
  typedef Vector4 Result;
  static Result transform(Vector4 const & v, TransT const & tr) { return tr * v; }
};

template <>
struct TransformerImpl<AxisAlignedBox3, RigidTransform3>
{
  typedef Box3 Result;
  static Result transform(AxisAlignedBox3 const & aabb, RigidTransform3 const & tr)
  { return Box3(aabb, CoordinateFrame3(tr)); }
};

template <>
struct TransformerImpl<Box3, RigidTransform3>
{
  typedef Box3 Result;
  static Result transform(Box3 const & box, RigidTransform3 const & tr)
  {
    return Box3(box.getLocalAAB(), CoordinateFrame3(tr) * box.getLocalFrame());
  }
};

template <>
struct TransformerImpl<Ball3, RigidTransform3>
{
  typedef Ball3 Result;
  static Ball3 transform(Ball3 const & ball, RigidTransform3 const & tr)
  { return Ball3(tr * ball.getCenter(), ball.getRadius()); }
};

template <typename VertexTripleT, typename TransT>
struct TransformerImpl< Triangle3<VertexTripleT>, TransT >
{
  typedef LocalTriangle3 Result;
  static Result transform(Triangle3<VertexTripleT> const & tri, TransT const & tr)
  {
    return LocalTriangle3(Transformer::transform(tri.getVertex(0), tr),
                          Transformer::transform(tri.getVertex(1), tr),
                          Transformer::transform(tri.getVertex(2), tr));
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
