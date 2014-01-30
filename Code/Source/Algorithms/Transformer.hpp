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

#ifndef __Thea_Algorithms_Transformer_hpp__
#define __Thea_Algorithms_Transformer_hpp__

#include "../Common.hpp"
#include "../RigidTransformN.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BoxN.hpp"
#include "../BallN.hpp"
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
template <typename ObjT, typename TransT, long N, typename ScalarT, typename Enable = void>
struct /* THEA_API */ TransformerImpl
{
  typedef ObjT Result;
  static Result transform(ObjT const & obj, TransT const & tr) { return tr * obj; }

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
    template <long N, typename ScalarT, typename ObjT, typename TransT>
    static typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result transform(ObjT const & obj, TransT const & tr)
    {
      return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(obj, tr);
    }

}; // class Transformer

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename ObjT, typename TransT, long N, typename ScalarT>
struct /* THEA_API */ TransformerImpl<ObjT *, TransT *, N, ScalarT>
{
  typedef typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result Result;

  static Result transform(ObjT const * obj, TransT const * tr)
  {
    return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(*obj, *tr);
  }
};

template <typename ObjT, typename TransT, long N, typename ScalarT>
struct /* THEA_API */ TransformerImpl<ObjT, TransT *, N, ScalarT>
{
  typedef typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result Result;

  static Result transform(ObjT const & obj, TransT const * tr)
  {
    return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(obj, *tr);
  }
};

template <typename ObjT, typename TransT, long N, typename ScalarT>
struct /* THEA_API */ TransformerImpl<ObjT *, TransT, N, ScalarT>
{
  typedef typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result Result;

  static Result transform(ObjT const * obj, TransT const & tr)
  {
    return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(*obj, tr);
  }
};

//=============================================================================================================================
// Specializations
//=============================================================================================================================

template <typename ObjT, typename TransT, long N, typename ScalarT>
struct TransformerImpl<ObjT, TransT, N, ScalarT, typename boost::enable_if< IsPointN<ObjT, N> >::type>
{
  typedef VectorN<N, ScalarT> Result;
  static Result transform(ObjT const & obj, TransT const & tr) { return tr * PointTraitsN<ObjT, N, ScalarT>::getPosition(obj); }
};

template <long N, typename ScalarT>
struct TransformerImpl< AxisAlignedBoxN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BoxN<N, ScalarT> Result;
  static Result transform(AxisAlignedBoxN<N, ScalarT> const & aabb, RigidTransformN<N, ScalarT> const & tr)
  {
    return Result(aabb, CoordinateFrameN<N, ScalarT>(tr));
  }
};

template <long N, typename ScalarT>
struct TransformerImpl< BoxN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BoxN<N, ScalarT> Result;
  static Result transform(BoxN<N, ScalarT> const & box, RigidTransformN<N, ScalarT> const & tr)
  {
    return Result(box.getLocalAAB(), CoordinateFrameN<N, ScalarT>(tr) * box.getLocalFrame());
  }
};

template <long N, typename ScalarT>
struct TransformerImpl< BallN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BallN<N, ScalarT> Result;
  static Result transform(BallN<N, ScalarT> const & ball, RigidTransformN<N, ScalarT> const & tr)
  { return BallN<N, ScalarT>(tr * ball.getCenter(), ball.getRadius()); }
};

template <typename VertexTripleT, typename TransT, typename ScalarT>
struct THEA_API TransformerImpl< Triangle3<VertexTripleT>, TransT, 3, ScalarT >
{
  typedef LocalTriangle3 Result;
  static Result transform(Triangle3<VertexTripleT> const & tri, TransT const & tr)
  {
    return LocalTriangle3(Transformer::transform<3, ScalarT>(tri.getVertex(0), tr),
                          Transformer::transform<3, ScalarT>(tri.getVertex(1), tr),
                          Transformer::transform<3, ScalarT>(tri.getVertex(2), tr));
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
