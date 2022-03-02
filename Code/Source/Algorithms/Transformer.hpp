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

#ifndef __Thea_Algorithms_Transformer_hpp__
#define __Thea_Algorithms_Transformer_hpp__

#include "../Common.hpp"
#include "../Concept.hpp"
#include "../RigidTransformN.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../BoxN.hpp"
#include "../BallN.hpp"
#include "../LineSegmentN.hpp"
#include "../MatVec.hpp"
#include "../Triangle3.hpp"
#include "PointTraitsN.hpp"
#include <type_traits>

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
template <typename ObjT, typename TransT, int N, typename ScalarT, typename Enable = void>
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
    template <int N, typename ScalarT, typename ObjT, typename TransT>
    static typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result transform(ObjT const & obj, TransT const & tr)
    {
      return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(obj, tr);
    }

}; // class Transformer

//=============================================================================================================================
// Support for pointer types
//=============================================================================================================================

template <typename ObjT, typename TransT, int N, typename ScalarT>
struct /* THEA_API */ TransformerImpl<ObjT *, TransT *, N, ScalarT>
{
  typedef typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result Result;

  static Result transform(ObjT const * obj, TransT const * tr)
  {
    return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(*obj, *tr);
  }
};

template <typename ObjT, typename TransT, int N, typename ScalarT>
struct /* THEA_API */ TransformerImpl<ObjT, TransT *, N, ScalarT>
{
  typedef typename TransformerImpl<ObjT, TransT, N, ScalarT>::Result Result;

  static Result transform(ObjT const & obj, TransT const * tr)
  {
    return TransformerImpl<ObjT, TransT, N, ScalarT>::transform(obj, *tr);
  }
};

template <typename ObjT, typename TransT, int N, typename ScalarT>
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

namespace TransformerInternal {

template <typename T, int N, typename Enable = void>
struct IsHomMatrix
{
  static bool const value = false;
};

template <typename T, int N>
struct IsHomMatrix< T, N, typename std::enable_if< T::RowsAtCompileTime == N + 1 && T::ColsAtCompileTime == N + 1 >::type >
{
  static bool const value = true;
};

} // namespace TransformerInternal

template <typename ObjT, typename TransT, int N, typename ScalarT>
struct TransformerImpl<ObjT, TransT, N, ScalarT,
                       typename std::enable_if< IsNonReferencedPointN<ObjT, N>::value
                                             && !TransformerInternal::IsHomMatrix<TransT, N>::value >::type>
{
  typedef Vector<N, ScalarT> Result;
  static Result transform(ObjT const & obj, TransT const & tr) { return tr * PointTraitsN<ObjT, N, ScalarT>::getPosition(obj); }
};

// Homogeneous multiplication
template <typename ObjT, typename TransT, int N, typename ScalarT>
struct TransformerImpl<ObjT, TransT, N, ScalarT,
                       typename std::enable_if< IsNonReferencedPointN<ObjT, N>::value
                                             && TransformerInternal::IsHomMatrix<TransT, N>::value >::type>
{
  typedef Vector<N, ScalarT> Result;
  static Result transform(ObjT const & obj, TransT const & tr)
  { return Math::hmul(tr, PointTraitsN<ObjT, N, ScalarT>::getPosition(obj)); }
};

template <int N, typename ScalarT>
struct TransformerImpl< AxisAlignedBoxN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BoxN<N, ScalarT> Result;
  static Result transform(AxisAlignedBoxN<N, ScalarT> const & aabb, RigidTransformN<N, ScalarT> const & tr)
  { return Result(aabb, CoordinateFrameN<N, ScalarT>(tr)); }
};

template <int N, typename ScalarT>
struct TransformerImpl< BoxN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BoxN<N, ScalarT> Result;
  static Result transform(BoxN<N, ScalarT> const & box, RigidTransformN<N, ScalarT> const & tr)
  { return Result(box.getLocalAAB(), CoordinateFrameN<N, ScalarT>(tr) * box.getLocalFrame()); }
};

template <int N, typename ScalarT>
struct TransformerImpl< BallN<N, ScalarT>, RigidTransformN<N, ScalarT>, N, ScalarT >
{
  typedef BallN<N, ScalarT> Result;
  static Result transform(BallN<N, ScalarT> const & ball, RigidTransformN<N, ScalarT> const & tr)
  { return Result(tr * ball.getCenter(), ball.getRadius()); }
};

template <typename VertexTripleT, typename TransT, typename ScalarT>
struct THEA_API TransformerImpl< Triangle3<VertexTripleT>, TransT, 3, ScalarT >
{
  typedef LocalTriangle3 Result;
  static Result transform(Triangle3<VertexTripleT> const & tri, TransT const & tr)
  {
    return Result(Transformer::transform<3, ScalarT>(tri.getVertex(0), tr),
                          Transformer::transform<3, ScalarT>(tri.getVertex(1), tr),
                          Transformer::transform<3, ScalarT>(tri.getVertex(2), tr));
  }
};

template <int N, typename ScalarT, typename TransT>
struct THEA_API TransformerImpl< LineSegmentN<N, ScalarT>, TransT, N, ScalarT >
{
  typedef LineSegmentN<N, ScalarT> Result;
  static Result transform(LineSegmentN<N, ScalarT> const & seg, TransT const & tr)
  {
    return Result(Transformer::transform<N, ScalarT>(seg.getEndpoint(0), tr),
                  Transformer::transform<N, ScalarT>(seg.getEndpoint(1), tr));
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
