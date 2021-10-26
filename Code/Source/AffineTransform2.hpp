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

#ifndef __Thea_AffineTransform2_hpp__
#define __Thea_AffineTransform2_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"

namespace Thea {

/** An affine transformation in 2-space, consisting of a linear transformation (2x2 matrix) plus a translation. */
template <typename T>
class /* THEA_API */ AffineTransformN<2, T> : public Internal::AffineTransformNBase<2, T>
{
  private:
    typedef Internal::AffineTransformNBase<2, T> BaseT;

  public:
    typedef typename BaseT::VectorT  VectorT;
    typedef typename BaseT::MatrixT  MatrixT;

    /** Default constructor. Does not initialize anything. */
    AffineTransformN() {}

    /** Construct from a linear transform, followed by a translation. */
    AffineTransformN(MatrixT const & linear_, VectorT const & translation_ = VectorT::Zero()) : BaseT(linear_, translation_) {}

    /** Construct from a 2 x 3 matrix. */
    template <typename AffineMatrixT> AffineTransformN(Eigen::DenseBase<AffineMatrixT> const & m) : BaseT(m) {}

    /** Copy constructor. */
    AffineTransformN(AffineTransformN const & src) : BaseT(src) {}

    /**
     * Construct from 2 basis vectors, specifying the linear transform, and a translation. The three arguments form the columns
     * of the 2x3 matrix specifying the transform.
     */
    AffineTransformN(VectorT const & x, VectorT const & y, VectorT const & translation_)
    : BaseT((MatrixT() << x, y).finished(), translation_) {}

    /** Construct from a 2x3 array. */
    AffineTransformN(T const &  m00, T const &  m01, T const &  m02,
                     T const &  m10, T const &  m11, T const &  m12)
    : BaseT((MatrixT() << m00, m01, m10, m11).finished(), VectorT(m02, m12)) {}

    using BaseT::scaling;

    /** Construct a scaling transform. */
    static AffineTransformN scaling(T const & sx, T const & sy)
    {
      return BaseT::scaling(VectorT(sx, sy));
    }

    using BaseT::translation;

    /** Construct a translation. */
    static AffineTransformN translation(T const & tx, T const & ty)
    {
      return BaseT::translation(VectorT(tx, ty));
    }

    /** Construct a rotation specified by an angle (in radians) around the origin. */
    static AffineTransformN rotation(Real radians)
    {
      return AffineTransformN(Math::rotation(radians), VectorT::Zero());
    }

}; // class AffineTransformN<2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AffineTransformN<2, Real>;
#endif

/** The default affine transform class in real 2-space. */
typedef AffineTransformN<2, Real> AffineTransform2;

} // namespace Thea

#endif
