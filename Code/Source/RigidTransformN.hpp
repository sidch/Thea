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
// First version: 2011
//
//============================================================================

#ifndef __Thea_RigidTransformN_hpp__
#define __Thea_RigidTransformN_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"
#include <sstream>

namespace Thea {

// Forward declarations
template <int N, typename T> class RigidTransformN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for a rigid transformation in N-space, consisting of a rotation followed by a translation.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 *
 * @note While this is technically an affine transform, it restricts enough functionality to make a separate implementation
 * preferable. It can be trivially converted to an AffineTransformN using the implicit conversion operator or toAffine().
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ RigidTransformNBase
{
  public:
    typedef RigidTransformN<N, T>   RigidTransformT;   ///< N-dimensional rigid transform type.
    typedef AffineTransformN<N, T>  AffineTransformT;  ///< N-dimensional affine transform type.
    typedef Vector<N, T>            VectorT;           ///< N-dimensional vector.
    typedef Matrix<N, N, T>         MatrixT;           ///< NxN matrix.

    THEA_DECL_SMART_POINTERS(RigidTransformT)

  public:
    /** Default constructor. Constructs the identity transform. */
    RigidTransformNBase() : aff(AffineTransformT::identity()) {}

    /** Construct a translation. */
    static RigidTransformT translation(VectorT const & translation_)
    {
      return _fromAffine(AffineTransformT::translation(translation_));
    }

    /**
     * Construct from an affine transform, assuming it is rigid (<b>use with caution</b> since it can break the rigidity
     * guarantee).
     */
    static RigidTransformT _fromAffine(AffineTransformT const & aff_)
    {
      RigidTransformT rt;
      rt.aff = aff_;
      return rt;
    }

    /** Cast the transform to a different scalar type. */
    template <typename U> RigidTransformN<N, U> cast() const
    {
      return RigidTransformN<N, U>::_fromAffine(aff.template cast<U>());
    }

    /** Get rotation component. */
    MatrixT const & getRotation() const { return aff.getLinear(); }

    /** Set the rotation matrix (<b>use with caution</b> since it can break the rigidity guarantee). */
    void _setRotation(MatrixT const & rot) { aff.setLinear(rot); }

    /** Get translation component. */
    VectorT const & getTranslation() const { return aff.getTranslation(); }

    /** Set translation component. */
    void setTranslation(VectorT const & translation_) { aff.setTranslation(translation_); }

    /** Convert to an N x (N + 1) transformation matrix. */
    Matrix<N, N + 1, T> toMatrix() const { return aff.toMatrix(); }

    /** Convert to an (N + 1) x (N + 1) transformation matrix in homogeneous coordinates (last row is identity). */
    Matrix<N + 1, N + 1, T> homogeneous() const { return aff.homogeneous(); }

    /**
     * Convert to a general affine transform. This <i>may</i> be more efficient than the conversion operator since it explicitly
     * omits any possible copying, but the code doesn't look as clean.
     */
    AffineTransformT const & toAffine() const { return aff; }

    /**
     * Convert to a general affine transform. This <i>may</i> be less efficient than toAffine(), but the code looks cleaner
     * (e.g. in <code>a * b * c</code>, each of the variables can be either AffineTransformT or RigidTransformT).
     */
    operator AffineTransformT() const { return aff; }

    /** Get the inverse transform. */
    RigidTransformT inverse() const
    {
      MatrixT inv = aff.getLinear().transpose();  // transpose == inverse for orthogonal matrix
      return _fromAffine(AffineTransformT(inv, inv * (-aff.getTranslation())));
    }

    /** Get an element of the N x (N + 1) matrix representing this transform. */
    Real operator()(int i, int j) const { return aff(i, j); }

    /** Compose this transform with another. The other is applied first. */
    RigidTransformT operator*(RigidTransformT const & rhs) const { return _fromAffine(aff * rhs.aff); }

    /** Apply this transform to a N-vector. */
    VectorT operator*(VectorT const & v) const { return aff * v; }

    /** Get a string representing the transform. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[R: " << Thea::toString(getRotation()) << ", T: " << Thea::toString(getTranslation()) << ']';
      return oss.str();
    }

    /** Get the identity transform. */
    static RigidTransformT const & identity()
    {
      static RigidTransformT const idty = _fromAffine(AffineTransformT::identity());
      return idty;
    }

  private:
    AffineTransformT aff;  ///< The wrapped affine transformation.

}; // class RigidTransformNBase

} // namespace Internal

/**
 * A rigid transformation in N-space, consisting of a rotation followed by a translation.
 *
 * @note While this is technically an affine transform, it restricts enough functionality to make a separate implementation
 * preferable. It can be trivially converted to an AffineTransformN using the implicit conversion operator or toAffine().
 */
template <int N, typename T = Real>
class /* THEA_API */ RigidTransformN : public Internal::RigidTransformNBase<N, T>
{
  private:
    typedef Internal::RigidTransformNBase<N, T> BaseT;

  public:
    typedef typename BaseT::AffineTransformT  AffineTransformT;
    typedef typename BaseT::VectorT           VectorT;
    typedef typename BaseT::MatrixT           MatrixT;

    /** Default constructor. Constructs the identity transform. */
    RigidTransformN() {}

}; // class RigidTransformN

/** Pipe a textual representation of a rigid transform to a <code>std::ostream</code>. */
template <int N, typename T>
std::ostream &
operator<<(std::ostream & os, RigidTransformN<N, T> const & tr)
{
  return os << tr.toString();
}

} // namespace Thea

#include "RigidTransform2.hpp"
#include "RigidTransform3.hpp"

#endif
