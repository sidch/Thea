//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_RigidTransformN_hpp__
#define __Thea_RigidTransformN_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"
#include <sstream>

namespace Thea {

// Forward declarations
template <long N, typename T> class RigidTransformN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for a rigid transformation in N-space, consisting of a rotation followed by a translation.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 *
 * @note While this is technically an affine transform, it restricts enough functionality to make a separate implementation
 * preferable. It can be trivially converted to an AffineTransformN using the implicit conversion operator or toAffine().
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ RigidTransformNBase
{
  public:
    typedef RigidTransformN<N, T>   RigidTransformT;  ///< N-dimensional rigid transform type.
    typedef AffineTransformN<N, T>  AffineTransformT;  ///< N-dimensional affine transform type.
    typedef VectorN<N, T>           VectorT;  ///< N-dimensional vector.
    typedef MatrixMN<N, N, T>       MatrixT;  ///< NxN matrix.

    THEA_DEF_POINTER_TYPES(RigidTransformT, shared_ptr, weak_ptr)

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

    /** Get rotation component. */
    MatrixT const & getRotation() const { return aff.getLinear(); }

    /** Set the rotation matrix (<b>use with caution</b> since it can break the rigidity guarantee). */
    void _setRotation(MatrixT const & rot) { aff.setLinear(rot); }

    /** Get translation component. */
    VectorT const & getTranslation() const { return aff.getTranslation(); }

    /** Set translation component. */
    void setTranslation(VectorT const & translation_) { aff.setTranslation(translation_); }

    /** Convert to an N x (N + 1) transformation matrix. */
    MatrixMN<N, N + 1, T> toMatrix() const { return aff.toMatrix(); }

    /** Convert to an (N + 1) x (N + 1) transformation matrix in homogenous coordinates (last row is identity). */
    MatrixMN<N + 1, N + 1, T> toHomMatrix() const { return aff.toHomMatrix(); }

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
      oss << "[R: " << getRotation() << ", T: " << getTranslation() << ']';
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
template <long N, typename T = Real>
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
template <long N, typename T>
std::ostream &
operator<<(std::ostream & os, RigidTransformN<N, T> const & tr)
{
  return os << tr.toString();
}

} // namespace Thea

#include "RigidTransform2.hpp"
#include "RigidTransform3.hpp"

#endif
