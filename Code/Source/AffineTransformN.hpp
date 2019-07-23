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

#ifndef __Thea_AffineTransformN_hpp__
#define __Thea_AffineTransformN_hpp__

#include "Common.hpp"
#include "MatVec.hpp"
#include <sstream>

namespace Thea {

// Forward declarations
template <int N, typename T> class AffineTransformN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for affine transform in N-dimensional space, where N is any <b>positive</b> (non-zero) integer
 * and T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ AffineTransformNBase
{
  public:
    typedef AffineTransformN<N, T>  AffineTransformT;  ///< N-dimensional affine transform type.
    typedef Vector<N, T>            VectorT;           ///< N-dimensional vector.
    typedef Matrix<N, N, T>         MatrixT;           ///< NxN matrix.

    THEA_DECL_SMART_POINTERS(AffineTransformT)

    /** Default constructor. Does not initialize anything. */
    AffineTransformNBase() {}

    /** Construct from a linear transform, followed by a translation. */
    AffineTransformNBase(MatrixT const & linear_, VectorT const & translation_) : linear(linear_), trans(translation_) {}

    /** Construct a scaling transform. */
    static AffineTransformT scaling(VectorT const & s)
    {
      return AffineTransformT(MatrixT(s.asDiagonal()), VectorT::Zero());
    }

    /** Construct a uniform scaling transform. */
    static AffineTransformT scaling(T const & s)
    {
      return scaling(VectorT::Constant(s));
    }

    /** Construct a translation. */
    static AffineTransformT translation(VectorT const & v)
    {
      return AffineTransformT(MatrixT::Identity(), v);
    }

    /** Get linear transform component. */
    MatrixT const & getLinear() const { return linear; }

    /** Get linear transform component. */
    MatrixT & getLinear() { return linear; }

    /** Set linear transform component. */
    void setLinear(MatrixT const & linear_) { linear = linear_; }

    /** Get translation component. */
    VectorT const & getTranslation() const { return trans; }

    /** Get translation component. */
    VectorT & getTranslation() { return trans; }

    /** Set translation component. */
    void setTranslation(VectorT const & translation_) { trans = translation_; }

    /** Convert to an (N + 1) x (N + 1) transformation matrix in homogeneous coordinates (last row is identity). */
    Matrix<N + 1, N + 1, T> homogeneous() const
    {
      Matrix<N + 1, N + 1, T> m = Matrix<N + 1, N + 1, T>::Identity();
      for (intx i = 0; i < N; ++i)
      {
        for (intx j = 0; j < N; ++j)
          m(i, j) = linear(i, j);

        m(i, N) = trans[i];
      }

      return m;
    }

    /** Convert to an N x (N + 1) transformation matrix. */
    Matrix<N, N + 1, T> toMatrix() const
    {
      Matrix<N, N + 1, T> m;
      for (intx i = 0; i < N; ++i)
      {
        for (intx j = 0; j < N; ++j)
          m(i, j) = linear(i, j);

        m(i, N) = trans[i];
      }

      return m;
    }

    /** Get the inverse transform. */
    AffineTransformT inverse() const
    {
      MatrixT inv = linear.inverse();
      return AffineTransformT(inv, inv * (-trans));
    }

    /** Get an element of the N x (N + 1) matrix representing this transform. */
    T operator()(intx i, intx j) const
    {
      debugAssertM(i >= 0 && i < N && j >= 0 && j <= N, "AffineTransformT: Index out of bounds");
      return j == N ? trans[i] : linear(i, j);
    }

    /** Get an element of the N x (N + 1) matrix representing this transform. */
    T & operator()(intx i, intx j)
    {
      debugAssertM(i >= 0 && i < N && j >= 0 && j <= N, "AffineTransformT: Index out of bounds");
      return j == N ? trans[i] : linear(i, j);
    }

    /** Compose this transform with another. The other is applied first. */
    AffineTransformT operator*(AffineTransformT const & rhs) const
    {
      return AffineTransformT(linear * rhs.linear, linear * rhs.trans + trans);
    }

    /** Apply this transform to a vector. */
    VectorT operator*(VectorT const & v) const { return linear * v + trans; }

    /** Get a string representing the transform. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[L: " << Thea::toString(linear) << ", T: " << Thea::toString(trans) << ']';
      return oss.str();
    }

    /** Get the identity transform. */
    static AffineTransformT const & identity()
    {
      static AffineTransformT const idty(MatrixT::Identity(), VectorT::Zero());
      return idty;
    }

  private:
    MatrixT linear;  ///< Linear component.
    VectorT trans;   ///< Translation component.

}; // class AffineTransformNBase

} // namespace Internal

/** An affine transform in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <int N, typename T = Real>
class /* THEA_API */ AffineTransformN : public Internal::AffineTransformNBase<N, T>
{
  private:
    typedef Internal::AffineTransformNBase<N, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;
    typedef typename BaseT::MatrixT MatrixT;

    /** Default constructor. Does not initialize anything. */
    AffineTransformN() {}

    /** Construct from a linear transform, followed by a translation. */
    AffineTransformN(MatrixT const & linear_, VectorT const & translation_ = VectorT::Zero()) : BaseT(linear_, translation_) {}

    /** Copy constructor. */
    AffineTransformN(AffineTransformN const & src) : BaseT(src) {}

}; // class AffineTransformN

/** Pipe a textual representation of an affine transform to a <code>std::ostream</code>. */
template <int N, typename T>
std::ostream &
operator<<(std::ostream & os, AffineTransformN<N, T> const & tr)
{
  return os << tr.toString();
}

} // namespace Thea

#include "AffineTransform2.hpp"
#include "AffineTransform3.hpp"

#endif
