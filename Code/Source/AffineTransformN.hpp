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

    /** Cast the transform to a different scalar type. */
    template <typename U> AffineTransformN<N, U> cast() const
    {
      return AffineTransformN<N, U>(linear.template cast<U>(), trans.template cast<U>());
    }

    /** Construct a scaling transform. */
    static AffineTransformT scaling(VectorT const & s)
    {
      return AffineTransformT(Math::scaling(s), VectorT::Zero());
    }

    /** Construct a uniform scaling transform. */
    static AffineTransformT scaling(T const & s)
    {
      return AffineTransformT(Math::scaling<N>(s), VectorT::Zero());;
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
