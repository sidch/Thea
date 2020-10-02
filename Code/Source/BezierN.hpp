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
// First version: 2018
//
//============================================================================
//
// Portions of this source file were derived from the Bezier curve
// implementation of:
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
//
//============================================================================

#ifndef __Thea_BezierN_hpp__
#define __Thea_BezierN_hpp__

#include "Common.hpp"
#include "MatVec.hpp"
#include "SplineN.hpp"

namespace Thea {

/**
 * A Bezier curve segment in N-dimensional space, parametrized by a scalar in the range [0, 1].
 *
 * This implementation uses code from Dave Eberly's Geometric Tools library.
 */
template <int N, typename T = Real>
class /* THEA_API */ BezierN : public SplineN<N, T>
{
  private:
    typedef SplineN<N, T> BaseT;  ///< Base curve type.

  public:
    THEA_DECL_SMART_POINTERS(BezierN)

    typedef typename BaseT::VectorT VectorT;  ///< N-dimensional vector type.

    /**
     * Construct an (initially zero length) Bezier curve of a given order (2 for quadratic Bezier, 3 for cubic Bezier, etc). The
     * segment will be initialized with \a order + 1 control vectors, all initially zero.
     */
    BezierN(intx order_ = 3) : BaseT(0, 1)
    {
      alwaysAssertM(order_ >= 1, "BezierN: Order must be non-negative");

      ctrl[0].resize((size_t)order_ + 1, VectorT::Zero());
      this->setChanged(true);
    }

    intx getOrder() const { return (intx)ctrl[0].size() - 1; }

    intx numControls() const { return (intx)ctrl[0].size(); }

    VectorT const & getControl(intx index) const
    {
      alwaysAssertM(index >= 0 && index < (intx)ctrl[0].size(), "BezierN: Control point index out of range");
      return ctrl[0][(size_t)index];
    }

    void setControl(intx index, VectorT const & pos)
    {
      alwaysAssertM(index >= 0 && index < (intx)ctrl[0].size(), "BezierN: Control point index out of range");

      size_t i = (size_t)index;
      ctrl[0][i] = pos;

      this->setChanged(true);
    }

    bool hasDeriv(intx deriv_order) const { return (deriv_order >= 0); }

  private:
    mutable Array<VectorT>   ctrl[4];   ///< Arrays of curve control vectors and first, second and third-order differences.
    mutable MatrixX<float64>  binom;    ///< Cached binomial coefficients.

    /** Cache binomial coefficients for computing Bernstein polynomials. */
    void cacheBinom() const
    {
      if (binom.rows() > 0)
        return;

      intx n = (intx)ctrl[0].size();
      binom.resize(n, n);

      // From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteBezierCurve.h
      //
      // Compute combinatorial values binom(n, k). The values binom(r, c) are invalid for r < c; that is, we use only the
      // entries for r >= c.
      binom(0, 0) = 1;
      binom(1, 0) = 1;
      binom(1, 1) = 1;
      for (intx i = 2; i < n; ++i)
      {
        binom(i, 0) = 1;
        binom(i, i) = 1;
        for (intx j = 1; j < i; ++j)
          binom(i, j) = binom(i - 1, j - 1) + binom(i - 1, j);
      }
    }

    void update() const
    {
      if (!this->isChanged()) return;

      size_t n = ctrl[0].size();

      // Compute first-order differences
      ctrl[1].resize(n - 1);
      for (size_t i = 0; i < n - 1; ++i)
        ctrl[1][i] = ctrl[0][i + 1] - ctrl[0][i];

      // Compute second-order differences
      if (n > 2)
      {
        ctrl[2].resize(n - 2);
        for (size_t i = 0; i < n - 2; ++i)
          ctrl[2][i] = ctrl[1][i + 1] - ctrl[1][i];
      }

      // Compute third-order differences
      if (n > 3)
      {
        ctrl[3].resize(n - 3);
        for (size_t i = 0; i < n - 3; ++i)
          ctrl[3][i] = ctrl[2][i + 1] - ctrl[2][i];
      }

      cacheBinom();
      this->setChanged(false);
    }

    VectorT eval(T const & t, intx deriv_order) const
    {
      alwaysAssertM(t >= -0.00001 && t <= 1.00001, format("BezierN: Curve parameter %lf out of range", static_cast<double>(t)));
      alwaysAssertM(deriv_order >= 0, format("BezierN: Invalid derivative order %ld", deriv_order));

      intx order = getOrder();
      if (deriv_order > order) return VectorT::Zero();

      update();

      T omt = 1 - t;
      VectorT result = omt * ctrl[deriv_order][0];

      T tpow = t;
      intx isup = order - deriv_order;
      for (intx i = 1; i < isup; ++i)
      {
        T c = static_cast<T>(binom(isup, i) * tpow);
        result = (result + c * ctrl[deriv_order][i]) * omt;
        tpow *= t;
      }
      result = (result + tpow * ctrl[deriv_order][isup]);

      intx multiplier = 1;
      for (intx i = 0; i < deriv_order; ++i)
        multiplier *= (order - i);

      result *= multiplier;

      return result;
    }

    void getBasisFunctions(float64 t, VectorX<float64> & b) const
    {
      // Bernstein basis

      cacheBinom();

      intx n = getOrder();
      b.resize(n + 1);

      b[0] = 1.0;
      float64 tpow = t;
      for (intx i = 1; i <= n; ++i, tpow *= t)
        b[i] = binom(n, i) * tpow;

      float64 omt = 1 - t;
      float64 omt_pow = omt;
      for (intx i = n - 1; i >= 0; --i, omt_pow *= omt)
        b[i] *= omt_pow;
    }

    bool firstAndLastControlsArePositions() const { return true; }

}; // class BezierN

} // namespace Thea

#endif
