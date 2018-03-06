//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2018, Siddhartha Chaudhuri
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
#include "Matrix.hpp"
#include "SplineN.hpp"

namespace Thea {

/**
 * A Bezier curve segment in N-dimensional space, parametrized by a scalar in the range [0, 1].
 *
 * This implementation uses code from Dave Eberly's Geometric Tools library.
 */
template <long N, typename T>
class /* THEA_API */ BezierN : public SplineN<N, T>
{
  private:
    typedef SplineN<N, T> BaseT;  ///< Base curve type.

  public:
    typedef typename BaseT::VectorT VectorT;  ///< N-dimensional vector type.

    /**
     * Construct an (initially zero length) Bezier curve of a given order (2 for quadratic Bezier, 3 for cubic Bezier, etc). The
     * segment will be initialized with \a order + 1 control vectors, all initially zero.
     */
    BezierN(long order_ = 3) : BaseT(0, 1)
    {
      alwaysAssertM(order_ >= 1, "BezierN: Order must be non-negative");

      ctrl[0].resize((array_size_t)order_ + 1, VectorT::zero());
      this->setChanged(true);
    }

    long getOrder() const { return (long)ctrl[0].size() - 1; }

    long numControls() const { return (long)ctrl[0].size(); }

    VectorT const & getControl(long index) const
    {
      alwaysAssertM(index >= 0 && index < (long)ctrl[0].size(), "BezierN: Control point index out of range");
      return ctrl[0][(array_size_t)index];
    }

    void setControl(long index, VectorT const & pos)
    {
      alwaysAssertM(index >= 0 && index < (long)ctrl[0].size(), "BezierN: Control point index out of range");

      array_size_t i = (array_size_t)index;
      ctrl[0][i] = pos;

      this->setChanged(true);
    }

    bool hasDeriv(long deriv_order) const { return (deriv_order >= 0 && deriv_order <= 3); }

  private:
    mutable TheaArray<VectorT>  ctrl[4];  ///< Arrays of curve control vectors and first, second and third-order differences.
    mutable Matrix<double>      binom;    ///< Cached binomial coefficients.

    /** Cache binomial coefficients for computing Bernstein polynomials. */
    void cacheBinom() const
    {
      if (binom.numRows() > 0)
        return;

      long n = (long)ctrl[0].size();
      binom.resize(n, n);

      // From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteBezierCurve.h
      //
      // Compute combinatorial values binom(n, k). The values binom(r, c) are invalid for r < c; that is, we use only the
      // entries for r >= c.
      binom(0, 0) = 1;
      binom(1, 0) = 1;
      binom(1, 1) = 1;
      for (long i = 2; i < n; ++i)
      {
        binom(i, 0) = 1;
        binom(i, i) = 1;
        for (long j = 1; j < i; ++j)
          binom(i, j) = binom(i - 1, j - 1) + binom(i - 1, j);
      }
    }

    void update() const
    {
      if (!this->isChanged()) return;

      array_size_t n = ctrl[0].size();

      // Compute first-order differences
      ctrl[1].resize(n - 1);
      for (array_size_t i = 0; i < n - 1; ++i)
        ctrl[1][i] = ctrl[0][i + 1] - ctrl[0][i];

      // Compute second-order differences
      if (n > 2)
      {
        ctrl[2].resize(n - 2);
        for (array_size_t i = 0; i < n - 2; ++i)
          ctrl[2][i] = ctrl[1][i + 1] - ctrl[1][i];
      }

      // Compute third-order differences
      if (n > 3)
      {
        ctrl[3].resize(n - 3);
        for (array_size_t i = 0; i < n - 3; ++i)
          ctrl[3][i] = ctrl[2][i + 1] - ctrl[2][i];
      }

      cacheBinom();
      this->setChanged(false);
    }

    VectorT eval(T const & t, long deriv_order) const
    {
      alwaysAssertM(t >= -0.00001 && t <= 1.00001, format("BezierN: Curve parameter %lf out of range", static_cast<double>(t)));
      alwaysAssertM(deriv_order >= 0 && deriv_order <= 3, format("BezierN: Invalid derivative order %ld", deriv_order));

      update();

      T omt = 1 - t;
      VectorT result = omt * ctrl[deriv_order][0];

      long order = getOrder();
      T tpow = t;
      long isup = order - deriv_order;
      for (long i = 1; i < isup; ++i)
      {
        T c = static_cast<T>(binom(isup, i) * tpow);
        result = (result + c * ctrl[deriv_order][i]) * omt;
        tpow *= t;
      }
      result = (result + tpow * ctrl[deriv_order][isup]);

      long multiplier = 1;
      for (long i = 0; i < deriv_order; ++i)
        multiplier *= (order - i);

      result *= multiplier;

      return result;
    }

    void getBasisFunctions(double t, TheaArray<double> & b) const
    {
      // Bernstein basis

      cacheBinom();

      long n = getOrder();
      b.resize((array_size_t)n + 1);

      b[0] = 1.0;
      double tpow = t;
      for (long i = 1; i <= n; ++i, tpow *= t)
        b[(array_size_t)i] = binom(n, i) * tpow;

      double omt = 1 - t;
      double omt_pow = omt;
      for (long i = n - 1; i >= 0; --i, omt_pow *= omt)
        b[(array_size_t)i] *= omt_pow;
    }

    bool firstAndLastControlsArePositions() const { return true; }

}; // class BezierN

} // namespace Thea

#endif
