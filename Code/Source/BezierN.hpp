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
// implementations of:
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
//
// and
//
// An Algorithm for Automatically Fitting Digitized Curves
// by Philip J. Schneider
// from "Graphics Gems", Academic Press, 1990
//
//============================================================================

#ifndef __Thea_BezierN_hpp__
#define __Thea_BezierN_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "Matrix.hpp"
#include "VectorN.hpp"
#include "Algorithms/PointTraitsN.hpp"
#include <iterator>

namespace Thea {

/**
 * A Bezier curve segment in N-dimensional space.
 *
 * This implementation uses code from Dave Eberly's Geometric Tools library, and Phil Schneider's Graphics Gems Bezier-fitting
 * module.
 */
template <long N, typename T>
class /* THEA_API */ BezierN
{
  public:
    typedef VectorN<N, T> VectorT;  ///< N-dimensional vector type.

    /**
     * Construct an (initially zero length) Bezier curve of a given order (2 for quadratic Bezier, 3 for cubic Bezier, etc). The
     * segment will be initialized with \a order + 1 control points, all initially zero.
     */
    BezierN(int order_ = 3) : changed(true)
    {
      alwaysAssertM(order_ >= 1, "BezierN: Order must be non-negative");

      ctrl[0].resize((array_size_t)order_ + 1, VectorT::zero());
    }

    /** Get the order of the curve. */
    int getOrder() const { return (int)ctrl[0].size() - 1; }

    /**
     * Get a control point of the curve.
     *
     * @param index The index of the control point, in the range 0 to getOrder() (inclusive).
     */
    VectorT const & getControlPoint(int index) const
    {
      debugAssertM(index >= 0 && index < (int)ctrl[0].size(), "BezierN: Control point index out of range");
      return ctrl[0][(array_size_t)index];
    }

    /**
     * Set a control point of the curve.
     *
     * @param index The index of the control point to be set, in the range 0 to getOrder() (inclusive).
     * @param pos The new position of the control point.
     */
    void setControlPoint(int index, VectorT const & pos)
    {
      debugAssertM(index >= 0 && index < (int)ctrl[0].size(), "BezierN: Control point index out of range");

      array_size_t i = (array_size_t)index;
      ctrl[0][i] = pos;
      changed = true;
    }

    /** Get the point on the curve with parameter value \a t, in the range [0, 1]. */
    VectorT getPoint(T const & t) const
    {
      return eval(t, 0);
    }

    /**
     * Get the first, second, or third derivative (specified by \a deriv_order = 1, 2 or 3) of the curve at parameter value
     * \a t, in the range [0, 1].
     */
    VectorT getDeriv(T const & t, int deriv_order = 1) const
    {
      return eval(t, deriv_order);
    }

    /**
     * Fit the Bezier curve segment to a sequence of points [begin, end), where the curve begins at the first point and
     * ends at the last one. Currently only quadratic and cubic Beziers are supported.
     *
     * @return The non-negative fitting error on success, a negative value on error.
     */
    template <typename InputIterator>
    T fitToPoints(InputIterator begin, InputIterator end,
                  typename boost::enable_if< Algorithms::IsPointN<typename std::iterator_traits<InputIterator>::value_type, N>
                                           >::type * dummy = NULL)
    {
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      switch (getOrder())
      {
        case 2  : return fitQuadratic(begin, end);
        case 3  : return fitCubic(begin, end);

        default :
        {
          THEA_ERROR << "BezierN: Currently only Bezier curves of order 2 or 3 can be fitted to points";
          return -1;
        }
      }
    }

  private:
    mutable TheaArray<VectorT>  ctrl[4];  ///< Arrays of curve control points and first, second and third-order differences.
    mutable Matrix<T>           binom;    ///< Cached binomial coefficients.
    mutable bool                changed;  ///< Was the curve changed?

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

    /** Updated cached data like the finite differences of control points and binomial coefficients. */
    void update() const
    {
      if (!changed) return;

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
      changed = false;
    }

    /** Evaluate the curve, or one of its derivatives, at a given parameter value. */
    VectorT eval(T const & t, int deriv_order) const
    {
      alwaysAssertM(t >= -0.00001 && t <= 1.00001, "BezierN: Curve parameter out of range");
      alwaysAssertM(deriv_order >= 0 && deriv_order <= 3, "BezierN: Invalid derivative order");

      update();

      T omt = 1 - t;
      VectorT result = omt * ctrl[deriv_order][0];

      int order = getOrder();
      T tpow = t;
      int isup = order - deriv_order;
      for (int i = 1; i < isup; ++i)
      {
        T c = binom(isup, i) * tpow;
        result = (result + c * ctrl[deriv_order][i]) * omt;
        tpow *= t;
      }
      result = (result + tpow * ctrl[deriv_order][isup]);

      int multiplier = 1;
      for (int i = 0; i < deriv_order; ++i)
        multiplier *= (order - i);

      result *= multiplier;

      return result;
    }

    template <typename InputIterator> T fitQuadratic(InputIterator begin, InputIterator end)
    {
      return -1;
    }

    template <typename InputIterator> T fitCubic(InputIterator begin, InputIterator end)
    {
      return -1;
    }

}; // class BezierN

} // namespace Thea

#endif
