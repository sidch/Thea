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
#include "Algorithms/FastCopy.hpp"
#include "Algorithms/LinearLeastSquares.hpp"
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
    BezierN(long order_ = 3) : changed(true)
    {
      alwaysAssertM(order_ >= 1, "BezierN: Order must be non-negative");

      ctrl[0].resize((array_size_t)order_ + 1, VectorT::zero());
    }

    /** Get the order of the curve. */
    long getOrder() const { return (long)ctrl[0].size() - 1; }

    /** Get the number of control points of the curve. Necessarily equal to getOrder() + 1. */
    long numControlPoints() const { return (long)ctrl[0].size(); }

    /**
     * Get a control point of the curve.
     *
     * @param index The index of the control point, in the range 0 to getOrder() (inclusive).
     */
    VectorT const & getControlPoint(long index) const
    {
      alwaysAssertM(index >= 0 && index < (long)ctrl[0].size(), "BezierN: Control point index out of range");
      return ctrl[0][(array_size_t)index];
    }

    /**
     * Set a control point of the curve.
     *
     * @param index The index of the control point to be set, in the range 0 to getOrder() (inclusive).
     * @param pos The new position of the control point.
     */
    void setControlPoint(long index, VectorT const & pos)
    {
      alwaysAssertM(index >= 0 && index < (long)ctrl[0].size(), "BezierN: Control point index out of range");

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
    VectorT getDeriv(T const & t, long deriv_order = 1) const
    {
      return eval(t, deriv_order);
    }

    /**
     * Fit the Bezier curve segment to a sequence of points [begin, end), where the curve begins at the first point and
     * ends at the last one.
     *
     * @param begin First point in the sequence to be fitted. This will be the position of the first control point of the curve.
     * @param end One past the last point in the sequence to be fitted. The last point will be the position of the last control
     *   point of the curve.
     * @param initial_params The curve parameters of the points, if known in advance. The first and last entries will be
     *   ignored, since these are always assumed to be 0 and 1 respectively.
     * @param max_reparam_iters If greater than zero, the parameter values of the points will be re-estimated (at most) this
     *  many times, guided by initial values (\a initial_params) if any.
     * @param num_reparam_steps_per_iter The number of successive Newton-Raphson steps in each iteration of reparametrization.
     * @param final_params If non-null, used to return the final parameter values of the point sequence. Must be pre-allocated
     *  to have at least as many entries as the number of points. The first and last values will always be 0 and 1,
     *  respectively, if the function succeeds.
     *
     * @return The non-negative squared fitting error on success, or a negative value on error.
     */
    template <typename InputIterator>
    double fitToPoints(InputIterator begin, InputIterator end,
                       T const * initial_params = NULL,
                       T * final_params = NULL,
                       long max_reparam_iters = 0,
                       long num_reparam_steps_per_iter = 1,  // conservative choice, following Graphics Gems
                       typename boost::enable_if<
                                  Algorithms::IsPointN<typename std::iterator_traits<InputIterator>::value_type, N>
                                                >::type * dummy = NULL)
    {
      // Initialize parameter values for the points
      TheaArray<double> u;
      if (initial_params)
      {
        T const * up = initial_params;
        for (InputIterator pi = begin; pi != end; ++pi, ++up)
          u.push_back(static_cast<double>(*up));
      }
      else
        chordLengthParametrize(begin, end, u);

      if (u.size() < ctrl[0].size())
      {
        THEA_ERROR << "BezierN: Cannot fit to fewer data points than control points";
        return -1;
      }

      // Fit the curve iteratively, alternating between a linear least squares fit with known parameter values, and
      // re-estimation of parameters via Newton-Raphson
      double sqerr = -1;
      while (true)
      {
        double e = llsqFit(begin, end, u);
        if (e < 0)  // could not fit, revert to last solution
          break;

        if (sqerr >= 0 && e > sqerr)  // error increased, stop iterating
          break;

        sqerr = e;

        // A bit wasteful to save it every iteration instead of outside the loop, but do it in case llsqFit fails and we have to
        // revert to the last solution
        if (final_params)
          Algorithms::fastCopy(u.begin(), u.end(), final_params);

        if (--max_reparam_iters < 0)
          break;

        refineParameters(begin, end, u, num_reparam_steps_per_iter);
      }

      return sqerr;
    }

    /** Get a textual representation of the curve. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[order = " << getOrder() << ", ctrl = [" << stringJoin(ctrl[0], ", ") << "]]";
      return oss.str();
    }

  private:
    mutable TheaArray<VectorT>  ctrl[4];  ///< Arrays of curve control points and first, second and third-order differences.
    mutable Matrix<double>      binom;    ///< Cached binomial coefficients.
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
    VectorT eval(T const & t, long deriv_order) const
    {
      alwaysAssertM(t >= -0.00001 && t <= 1.00001, "BezierN: Curve parameter out of range");
      alwaysAssertM(deriv_order >= 0 && deriv_order <= 3, "BezierN: Invalid derivative order");

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

    /** Estimate curve parameters for a sequence of points, by accumulating pairwise segment lengths along the sequence. */
    template <typename InputIterator>
    static void chordLengthParametrize(InputIterator begin, InputIterator end, TheaArray<double> & u)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      u.clear();

      if (begin == end)
        return;

      u.push_back(0.0);

      VectorT last = PointTraitsN<PointT, N, T>::getPosition(*begin);
      InputIterator pi = begin;
      double sum_u = 0.0;
      for (++pi; pi != end; ++pi)
      {
        VectorT curr = PointTraitsN<PointT, N, T>::getPosition(*pi);
        sum_u += static_cast<double>((curr - last).fastLength());
        u.push_back(sum_u);
        last = curr;
      }

      double umax = u[u.size() - 1];
      if (umax <= 0)
        return;

      for (array_size_t i = 0; i < u.size(); ++i)
        u[i] /= umax;
    }

    /** Get the Bernstein basis coefficients for a given curve parameter \a t. */
    void getBernsteinBasis(double t, TheaArray<double> & b) const
    {
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

    /**
     * Fit the curve to a sequence of points using linear least-squares. Control point positions are estimated by minimizing the
     * sum of squared errors between the points and their corresponding curve points with known parameters \a u.
     *
     * @return The sum of squared fitting errors, or a negative value on error.
     */
    template <typename InputIterator> double llsqFit(InputIterator begin, InputIterator end, TheaArray<double> const & u)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      array_size_t num_ctrls = ctrl[0].size();
      array_size_t num_unknown_ctrls = num_ctrls - 2;  // first and last points are known
      array_size_t num_unknowns = (array_size_t)(N * num_unknown_ctrls);

      LinearLeastSquares llsq((long)num_unknowns);
      TheaArray<double> basis;
      TheaArray<double> coeffs(num_unknowns, 0.0);

      // Cache the first and last points of the sequence
      InputIterator first = begin;
      InputIterator last = begin;
      for (InputIterator pi = begin; pi != end; ++pi)
        last = pi;

      VectorT first_pos  =  PointTraitsN<PointT, N, T>::getPosition(*first);
      VectorT last_pos   =  PointTraitsN<PointT, N, T>::getPosition(*last);

      // Try to make each point in the sequence match the point on the curve with the corresponding parameters
      array_size_t i = 0;
      for (InputIterator pi = begin; pi != end; ++pi, ++i)
      {
        getBernsteinBasis(u[i], basis);
        VectorT d = PointTraitsN<PointT, N, T>::getPosition(*pi)
                  - basis[0] * first_pos  // first control point of curve is fixed at starting point of sequence
                  - basis[basis.size() - 1] * last_pos;  // last control point is also fixed to last point of sequence

        // One scalar objective for each dimension
        for (long j = 0; j < N; ++j)
        {
          array_size_t offset = (array_size_t)j * num_unknown_ctrls;
          fastCopy(&basis[1], &basis[basis.size() - 1], coeffs.begin() + offset);

          llsq.addObjective(&coeffs[0], d[j]);

          // Reset the range
          fill(coeffs.begin() + offset, coeffs.begin() + offset + num_unknown_ctrls, 0.0);
        }
      }

      if (!llsq.solve(LinearLeastSquares::Constraint::UNCONSTRAINED))
      {
        THEA_ERROR << "BezierN: Could not solve linear least-squares curve fitting problem";
        return -1;
      }

      TheaArray<double> const & sol = llsq.getSolution();
      alwaysAssertM(sol.size() == num_unknowns,
                    "BezierN: Linear least-squares solution length doesn't match number of unknowns");

      ctrl[0][0] = first_pos;
      for (long i = 0; i < N; ++i)
      {
        array_size_t offset = (array_size_t)(i * num_unknown_ctrls);

        for (array_size_t j = 1; j <= num_unknown_ctrls; ++j)
          ctrl[0][j][i] = static_cast<T>(sol[offset + (j - 1)]);
      }
      ctrl[0][num_ctrls - 1] = last_pos;

      changed = true;
      return llsq.squaredError();
    }

    /**
     * Optimize each parameter value to bring the curve closer to the corresponding target point, using Newton-Raphson steps.
     *
     * See "An Algorithm for Automatically Fitting Digitized Curves", Philip J. Schneider, "Graphics Gems", 1990.
     */
    template <typename InputIterator>
    bool refineParameters(InputIterator begin, InputIterator end, TheaArray<double> & u, long num_newton_iters)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      array_size_t i = 0;
      for (InputIterator pi = begin; pi != end; ++pi, ++i)
      {
        // Target point
        VectorT p = PointTraitsN<PointT, N, T>::getPosition(*pi);

        for (long j = 0; j < num_newton_iters; ++j)
        {
          // We're trying to minimize (Q(t) - P)^2, i.e. finding the root of (Q(t) - P).Q'(t) = 0. The corresponding
          // Newton-Raphson step is t <-- t - f(t)/f'(t), where f(t) = (Q(t) - P).Q'(t). Differentiating, we get
          // f'(t) = Q'(t).Q'(t) + (Q(t) - P).Q''(t)

          double t = static_cast<T>(u[i]);
          VectorT q   =  eval(t, 0);
          VectorT q1  =  eval(t, 1);
          VectorT q2  =  eval(t, 2);

          double numer = static_cast<double>((q - p).dot(q1));
          double denom = static_cast<double>(q1.dot(q1) + (q - p).dot(q2));
          if (std::fabs(denom) >= Math::eps<double>())
            u[i] -= numer / denom;
        }
      }

      return true;
    }

}; // class BezierN

} // namespace Thea

#endif
