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

#ifndef __Thea_ParametricCurveN_hpp__
#define __Thea_ParametricCurveN_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "Algorithms/PointTraitsN.hpp"
#include <algorithm>
#include <iterator>

namespace Thea {

// Forward declarations
template <int N, typename T> class ParametricCurveN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for parametric curve segment in N-dimensional space.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ ParametricCurveNBase
{
  public:
    typedef ParametricCurveN<N, T>  ParametricCurveT;  ///< Parametric curve segment in N dimensions.
    typedef Vector<N, T>            VectorT;           ///< N-dimensional vector type.

    THEA_DECL_SMART_POINTERS(ParametricCurveT)

    /** Constructor, sets parameter limits. */
    ParametricCurveNBase(T const & min_param_ = 0, T const & max_param_ = 1)
    : min_param(min_param_), max_param(max_param_) {}

    /** Destructor. */
    virtual ~ParametricCurveNBase() = 0;

    /** Get the minimum possible parameter value for the segment, which is the value at the beginning of the curve. */
    T const & minParam() const { return min_param; }

    /** Get the maximum possible parameter value for the segment, which is the value at the end of the curve. */
    T const & maxParam() const { return max_param; }

    /** Get the order of the curve, or a negative value if the curve has no well-defined order. */
    virtual intx getOrder() const { return -1; }

    /** Get the number of control vectors of the curve, or a negative value if the curve is not defined by control vectors. */
    virtual intx numControls() const { return -1; }

    /** Get the point on the curve with parameter value \a t, in the range [minParam(), maxParam()]. */
    VectorT getPoint(T const & t) const
    {
      return eval(t, 0);
    }

    /**
     * Check if the curve's getDeriv() function supports evaluating derivatives upto and including a given order (1 for first
     * derivative, 2 for second, and so on). Note that if the function returns false for some \a deriv_order, then it must also
     * return false for all higher orders.
     */
    virtual bool hasDeriv(intx deriv_order) const = 0;

    /**
     * Get the first, second, or higher derivative (specified by \a deriv_order = 1, 2, ...) of the curve at parameter value
     * \a t, in the range [minParam(), maxParam()].
     */
    VectorT getDeriv(T const & t, intx deriv_order = 1) const
    {
      return eval(t, deriv_order);
    }

    /**
     * Get the unit tangent vector (first Frenet vector) to the curve at parameter value \a t. This requires the first
     * derivative, and the return value is zero if hasDeriv(1) returns false (or N < 1, in which case the tangent is undefined).
     */
    VectorT getTangent(T const & t) const
    {
      if (N < 1 || !hasDeriv(1)) return VectorT::Zero();

      return eval(t, 1).normalized();
    }

    /**
     * Get the unit normal vector (second Frenet vector) to the curve at parameter value \a t. This requires the second
     * derivative, and the return value is zero if hasDeriv(2) returns false (or N < 2, in which case the normal is undefined).
     */
    VectorT getNormal(T const & t) const
    {
      if (N < 2 || !hasDeriv(2)) return VectorT::Zero();

      VectorT d1 = eval(t, 1);
      VectorT d2 = eval(t, 2);

      T d1_sqlen = d1.squaredNorm();
      if (!Math::fuzzyEq(d1_sqlen, static_cast<T>(0), Math::square(Math::eps<T>())))
        d2 -= (d2.dot(d1) / d1_sqlen) * d1;  // sqrt in normalizing d1 avoided because of repeated d1

      return d2.normalized();
    }

    /**
     * Get the Frenet vector of the curve at parameter value \a t, for a given order. The first Frenet vector is the unit
     * tangent vector returned by getTangent(), the second is the unit normal vector returned by getNormal(), the third the unit
     * binormal, and so on. The return value is zero if hasDeriv(deriv_order) returns false (or \a deriv_order > N, in which
     * case the Frenet vector is undefined).
     *
     * @warning This function is currently strongly susceptible to numerical error, when the <i>unnormalized</i> Frenet vector
     *   of any order (upto \a deriv_order) has length nearly zero, but normalization rescales it to unit length. Use
     *   higher-order Frenet vectors with caution.
     */
    VectorT getFrenetVector(T const & t, intx deriv_order) const
    {
      alwaysAssertM(deriv_order > 0, format("ParametricCurveN: Frenet vector of order %ld does not exist", (long)deriv_order));

      if (deriv_order > N || !hasDeriv(deriv_order)) return VectorT::Zero();

      if (deriv_order == 1) return getTangent(t);
      if (deriv_order == 2) return getNormal(t);
      if (deriv_order == 3) return static_cast<ParametricCurveT>(this)->getBinormal(t);

      // TODO: Is there a faster alternative to this recursion? Probably not very important since deriv_order is unlikely to be
      // more than 2 or 3.
      VectorT d = eval(t, deriv_order);
      VectorT f = d;
      for (intx i = 1; i < deriv_order; ++i)
      {
        VectorT g = getFrenetVector(t, i);
        f -= d.dot(g) * g;
      }

      return f.normalized();
    }

    /**
     * Get a sequence of points roughly evenly spaced by arc length along the curve, with associated parameter values. The
     * beginning and end of the curve (parameter values minParam() and maxParam() respectively) are always included.
     *
     * @param num_points The number of evenly spaced points to be returned (must be at least 2).
     * @param pts_begin If non-null, used to return the point sequence. Must have capacity at least \a num_points.
     * @param params_begin If non-null, used to return the corresponding curve parameter sequence. Must have capacity at least
     *   \a num_points.
     * @param num_arc_samples If non-negative, specifies the number of samples for approximating arc lengths along the curve
     *   (must be at least 2). This value should normally be left at the default setting.
     */
    virtual void getEvenlySpacedPoints(intx num_points, VectorT * pts_begin = nullptr, T * params_begin = nullptr,
                                       intx num_arc_samples = -1) const
    {
      alwaysAssertM(num_points >= 2, "ParametricCurveNBase: At least two evenly-spaced points must be sampled");

      if (num_arc_samples < 0)
      {
        intx curve_complexity = std::max(2L, std::max(numControls() - 1, getOrder()));
        num_arc_samples = curve_complexity * 50;
      }

      alwaysAssertM(num_arc_samples >= 2, "ParametricCurveNBase: At least two samples must be used to approximate arc lengths");

      // Generate num_points samples with uniform parameter spacing
      Array<VectorT> arc_samples((size_t)num_arc_samples);
      T arc_scaling = (max_param - min_param) / static_cast<T>(num_arc_samples - 1);
      T t;
      for (intx i = 0; i < num_arc_samples; ++i)
      {
        t = i * arc_scaling + min_param;
        arc_samples[(size_t)i] = eval(t, 0);
      }

      // Compute the normalized arc length for each sampled point, by a piecewise linear approximation. This is a monotonically
      // increasing sorted array with first element 0 and last element 1.
      Array<double> arclen;
      chordLengthParametrize(arc_samples.begin(), arc_samples.end(), arclen, 0.0, 1.0);

      // Generate a uniform distribution of normalized arc lengths, map each arc length to the corresponding interval in the
      // generated sequence, and use linear interpolation within the interval to find the corresponding curve parameter
      typedef Array<double>::iterator DoubleIterator;
      DoubleIterator last = arclen.begin();
      for (size_t i = 0; i < (size_t)num_points; ++i)
      {
        if (i == 0)  // short-circuit for first and last points and save two binary searches
          t = min_param;
        else if (i == (size_t)num_points - 1)
          t = max_param;
        else
        {
          double s = i / (double)(num_points - 1);  // target arc-length
          DoubleIterator seg_stop = std::upper_bound(last, arclen.end(), s);
          if (seg_stop == arclen.end())  // in degenerate cases
            t = max_param;
          else if (seg_stop != arclen.begin())
          {
            DoubleIterator seg_start = seg_stop - 1;
            double f = (s - *seg_start) / (*seg_stop - *seg_start);  // denom won't be zero because of upper_bound spec
            size_t index = (seg_start - arclen.begin());
            t = static_cast<T>(index + f) * arc_scaling + min_param;
          }

          last = seg_stop;
        }

        if (pts_begin)     pts_begin[i]     =  eval(t, 0);
        if (params_begin)  params_begin[i]  =  t;
      }
    }

    /** Get a textual representation of the curve. */
    virtual std::string toString() const
    {
      std::ostringstream oss;
      oss << "[order = " << getOrder() << ", param-range = [" << min_param << ", " << max_param << ']';
      return oss.str();
    }

  protected:
    /**
     * Evaluate the curve position, or one of its derivatives, at a given parameter value.
     *
     * @param t The curve parameter value, in the range [minParam(), maxParam()].
     * @param deriv_order 0 to get the point position, 1 to get the first derivative, 2 for the second, and so on.
     *
     * @return The desired position or derivative vector.
     */
    virtual VectorT eval(T const & t, intx deriv_order) const = 0;

    /** Estimate curve parameters for a sequence of points, by accumulating pairwise segment lengths along the sequence. */
    template <typename InputIterator>
    static void chordLengthParametrize(InputIterator begin, InputIterator end, Array<double> & u,
                                       double min_param_ = 0, double max_param_ = 1)
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
        sum_u += static_cast<double>((curr - last).norm());
        u.push_back(sum_u);
        last = curr;
      }

      double umax = u[u.size() - 1];
      if (umax <= 0)
        return;

      double scaling = (max_param_ - min_param_) / umax;
      for (size_t i = 0; i < u.size(); ++i)
        u[i] = u[i] * scaling + min_param_;
    }

  private:
    T min_param;  ///< The parameter value of the beginning of the curve.
    T max_param;  ///< The parameter value of the end of the curve.

}; // class ParametricCurveNBase

template <int N, typename T>
ParametricCurveNBase<N, T>::~ParametricCurveNBase()
{
  // Pure virtual destructor should have a body
  // http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
}

} // namespace Internal

/** A parametric curve segment in N-dimensional space. */
template <int N, typename T = Real>
class /* THEA_API */ ParametricCurveN : public Internal::ParametricCurveNBase<N, T>
{
  private:
    typedef Internal::ParametricCurveNBase<N, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    /** Constructor, sets parameter limits. */
    ParametricCurveN(T const & min_param_ = 0, T const & max_param_ = 1) : BaseT(min_param_, max_param_) {}

    /** Destructor. */
    virtual ~ParametricCurveN() = 0;

    /**
     * Get the unit binormal vector (third Frenet vector) to the curve at parameter value \a t. This requires the third
     * derivative, and the return value is zero if hasDeriv(3) returns false (or N < 3, in which case the binormal is
     * undefined).
     *
     * @note This function has an optimized implementation in 3 dimensions, where only the second derivative is required since
     *   the binormal can be computed as the cross-product of the unit tangent and normal vectors.
     *
     * @warning This function is currently strongly susceptible to numerical error, when the <i>unnormalized</i>
     *   tangent/normal/binormal has length nearly zero, but normalization rescales it to unit length. Use with caution.
     */
    VectorT getBinormal(T const & t) const
    {
      if (N < 3 || !this->hasDeriv(3)) return VectorT::Zero();

      VectorT d1 = this->eval(t, 1).normalized();
      VectorT d2 = this->eval(t, 2);
      d2 = (d2 - (d2.dot(d1) * d1)).normalized();

      VectorT d3 = this->eval(t, 3);
      return (d3 - (d3.dot(d1) * d1) - (d3.dot(d2) * d2)).normalized();
    }

}; // class ParametricCurveN

template <int N, typename T> ParametricCurveN<N, T>::~ParametricCurveN() {}  // pure virtual destructor should have body

} // namespace Thea

#include "ParametricCurve3.hpp"

#endif
