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
#include "VectorN.hpp"
#include "Algorithms/PointTraitsN.hpp"
#include <algorithm>
#include <iterator>

namespace Thea {

/** A parametric curve segment in N-dimensional space. */
template <long N, typename T>
class /* THEA_API */ ParametricCurveN
{
  public:
    typedef VectorN<N, T> VectorT;  ///< N-dimensional vector type.

    /** Constructor, sets parameter limits. */
    ParametricCurveN(T const & min_param_ = 0, T const & max_param_ = 1)
    : min_param(min_param_), max_param(max_param_) {}

    /** Destructor. */
    virtual ~ParametricCurveN() = 0;

    /** Get the minimum possible parameter value for the segment, which is the value at the beginning of the curve. */
    T const & minParam() const { return min_param; }

    /** Get the maximum possible parameter value for the segment, which is the value at the end of the curve. */
    T const & maxParam() const { return max_param; }

    /** Get the order of the curve, or a negative value if the curve has no well-defined order. */
    virtual long getOrder() const { return -1; }

    /** Get the number of control vectors of the curve, or a negative value if the curve is not defined by control vectors. */
    virtual long numControls() const { return -1; }

    /** Get the point on the curve with parameter value \a t, in the range [minParam(), maxParam()]. */
    VectorT getPoint(T const & t) const
    {
      return eval(t, 0);
    }

    /**
     * Check if the curve's getDeriv() function supports evaluating derivatives of a given order (1 for first derivative, 2 for
     * second, and so on).
     */
    virtual bool hasDeriv(long deriv_order) const = 0;

    /**
     * Get the first, second, or higher derivative (specified by \a deriv_order = 1, 2, ...) of the curve at parameter value
     * \a t, in the range [minParam(), maxParam()].
     */
    VectorT getDeriv(T const & t, long deriv_order = 1) const
    {
      return eval(t, deriv_order);
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
    virtual void getEvenlySpacedPoints(long num_points, VectorT * pts_begin = NULL, T * params_begin = NULL,
                                       long num_arc_samples = -1) const
    {
      alwaysAssertM(num_points >= 2, "ParametricCurveN: At least two evenly-spaced points must be sampled");

      if (num_arc_samples < 0)
      {
        long curve_complexity = std::max(2L, std::max(numControls() - 1, getOrder()));
        num_arc_samples = curve_complexity * 50;
      }

      alwaysAssertM(num_arc_samples >= 2, "ParametricCurveN: At least two samples must be used to approximate arc lengths");

      // Generate num_points samples with uniform parameter spacing
      TheaArray<VectorT> arc_samples((array_size_t)num_arc_samples);
      T arc_scaling = (max_param - min_param) / static_cast<T>(num_arc_samples - 1);
      T t;
      for (long i = 0; i < num_arc_samples; ++i)
      {
        t = i * arc_scaling + min_param;
        arc_samples[(array_size_t)i] = eval(t, 0);
      }

      // Compute the normalized arc length for each sampled point, by a piecewise linear approximation. This is a monotonically
      // increasing sorted array with first element 0 and last element 1.
      TheaArray<double> arclen;
      chordLengthParametrize(arc_samples.begin(), arc_samples.end(), arclen, 0.0, 1.0);

      // Generate a uniform distribution of normalized arc lengths, map each arc length to the corresponding interval in the
      // generated sequence, and use linear interpolation within the interval to find the corresponding curve parameter
      typedef TheaArray<double>::iterator DoubleIterator;
      DoubleIterator last = arclen.begin();
      for (array_size_t i = 0; i < (array_size_t)num_points; ++i)
      {
        if (i == 0)  // short-circuit for first and last points and save two binary searches
          t = min_param;
        else if (i == (array_size_t)num_points - 1)
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
            array_size_t index = (seg_start - arclen.begin());
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
    virtual VectorT eval(T const & t, long deriv_order) const = 0;

    /** Estimate curve parameters for a sequence of points, by accumulating pairwise segment lengths along the sequence. */
    template <typename InputIterator>
    static void chordLengthParametrize(InputIterator begin, InputIterator end, TheaArray<double> & u,
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
        sum_u += static_cast<double>((curr - last).fastLength());
        u.push_back(sum_u);
        last = curr;
      }

      double umax = u[u.size() - 1];
      if (umax <= 0)
        return;

      double scaling = (max_param_ - min_param_) / umax;
      for (array_size_t i = 0; i < u.size(); ++i)
        u[i] = u[i] * scaling + min_param_;
    }

  private:
    T min_param;  ///< The parameter value of the beginning of the curve.
    T max_param;  ///< The parameter value of the end of the curve.

}; // class ParametricCurveN

template <long N, typename T>
ParametricCurveN<N, T>::~ParametricCurveN()
{
  // Pure virtual destructor should have a body
  //
  // http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
}

} // namespace Thea

#endif
