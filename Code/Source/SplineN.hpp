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
// Portions of this source file were derived from:
//
// An Algorithm for Automatically Fitting Digitized Curves
// by Philip J. Schneider
// from "Graphics Gems", Academic Press, 1990
//
//============================================================================

#ifndef __Thea_SplineN_hpp__
#define __Thea_SplineN_hpp__

#include "Common.hpp"
#include "ParametricCurveN.hpp"
#include "Algorithms/FastCopy.hpp"
#include "Algorithms/LinearLeastSquares.hpp"
#include "Algorithms/PointTraitsN.hpp"

namespace Thea {

/**
 * A spline curve segment in N-dimensional space. The curve is assumed to be the weighted sum of a set of <i>control
 * vectors</i>, where the weights are given by (typically polynomial or rational) differentiable <i>basis functions</i> of a
 * scalar curve parameter. If the functions are polynomial, the maximum degree of these polynomials is the <i>order</i> of the
 curve.
 *
 * This class contains code for fitting the curve to a sequence of points, extending the code from
 * "An Algorithm for Automatically Fitting Digitized Curves", Philip J. Schneider, <i>%Graphics Gems</i>, Academic Press, 1990.
 */
template <long N, typename T>
class /* THEA_API */ SplineN : public ParametricCurveN<N, T>
{
  private:
    typedef ParametricCurveN<N, T> BaseT;  ///< Base curve type.

  public:
    typedef typename BaseT::VectorT VectorT;  ///< N-dimensional vector type.

    /** Constructor, sets parameter limits. */
    SplineN(T const & min_param_ = 0, T const & max_param_ = 1)
    : BaseT(min_param_, max_param_), changed(true) {}

    /**
     * Get a control vector of the curve.
     *
     * @param index The index of the control vector, in the range 0 to numControls() - 1.
     */
    virtual VectorT const & getControl(long index) const = 0;

    /**
     * Set a control vector of the curve.
     *
     * @param index The index of the control vector to be set, in the range 0 to numControls() - 1.
     * @param pos The new position of the control vector.
     *
     * @note Subclasses that implement this function should remember to call <tt>setChanged(true)</tt>, to trigger a call to
     *   update() to refresh cached data when required.
     */
    virtual void setControl(long index, VectorT const & pos) = 0;

    /**
     * Fit the curve to a sequence of points [begin, end). The algorithm alternates between least-squares fitting (with known
     * parameters) and Newton-Raphson parameter optimization, following %Graphics Gems.
     *
     * @param begin First point in the sequence to be fitted.
     * @param end One past the last point in the sequence to be fitted.
     * @param initial_params The curve parameters of the points, if known in advance. If this argument is not null, no
     *   reparametrization will be done by default, unless \a max_reparam_iters is explicitly set to a positive number.
     * @param fix_first_and_last If true, the first and last control vectors will be set to the positions of the first and last
     *   points, respectively, in the sequence. Note that curves whose first and last control vectors are not precisely the
     *   beginning and end positions of the curve will have this feature automatically disabled, with a warning.
     * @param max_reparam_iters If greater than zero, the parameter values of the points will be re-estimated (at most) this
     *  many times, guided by initial values (\a initial_params) if any. Pass a negative value to pick a suitable default.
     * @param num_reparam_steps_per_iter The number of successive Newton-Raphson steps in each iteration of reparametrization.
     *  Pass a negative value to pick a suitable default.
     * @param final_params If non-null, used to return the final parameter values of the point sequence. Must be pre-allocated
     *  to have at least as many entries as the number of points.
     *
     * @return The non-negative squared fitting error on success, or a negative value on error.
     */
    template <typename InputIterator>
    double fitToPoints(InputIterator begin, InputIterator end,
                       T const * initial_params = NULL,
                       T * final_params = NULL,
                       bool fix_first_and_last = false,
                       long max_reparam_iters = -1,
                       long num_reparam_steps_per_iter = -1,

                       typename boost::enable_if<
                                  Algorithms::IsPointN<typename std::iterator_traits<InputIterator>::value_type, N>
                                                >::type * dummy = NULL)
    {
      if (max_reparam_iters < 0) max_reparam_iters = (initial_params ? 0 : 3);
      if (num_reparam_steps_per_iter < 0) num_reparam_steps_per_iter = 1;  // conservative choice, following Graphics Gems

      // Initialize parameter values for the points
      TheaArray<double> u;
      if (initial_params)
      {
        T const * up = initial_params;
        for (InputIterator pi = begin; pi != end; ++pi, ++up)
          u.push_back(static_cast<double>(*up));
      }
      else
        this->chordLengthParametrize(begin, end, u, (double)this->minParam(), (double)this->maxParam());

      if ((long)u.size() < this->numControls())
      {
        THEA_ERROR << "SplineN: Cannot fit curve to fewer data points than control vectors";
        return -1;
      }

      // Fit the curve iteratively, alternating between a linear least squares fit with known parameter values, and
      // re-estimation of parameters via Newton-Raphson
      double sqerr = -1;
      while (true)
      {
        double e = llsqFit(begin, end, u, fix_first_and_last);
        if (e < 0)  // could not fit, revert to last solution
          break;

        if (sqerr >= 0 && e > sqerr)  // error increased, stop iterating
          break;

        sqerr = e;

        // A bit wasteful to save it every iteration instead of outside the loop, but do it in case llsqFit fails and we have to
        // revert to the last solution
        if (final_params)
          Algorithms::fastCopy(&u[0], &u[0] + u.size(), final_params);

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
      oss << "[order = " << this->getOrder() << ", param-range = [" << this->minParam() << ", " << this->maxParam()
          << "], ctrl = [";
      for (long i = 0; i < this->numControls(); ++i)
      {
        if (i > 0) oss << ", ";
        oss << getControl(i);
      }
      oss << "]]";
      return oss.str();
    }

  protected:
    /**
     * Updated cached data for the curve, marked as invalid after control vectors (for instance) have changed.
     *
     * @note Subclasses implementing this function should call isChanged() at the beginning to check if the curve has actually
     *   been changed, and setChanged(false) at the end to notify that the cached data is now up-to-date. In other words, the
     *   structure of the subclass function should be:
     * \code
     * void update() const
     * {
     *   if (!this->isChanged()) return;
     *   ...
     *   // update cached data
     *   ...
     *   this->setChanged(false);
     * }
     * \endcode
     */
    virtual void update() const = 0;

    /**
     * Get the numControls() basis functions for the curve, evaluated at curve parameter \a t. A point on the curve is the sum
     * of the control vectors, weighted by these basis functions as coefficients. (This need not be the most efficient way to
     * evaluate points on the curve, and hence the implementation of eval() is left to the subclass.)
     */
    virtual void getBasisFunctions(double t, TheaArray<double> & b) const = 0;

    /**
     * Check if the first and last control vectors of the curve can be interpreted as the positions of the beginning and end
     * of the curve respectively.
     */
    virtual bool firstAndLastControlsArePositions() const = 0;

    /**
     * Fit the curve to a sequence of points using linear least-squares. Control vector positions are estimated by minimizing
     * the sum of squared errors between the points and their corresponding curve points with known parameters \a u.
     *
     * @param begin The beginning of the point sequence.
     * @param end One past the end of the point sequence.
     * @param u The curve parameter for each point.
     * @param fix_first_and_last If true, the first and last control vectors will be set to the positions of the first and last
     *   points, respectively, in the sequence. Note that curves whose first and last control vectors are not precisely the
     *   beginning and end positions of the curve will have this feature automatically disabled, with a warning.
     *
     * @return The sum of squared fitting errors, or a negative value on error.
     */
    template <typename InputIterator>
    double llsqFit(InputIterator begin, InputIterator end, TheaArray<double> const & u, bool fix_first_and_last)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      // Check if fixing first and last positions is even possible. Else, turn this feature off.
      if (fix_first_and_last && !firstAndLastControlsArePositions())
      {
        THEA_WARNING << "SplineN: The beginning and end of the curve are not the first and last control vectors, hence"
                        " they cannot be fixed: this feature will be disabled";
        fix_first_and_last = false;
      }

      size_t num_ctrls = (size_t)this->numControls();
      size_t num_fixed = (fix_first_and_last ? 2 : 0);
      size_t fixed_offset = (fix_first_and_last ? 1 : 0);
      size_t num_unknown_ctrls = num_ctrls - num_fixed;
      size_t num_unknowns = (size_t)(N * num_unknown_ctrls);

      LinearLeastSquares llsq((long)num_unknowns);
      TheaArray<double> basis;
      TheaArray<double> coeffs(num_unknowns, 0.0);

      // Cache the first and last points of the sequence
      InputIterator first = begin, last = begin;
      VectorT first_pos, last_pos;

      if (fix_first_and_last)
      {
        for (InputIterator pi = begin; pi != end; ++pi)
          last = pi;

        first_pos  =  PointTraitsN<PointT, N, T>::getPosition(*first);
        last_pos   =  PointTraitsN<PointT, N, T>::getPosition(*last);
      }

      // Try to make each point in the sequence match the point on the curve with the corresponding parameters
      size_t i = 0;
      for (InputIterator pi = begin; pi != end; ++pi, ++i)
      {
        getBasisFunctions(u[i], basis);
        VectorT d = PointTraitsN<PointT, N, T>::getPosition(*pi);
        if (fix_first_and_last)
        {
          d = d - basis[0] * first_pos                 // first control vector of curve is fixed at starting point of sequence
                - basis[basis.size() - 1] * last_pos;  // last control vector is also fixed to last point of sequence
        }

        // One scalar objective for each dimension
        for (long j = 0; j < N; ++j)
        {
          size_t offset = (size_t)j * num_unknown_ctrls;
          fastCopy(&basis[0] + fixed_offset, &basis[0] + basis.size() - fixed_offset, &coeffs[0] + offset);

          llsq.addObjective(&coeffs[0], d[j]);

          // Reset the range
          fill(coeffs.begin() + offset, coeffs.begin() + offset + num_unknown_ctrls, 0.0);
        }
      }

      // Solve the least-squares linear system
      if (!llsq.solve(LinearLeastSquares::Constraint::UNCONSTRAINED))
      {
        THEA_ERROR << "SplineN: Could not solve linear least-squares curve fitting problem";
        return -1;
      }

      TheaArray<double> const & sol = llsq.getSolution();
      alwaysAssertM(sol.size() == num_unknowns,
                    "SplineN: Linear least-squares solution length doesn't match number of unknowns");

      // Update the control vectors
      if (fix_first_and_last) setControl(0, first_pos);

      TheaArray<VectorT> new_ctrls(num_unknown_ctrls);
      for (long i = 0; i < N; ++i)
      {
        size_t offset = (size_t)(i * num_unknown_ctrls);

        for (size_t j = 0; j < num_unknown_ctrls; ++j)
          new_ctrls[j][i] = static_cast<T>(sol[offset + j]);
      }

      for (size_t i = 0; i < num_unknown_ctrls; ++i)
        setControl((long)(i + fixed_offset), new_ctrls[i]);

      if (fix_first_and_last) setControl((long)num_ctrls - 1, last_pos);

      return llsq.squaredError();
    }

    /**
     * Optimize each parameter value to bring the curve closer to the corresponding target point, using Newton-Raphson steps.
     *
     * See "An Algorithm for Automatically Fitting Digitized Curves", Philip J. Schneider, <i>%Graphics Gems</i>, 1990.
     */
    template <typename InputIterator>
    bool refineParameters(InputIterator begin, InputIterator end, TheaArray<double> & u, long num_newton_iters)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      if (!this->hasDeriv(1) || !this->hasDeriv(2))
      {
        THEA_ERROR << "SplineN: Reparametrization requires first and second derivatives";
        return false;
      }

      size_t i = 0;
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
          VectorT q   =  this->eval(t, 0);
          VectorT q1  =  this->eval(t, 1);
          VectorT q2  =  this->eval(t, 2);

          double numer = static_cast<double>((q - p).dot(q1));
          double denom = static_cast<double>(q1.dot(q1) + (q - p).dot(q2));
          if (std::fabs(denom) >= Math::eps<double>())
            u[i] = Math::clamp(u[i] - numer / denom, (double)this->minParam(), (double)this->maxParam());
        }
      }

      return true;
    }

    /** Check if the curve was changed and hence cached data should be recomputed. */
    bool isChanged() const
    {
      return changed;
    }

    /**
     * Mark the curve as having changed or not changed. Setting the value to true will cause update() to refresh cached data.
     * Can be called even from the const update() function to unset the flag.
     */
    void setChanged(bool value = true) const
    {
      changed = value;
    }

  private:
    mutable bool changed;  ///< Was the curve changed without its cached data being updated?

}; // class SplineN

} // namespace Thea

#endif
