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
#include "Array.hpp"
#include "MatrixWrapper.hpp"
#include "ParametricCurveN.hpp"
#include "Algorithms/FastCopy.hpp"
#include "Algorithms/StdLinearSolver.hpp"
#include "Algorithms/PointTraitsN.hpp"
#include <iterator>

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
template <int N, typename T = Real>
class /* THEA_API */ SplineN : public ParametricCurveN<N, T>
{
  private:
    typedef ParametricCurveN<N, T> BaseT;  ///< Base curve type.

  public:
    THEA_DECL_SMART_POINTERS(SplineN)

    typedef typename BaseT::VectorT VectorT;  ///< N-dimensional vector type.

    /** Constructor, sets parameter limits. */
    SplineN(T const & min_param_ = 0, T const & max_param_ = 1)
    : BaseT(min_param_, max_param_), changed(true) {}

    /**
     * Get a control vector of the curve.
     *
     * @param index The index of the control vector, in the range 0 to numControls() - 1.
     */
    virtual VectorT const & getControl(intx index) const = 0;

    /**
     * Set a control vector of the curve.
     *
     * @param index The index of the control vector to be set, in the range 0 to numControls() - 1.
     * @param pos The new position of the control vector.
     *
     * @note Subclasses that implement this function should remember to call <tt>setChanged(true)</tt>, to trigger a call to
     *   update() to refresh cached data when required.
     */
    virtual void setControl(intx index, VectorT const & pos) = 0;

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
    template < typename InputIterator,
               typename std::enable_if<
                   Algorithms::IsPointN<typename std::iterator_traits<InputIterator>::value_type, N>::value, int >::type = 0 >
    double fitToPoints(InputIterator begin, InputIterator end,
                       T const * initial_params = nullptr,
                       T * final_params = nullptr,
                       bool fix_first_and_last = false,
                       intx max_reparam_iters = -1,
                       intx num_reparam_steps_per_iter = -1)
    {
      if (max_reparam_iters < 0) max_reparam_iters = (initial_params ? 0 : 3);
      if (num_reparam_steps_per_iter < 0) num_reparam_steps_per_iter = 1;  // conservative choice, following Graphics Gems

      // Initialize parameter values for the points
      Array<float64> u;
      if (initial_params)
      {
        T const * up = initial_params;
        for (InputIterator pi = begin; pi != end; ++pi, ++up)
          u.push_back(static_cast<float64>(*up));
      }
      else
        this->chordLengthParametrize(begin, end, u, (float64)this->minParam(), (float64)this->maxParam());

      if ((intx)u.size() < this->numControls())
      {
        THEA_ERROR << "SplineN: Cannot fit curve to fewer data points than control vectors";
        return -1;
      }

      // Fit the curve iteratively, alternating between a linear least squares fit with known parameter values, and
      // re-estimation of parameters via Newton-Raphson
      float64 sqerr = -1;
      while (true)
      {
        float64 e = llsqFit(begin, end, u, fix_first_and_last);
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

      return (double)sqerr;
    }

    /** Get a textual representation of the curve. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[order = " << this->getOrder() << ", param-range = [" << this->minParam() << ", " << this->maxParam()
          << "], ctrl = [";
      for (intx i = 0; i < this->numControls(); ++i)
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
    virtual void getBasisFunctions(float64 t, VectorX<float64> & b) const = 0;

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
    float64 llsqFit(InputIterator begin, InputIterator end, Array<float64> const & u, bool fix_first_and_last)
    {
      using namespace Algorithms;
      typedef typename std::iterator_traits<InputIterator>::value_type PointT;

      // If the sequence is empty, fitting is not possible
      if (begin == end)
      {
        THEA_ERROR << "SplineN: Cannot fit curve to empty sequence of points";
        return -1;
      }

      // Check if fixing first and last positions is even possible. Else, turn this feature off.
      if (fix_first_and_last && !firstAndLastControlsArePositions())
      {
        THEA_WARNING << "SplineN: The beginning and end of the curve are not the first and last control vectors, hence"
                        " they cannot be fixed: this feature will be disabled";
        fix_first_and_last = false;
      }

      intx num_ctrls = (intx)this->numControls();
      intx num_fixed = (fix_first_and_last ? 2 : 0);
      intx fixed_offset = (fix_first_and_last ? 1 : 0);
      intx num_unknown_ctrls = num_ctrls - num_fixed;
      intx num_unknowns = (intx)(N * num_unknown_ctrls);
      intx num_objectives = (intx)(N * std::distance(begin, end));

      VectorX<float64> basis;
      MatrixX<float64> coeffs(num_objectives, num_unknowns);
      VectorX<float64> constants(num_objectives);

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
      coeffs.fill(0.0);
      size_t i = 0;
      intx obj = 0;
      for (InputIterator pi = begin; pi != end; ++pi, ++i)
      {
        getBasisFunctions(u[i], basis);
        Vector<N, float64> d = PointTraitsN<PointT, N, T>::getPosition(*pi).template cast<float64>();
        if (fix_first_and_last)
        {                                                                       // first control vector of curve is fixed at
          d = d - basis[0] * first_pos.template cast<float64>()                 // starting point of sequence, last control
                - basis[basis.size() - 1] * last_pos.template cast<float64>();  // vector is also fixed to last point
        }

        // One scalar objective for each dimension
        for (intx j = 0; j < N; ++j, ++obj)
        {
          intx offset = j * num_unknown_ctrls;
          coeffs.row(obj).segment(offset, basis.size() - num_fixed)
              = basis.segment(fixed_offset, basis.size() - num_fixed);  // assigning col vector to row vector should be ok

          constants[obj] = d[j];
        }
      }

      // Solve the least-squares linear system
      StdLinearSolver llsq(StdLinearSolver::Method::DEFAULT, StdLinearSolver::Constraint::UNCONSTRAINED);
      if (!llsq.solve(&asLvalue(Math::wrapMatrix(coeffs)), constants.data()))
      {
        THEA_ERROR << "SplineN: Could not solve linear least-squares curve fitting problem";
        return -1;
      }

      float64 const * sol = llsq.getSolution();

      // Update the control vectors
      if (fix_first_and_last) setControl(0, first_pos);

      Array<VectorT> new_ctrls(num_unknown_ctrls);
      for (intx i = 0; i < N; ++i)
      {
        intx offset = i * num_unknown_ctrls;
        for (intx j = 0; j < num_unknown_ctrls; ++j)
          new_ctrls[(size_t)j][i] = static_cast<T>(sol[offset + j]);
      }

      for (intx i = 0; i < num_unknown_ctrls; ++i)
        setControl(i + fixed_offset, new_ctrls[(size_t)i]);

      if (fix_first_and_last) setControl(num_ctrls - 1, last_pos);

      float64 err = 0;
      if (llsq.getSquaredError(&err))
        return err;

      VectorXConstMap<float64> sol_vec(sol, (intx)num_unknowns);
      return (coeffs * sol_vec - constants).squaredNorm();
    }

    /**
     * Optimize each parameter value to bring the curve closer to the corresponding target point, using Newton-Raphson steps.
     *
     * See "An Algorithm for Automatically Fitting Digitized Curves", Philip J. Schneider, <i>%Graphics Gems</i>, 1990.
     */
    template <typename InputIterator>
    bool refineParameters(InputIterator begin, InputIterator end, Array<float64> & u, intx num_newton_iters)
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

        for (intx j = 0; j < num_newton_iters; ++j)
        {
          // We're trying to minimize (Q(t) - P)^2, i.e. finding the root of (Q(t) - P).Q'(t) = 0. The corresponding
          // Newton-Raphson step is t <-- t - f(t)/f'(t), where f(t) = (Q(t) - P).Q'(t). Differentiating, we get
          // f'(t) = Q'(t).Q'(t) + (Q(t) - P).Q''(t)

          float64 t   = static_cast<T>(u[i]);
          VectorT q   =  this->eval(t, 0);
          VectorT q1  =  this->eval(t, 1);
          VectorT q2  =  this->eval(t, 2);

          float64 numer = static_cast<float64>((q - p).dot(q1));
          float64 denom = static_cast<float64>(q1.dot(q1) + (q - p).dot(q2));
          if (std::fabs(denom) >= Math::eps<float64>())
            u[i] = Math::clamp(u[i] - numer / denom, (float64)this->minParam(), (float64)this->maxParam());
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
