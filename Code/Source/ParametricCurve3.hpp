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

#ifndef __Thea_ParametricCurve3_hpp__
#define __Thea_ParametricCurve3_hpp__

#include "Common.hpp"
#include "ParametricCurveN.hpp"

namespace Thea {

/** A parametric curve segment in 3-dimensional space. */
template <typename T>
class /* THEA_API */ ParametricCurveN<3, T> : public Internal::ParametricCurveNBase<3, T>
{
  private:
    typedef Internal::ParametricCurveNBase<3, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    /** Constructor, sets parameter limits. */
    ParametricCurveN(T const & min_param_ = 0, T const & max_param_ = 1) : BaseT(min_param_, max_param_) {}

    /** Destructor. */
    virtual ~ParametricCurveN() = 0;

    /**
     * Get the unit binormal vector (third Frenet vector) to the curve at parameter value \a t. This requires the second
     * derivative (dimensions > 3 require the third derivative), and the return value is zero if hasDeriv(2) returns false (or
     * N < 3, in which case the binormal is undefined).
     */
    VectorT getBinormal(T const & t) const
    {
      if (!this->hasDeriv(2)) return VectorT::Zero();

      VectorT d1 = this->eval(t, 1);
      T d1_sqlen = d1.squaredNorm();
      if (Math::fuzzyEq(d1_sqlen, static_cast<T>(0), Math::square(Math::eps<T>())))
        return VectorT::Zero();

      VectorT d2 = this->eval(t, 2);
      d2 -= (d2.dot(d1) / d1_sqlen) * d1;  // sqrt in normalizing d1 avoided because of repeated d1

      return d1.cross(d2).normalized();
    }

}; // class ParametricCurveN

template <typename T> ParametricCurveN<3, T>::~ParametricCurveN() {}  // pure virtual destructor should have body

} // namespace Thea

#endif
