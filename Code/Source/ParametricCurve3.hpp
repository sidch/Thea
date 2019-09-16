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
     * derivative (dimensions > 3 require the third derivative), and the return value is zero if hasDeriv(2) returns false.
     */
    VectorT getBinormal(T const & t) const
    {
      if (!this->hasDeriv(2)) return VectorT::Zero();

      VectorT d1 = this->eval(t, 1);
      T d1_sqlen = d1.squaredNorm();
      if (Math::fuzzyEq(d1_sqlen, static_cast<T>(0)))
        return VectorT::Zero();

      VectorT d2 = this->eval(t, 2);
      d2 -= (d2.dot(d1) / d1_sqlen) * d1;  // sqrt in normalizing d1 avoided because of repeated d1

      return d1.cross(d2).normalized();
    }

}; // class ParametricCurveN

template <typename T> ParametricCurveN<3, T>::~ParametricCurveN() {}  // pure virtual destructor should have body

} // namespace Thea

#endif
