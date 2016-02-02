//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Hyperplane3_hpp__
#define __Thea_Hyperplane3_hpp__

#include "Common.hpp"
#include "HyperplaneN.hpp"
#include "Math.hpp"

namespace Thea {

/** A plane (2-flat) in 3-dimensional space. */
template <typename T>
class /* THEA_API */ HyperplaneN<3, T> : public Internal::HyperplaneNBase<3, T>
{
  private:
    typedef Internal::HyperplaneNBase<3, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    static HyperplaneN fromThreePoints(VectorT const & point1, VectorT const & point2, VectorT const & point3)
    {
      HyperplaneN hyperplane;

      hyperplane.normal = (point2 - point1).cross(point3 - point1).unit();
      hyperplane.dist = hyperplane.normal.dot(point1);
      return hyperplane;
    }

    /** Construct a hyperplane given coefficients a, b, c, d of the plane equation a * x + b * y + c * z + d = 0. */
    static HyperplaneN fromEquation(T const & a, T const & b, T const & c, T const & d)
    {
      HyperplaneN hyperplane;

      hyperplane.normal = VectorT(a, b, c);
      T sqlen = hyperplane.normal.squaredLength();
      if (Math::fuzzyEq(sqlen, static_cast<T>(0)))
      {
        hyperplane.normal = VectorT::zero();
        hyperplane.dist = 0;
      }
      else
      {
        T len = std::sqrt(sqlen);
        hyperplane.normal /= len;
        hyperplane.dist = -d / len;
      }

      return hyperplane;
    }

    using BaseT::getEquation;

    /** Get the coefficients a, b, c, d of the plane equation a * x + b * y + c * z + d = 0. */
    void getEquation(T & a, T & b, T & c, T & d) const
    {
      a = this->normal.x();
      b = this->normal.y();
      c = this->normal.z();
      d = -this->dist;
    }

}; // class HyperplaneN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API HyperplaneN<3, Real>;
#endif

/** The default plane class in 3-dimensional real space. */
typedef HyperplaneN<3, Real> Plane3;

} // namespace Thea

#endif
