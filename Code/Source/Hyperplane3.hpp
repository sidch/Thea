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
// First version: 2011
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

      hyperplane.normal = (point2 - point1).cross(point3 - point1).normalized();
      hyperplane.dist = hyperplane.normal.dot(point1);
      return hyperplane;
    }

    /** Construct a hyperplane given coefficients a, b, c, d of the plane equation a * x + b * y + c * z + d = 0. */
    static HyperplaneN fromEquation(T const & a, T const & b, T const & c, T const & d)
    {
      HyperplaneN hyperplane;

      hyperplane.normal = VectorT(a, b, c);
      T sqlen = hyperplane.normal.squaredNorm();
      if (Math::fuzzyEq(sqlen, static_cast<T>(0)))
      {
        hyperplane.normal = VectorT::Zero();
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
