//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2016, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_BestFitSphereN_hpp__
#define __Thea_Algorithms_BestFitSphereN_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../BallN.hpp"
#include "miniball/Seb.h"

namespace Thea {
namespace Algorithms {

/** Approximate best-fit sphere. */
template <long N, typename T = Real>
class /* THEA_API */ BestFitSphereN
{
  private:
    typedef VectorN<N, double> PointT;  ///< The N-dimensional vector type used to internally represent points.

  public:
    THEA_DEF_POINTER_TYPES(BestFitSphereN, shared_ptr, weak_ptr)

    typedef BallN   <N, T> BallT;    ///< An N-dimensional ball.
    typedef VectorN <N, T> VectorT;  ///< An N-dimensional vector.

    /** Default constructor. */
    BestFitSphereN() : ball(VectorT::zero(), 0), updated(true) {}

    /** Add a point to the set. */
    void addPoint(VectorT const & point) { points.push_back(point); updated = false; }

    /** Remove all data and (lazily) set the sphere to null. */
    void clear() { points.clear(); updated = false; }

    /** Get the radius of the sphere. */
    T getRadius() const { update(); return ball.getRadius(); }

    /** Get the diameter of the sphere. */
    T getDiameter() const { update(); return ball.getDiameter(); }

    /** Get the center of the sphere. */
    VectorT const & getCenter() const { update(); return ball.getCenter(); }

    /** Get the ball bounded by the sphere. */
    BallT const & getBall() const { update(); return ball; }

  private:
    /** Recompute the best-fit sphere. */
    void update() const
    {
      if (updated)
        return;

      if (points.empty())
      {
        ball = BallT(VectorT::zero(), 0);
      }
      else if (points.size() == 1)
      {
        ball = BallT(points[0], 0);
      }
      else
      {
        typedef Seb::Smallest_enclosing_ball< (unsigned int)N, double, PointT, TheaArray<PointT> > Miniball;

        Miniball mb(const_cast< TheaArray<PointT> &>(points));
        mb.invalidate();  // to schedule recompute; this does not get rid of data
        double r = mb.radius();
        double const * mbc = mb.center_begin();
        VectorT c;
        for (long i = 0; i < N; ++i)
          c[i] = mbc[i];

        ball = BallT(c, r);
      }

      updated = true;
    }

    TheaArray<PointT> points;
    mutable BallT ball;
    mutable bool updated;

}; // class BestFitSphereN

} // namespace Algorithms
} // namespace Thea

#endif
