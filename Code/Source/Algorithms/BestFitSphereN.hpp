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
// First version: 2016
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
template <int N, typename T = Real>
class /* THEA_API */ BestFitSphereN
{
  private:
    typedef Vector<N, double> PointT;  ///< The N-dimensional vector type used to internally represent points.

  public:
    THEA_DECL_SMART_POINTERS(BestFitSphereN)

    typedef BallN <N, T> BallT;    ///< An N-dimensional ball.
    typedef Vector<N, T> VectorT;  ///< An N-dimensional vector.

    /** Default constructor. */
    BestFitSphereN() : ball(VectorT::Zero(), 0), updated(true) {}

    /** Add a point to the set. */
    void addPoint(VectorT const & point) { points.push_back(point.template cast<double>()); updated = false; }

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
        ball = BallT(VectorT::Zero(), 0);
      }
      else if (points.size() == 1)
      {
        ball = BallT(points[0].template cast<T>(), 0);
      }
      else
      {
        typedef Seb::Smallest_enclosing_ball< (unsigned int)N, double, PointT, Array<PointT> > Miniball;

        Miniball mb(const_cast< Array<PointT> &>(points));
        mb.invalidate();  // to schedule recompute; this does not get rid of data
        double r = mb.radius();
        double const * mbc = mb.center_begin();
        VectorT c;
        for (intx i = 0; i < N; ++i)
          c[i] = mbc[i];

        ball = BallT(c, r);
      }

      updated = true;
    }

    Array<PointT> points;
    mutable BallT ball;
    mutable bool updated;

}; // class BestFitSphereN

} // namespace Algorithms
} // namespace Thea

#endif
