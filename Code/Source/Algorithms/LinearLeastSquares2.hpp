//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_LinearLeastSquares2_hpp__
#define __Thea_Algorithms_LinearLeastSquares2_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Line2.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

namespace Thea {
namespace Algorithms {

/** Fitting linear models to 2D data by minimizing sum-of-squared-errors. */
template <typename T, typename Enable = void>
class /* THEA_API */ LinearLeastSquares2
{
  public:
    /**
     * Linear least-squares fitting of a line to a set of 2D objects. InputIterator must dereference to type T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     * @param line Best-fit line.
     *
     * @return The fitting quality: 0 (worst) to 1 (perfect).
     */
    template <typename InputIterator> static double fitLine(InputIterator begin, InputIterator end, Line2 & line);

}; // class LinearLeastSquares2

// Fitting lines to 2D data passed as pointers.
template <typename T>
class /* THEA_API */ LinearLeastSquares2<T *>
{
  public:
    template <typename InputIterator> static double fitLine(InputIterator begin, InputIterator end, Line2 & line)
    {
      LinearLeastSquares2<T>::fitLine(PtrToRefIterator<T const, InputIterator>(begin), PtrToRefIterator<T, InputIterator>(end),
                                      line);
    }

}; // class LinearLeastSquares2<T *>

// Fitting lines to sets of objects that map to single points in 2-space.
template <typename T>
class /* THEA_API */ LinearLeastSquares2<T, typename boost::enable_if< IsPointN<T, 2> >::type>
{
  private:
    typedef CGAL::Cartesian<Real>    CGALKernel;
    typedef CGALKernel::Point_2      CGALPoint;
    typedef CGALKernel::Direction_2  CGALDirection;
    typedef CGALKernel::Line_2       CGALLine;

  public:
    template <typename InputIterator> static double fitLine(InputIterator begin, InputIterator end, Line2 & line)
    {
      TheaArray<CGALPoint> cgal_points;
      toCGALPoints(begin, end, cgal_points);

      CGALLine cgal_line;
      CGALPoint centroid;
      double quality = CGAL::linear_least_squares_fitting_2(cgal_points.begin(), cgal_points.end(), cgal_line, centroid,
                                                            CGAL::Dimension_tag<0>());

      CGALPoint lp = cgal_line.point(0);
      CGALDirection ld = cgal_line.direction();

      line = Line2::fromPointAndDirection(Vector2(lp.x(), lp.y()), Vector2(ld.dx(), ld.dy()));

      return quality;
    }

  private:
    // Convert a set of point objects in 2-space to CGAL points.
    template <typename InputIterator>
    static void toCGALPoints(InputIterator begin, InputIterator end, TheaArray<CGALPoint> & result)
    {
      for (InputIterator i = begin; i != end; ++i)
      {
        Vector2 const & p = PointTraitsN<T, 2>::getPosition(*i);
        result.push_back(CGALPoint(p.x(), p.y()));
      }
    }

}; // class LinearLeastSquares2<Point2>

} // namespace Thea
} // namespace Algorithms

#endif
