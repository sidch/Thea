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
#include "../MatrixMN.hpp"
#include "../VectorN.hpp"
#include "CentroidN.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"

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
     * @param line Used to return the best-fit line.
     * @param centroid If non-null, used to return the centroid of the objects, which is computed in the process of finding the
     *   best-fit line.
     *
     * @return The sum of squared fitting errors.
     */
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line2 & line, Vector2 * centroid = NULL);

}; // class LinearLeastSquares2

// Fitting lines to 2D data passed as pointers.
template <typename T>
class /* THEA_API */ LinearLeastSquares2<T *>
{
  public:
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line2 & line, Vector2 * centroid = NULL)
    {
      LinearLeastSquares2<T>::fitLine(PtrToRefIterator<T, InputIterator>(begin),
                                      PtrToRefIterator<T, InputIterator>(end), line, centroid);
    }

}; // class LinearLeastSquares2<T *>

// Fitting lines to sets of objects that map to single points in 2-space.
template <typename T>
class /* THEA_API */ LinearLeastSquares2<T, typename boost::enable_if< IsPointN<T, 2> >::type>
{
  public:
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line2 & line, Vector2 * centroid = NULL)
    {
      typedef VectorN<2, double> DVec2;
      typedef MatrixMN<2, 2, double> DMat2;

      DVec2 center = CentroidN<T, 2>::compute(begin, end);
      DMat2 m = DMat2::zero();
      for (InputIterator iter = begin; iter != end; ++iter)
      {
        DVec2 diff = DVec2(PointTraitsN<T, 2>::getPosition(*iter)) - center;

        m(0, 0) += (diff.y() * diff.y());
        m(0, 1) -= (diff.x() * diff.y());
        m(1, 1) += (diff.x() * diff.x());
      }

      m(1, 0) = m(0, 1);

      double eigenvalues[2];
      DVec2 eigenvectors[2];

      int num_eigen = m.eigenSolve((double *)eigenvalues, (DVec2 *)eigenvectors);
      if (num_eigen == 2)
      {
        if (eigenvalues[1] < eigenvalues[0])
        {
          std::swap(eigenvalues[0], eigenvalues[1]);
          std::swap(eigenvectors[0], eigenvectors[1]);
        }
      }
      else if (num_eigen <= 0)
      {
        throw Error("LinearLeastSquares2: Could not eigensolve matrix");
      }

      line = Line2::fromPointAndDirection(Vector2(center), Vector2(eigenvectors[0]));

      if (centroid) *centroid = center;
      return eigenvalues[0];
    }

}; // class LinearLeastSquares2<Point2>

} // namespace Thea
} // namespace Algorithms

#endif
