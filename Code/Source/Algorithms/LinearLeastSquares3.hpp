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

#ifndef __Thea_Algorithms_LinearLeastSquares3_hpp__
#define __Thea_Algorithms_LinearLeastSquares3_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Line3.hpp"
#include "../MatrixMN.hpp"
#include "../Plane3.hpp"
#include "../VectorN.hpp"
#include "CentroidN.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"

namespace Thea {
namespace Algorithms {

/** Fitting linear models to 3D data by minimizing sum-of-squared-errors. */
template <typename T, typename Enable = void>
class /* THEA_API */ LinearLeastSquares3
{
  public:
    /**
     * Linear least-squares fitting of a line to a set of 3D objects. InputIterator must dereference to type T.
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
    static double fitLine(InputIterator begin, InputIterator end, Line3 & line, Vector3 * centroid = NULL);

    /**
     * Linear least-squares fitting of a plane to a set of 3D objects. InputIterator must dereference to type T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     * @param plane Used to return the best-fit plane.
     * @param centroid If non-null, used to return the centroid of the objects, which is computed in the process of finding
     *   the best-fit plane.
     *
     * @return The fitting quality: 0 (worst) to 1 (perfect).
     */
    template <typename InputIterator>
    static double fitPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 * centroid = NULL);

}; // class LinearLeastSquares3

// Fitting linear models to 3D data passed as pointers.
template <typename T>
class /* THEA_API */ LinearLeastSquares3<T *>
{
  public:
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line3 & line, Vector3 * centroid = NULL)
    {
      return LinearLeastSquares3<T>::fitLine(PtrToRefIterator<T, InputIterator>(begin),
                                             PtrToRefIterator<T, InputIterator>(end), line, centroid);
    }

    template <typename InputIterator>
    static double fitPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 * centroid = NULL)
    {
      return LinearLeastSquares3<T>::fitPlane(PtrToRefIterator<T, InputIterator>(begin),
                                              PtrToRefIterator<T, InputIterator>(end), plane, centroid);
    }

}; // class LinearLeastSquares3

// Fitting linear models to sets of objects that map to single points in 3-space.
template <typename T>
class LinearLeastSquares3<T, typename boost::enable_if< IsNonReferencedPointN<T, 3> >::type>
{
  private:
    typedef VectorN<3, double> DVec3;
    typedef MatrixMN<3, 3, double> DMat3;

  public:
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line3 & line, Vector3 * centroid = NULL)
    {
      DVec3 center;
      DMat3 cov = covMatrix(begin, end, center);

      double sum = 0;
      for (InputIterator iter = begin; iter != end; ++iter)
      {
        DVec3 diff = DVec3(PointTraitsN<T, 3>::getPosition(*iter)) - center;
        sum += (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      }

      DMat3 m = sum * DMat3::identity() - cov;
      DVec3 eigenvector;
      double eigenvalue = eigenSolveSmallest(m, eigenvector);

      line = Line3::fromPointAndDirection(Vector3(center), Vector3(eigenvector));

      if (centroid) *centroid = center;
      return eigenvalue;
    }

    template <typename InputIterator>
    static double fitPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 * centroid = NULL)
    {
      DVec3 center;
      DMat3 cov = covMatrix(begin, end, center);

      DVec3 eigenvector;
      double eigenvalue = eigenSolveSmallest(cov, eigenvector);

      plane = Plane3::fromPointAndNormal(Vector3(center), Vector3(eigenvector));

      if (centroid) *centroid = center;
      return eigenvalue;
    }

  private:
    /** Compute the covariance matrix between the coordinates of a set of 3D points. */
    template <typename InputIterator>
    static DMat3 covMatrix(InputIterator begin, InputIterator end, DVec3 & centroid)
    {
      centroid = CentroidN<T, 3>::compute(begin, end);

      DMat3 m = DMat3::zero();
      for (InputIterator iter = begin; iter != end; ++iter)
      {
        DVec3 diff = DVec3(PointTraitsN<T, 3>::getPosition(*iter)) - centroid;

        m(0, 0) += (diff[0] * diff[0]);
        m(0, 1) += (diff[0] * diff[1]);
        m(0, 2) += (diff[0] * diff[2]);

        m(1, 1) += (diff[1] * diff[1]);
        m(1, 2) += (diff[1] * diff[2]);

        m(2, 2) += (diff[2] * diff[2]);
      }

      m(1, 0) = m(0, 1);
      m(2, 0) = m(0, 2);
      m(2, 1) = m(1, 2);

      return m;
    }

    /** Compute the smallest eigenvalue (returned) and corresponding eigenvector of a symmetric 3x3 matrix. */
    static double eigenSolveSmallest(DMat3 const & m, DVec3 & eigenvector)
    {
      double eigenvalues[3];
      DVec3 eigenvectors[3];

      int num_eigen = m.eigenSolveSymmetric((double *)eigenvalues, (DVec3 *)eigenvectors);
      if (num_eigen <= 0)
        throw Error("LinearLeastSquares2: Could not eigensolve matrix");

      int min_index = 0;
      for (int i = 1; i < num_eigen; ++i)
        if (eigenvalues[i] < eigenvalues[min_index])
          min_index = i;

      eigenvector = eigenvectors[min_index];
      return eigenvalues[min_index];
    }

}; // class LinearLeastSquares3<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
