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
#include "../MatVec.hpp"
#include "../Plane3.hpp"
#include "CentroidN.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"
#include <Eigen/Eigenvalues>

namespace Thea {
namespace Algorithms {

/** Fitting linear models to 3D data by minimizing sum-of-squared-errors. */
template <typename T, typename Enable = void>
class /* THEA_API */ LinearLeastSquares3
{
  public:
    /**
     * Linear least-squares fitting of a line to a set of 3D objects. InputIterator must dereference to type T or pointer-to-T.
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
    static double fitLine(InputIterator begin, InputIterator end, Line3 & line, Vector3 * centroid = nullptr);

    /**
     * Linear least-squares fitting of a plane to a set of 3D objects. InputIterator must dereference to type T or pointer-to-T.
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
    static double fitPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 * centroid = nullptr);

}; // class LinearLeastSquares3

// Fitting linear models to sets of objects that map to single points in 3-space.
template <typename T>
class LinearLeastSquares3<T, typename std::enable_if< IsNonReferencedPointN<T, 3>::value >::type>
{
  public:
    template <typename InputIterator>
    static double fitLine(InputIterator begin, InputIterator end, Line3 & line, Vector3 * centroid = nullptr)
    {
      Vector3d center;
      Matrix3d cov = covMatrix(begin, end, center);

      double sum = 0;
      for (auto iter = makeRefIterator(begin); iter != makeRefIterator(end); ++iter)
      {
        Vector3d diff = PointTraitsN<T, 3>::getPosition(*iter).template cast<double>() - center;
        sum += (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      }

      Matrix3d m = sum * Matrix3d::Identity() - cov;
      Eigen::SelfAdjointEigenSolver<Matrix3d> eigensolver;
      if (eigensolver.computeDirect(m).info() != Eigen::Success)
        throw Error("LinearLeastSquares3: Could not eigensolve matrix");

      // Eigenvalues are in increasing order
      line = Line3::fromPointAndDirection(center.cast<Real>(), eigensolver.eigenvectors().col(0).cast<Real>());
      if (centroid) *centroid = center.cast<Real>();
      return eigensolver.eigenvalues()[0];
    }

    template <typename InputIterator>
    static double fitPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 * centroid = nullptr)
    {
      Vector3d center;
      Matrix3d cov = covMatrix(begin, end, center);

      Eigen::SelfAdjointEigenSolver<Matrix3d> eigensolver;
      if (eigensolver.computeDirect(cov).info() != Eigen::Success)
        throw Error("LinearLeastSquares3: Could not eigensolve covariance matrix");

      // Eigenvalues are non-negative (covariance matrix is non-negative definite) and in increasing order
      plane = Plane3::fromPointAndNormal(center.cast<Real>(), eigensolver.eigenvectors().col(0).cast<Real>());
      if (centroid) *centroid = center.cast<Real>();
      return eigensolver.eigenvalues()[0];
    }

  private:
    /** Compute the covariance matrix between the coordinates of a set of 3D points. */
    template <typename InputIterator>
    static Matrix3d covMatrix(InputIterator begin, InputIterator end, Vector3d & centroid)
    {
      centroid = CentroidN<T, 3>::compute(begin, end).template cast<double>();

      Matrix3d m = Matrix3d::Zero();
      for (auto iter = makeRefIterator(begin); iter != makeRefIterator(end); ++iter)
      {
        Vector3d diff = PointTraitsN<T, 3>::getPosition(*iter).template cast<double>() - centroid;
        m += diff * diff.transpose();  // outer product
      }

      return m;
    }

}; // class LinearLeastSquares3<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
