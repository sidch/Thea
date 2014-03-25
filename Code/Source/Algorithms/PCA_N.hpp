//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_PCA_N_hpp__
#define __Thea_Algorithms_PCA_N_hpp__

#include "../Common.hpp"
#include "../Math.hpp"
#include "../Matrix2.hpp"
#include "../Matrix3.hpp"
#include "../MatrixMN.hpp"
#include "../VectorN.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"
#include "CentroidN.hpp"

namespace Thea {
namespace Algorithms {

/** Principal component analysis of N-D data. */
template <typename T, long N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ PCA_N
{
  public:
    typedef VectorN<N, ScalarT> VectorT;  ///< N-D vector used to represent points and directions.

    /**
     * Compute the PCA axes of a set of N-D objects. InputIterator must dereference to type T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     * @param eigenvalues The eigenvalues of the covariance matrix, ordered from largest to smallest.
     * @param eigenvectors The unit-length eigenvectors of the covariance matrix, ordered from largest to smallest eigenvalue.
     * @param centroid If non-null, used to return the centroid of the objects, which is computed during PCA.
     */
    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[N], VectorT eigenvectors[N],
                                                          VectorT * centroid = NULL);

}; // class PCA_N

// Principal component analysis of N-D data passed as pointers
template <typename T, long N, typename ScalarT>
class /* THEA_API */ PCA_N<T *, N, ScalarT>
{
  public:
    typedef VectorN<N, ScalarT> VectorT;

    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[N], VectorT eigenvectors[N],
                                                          VectorT * centroid = NULL)
    {
      PCA_N<T, N, ScalarT>::compute(PtrToRefIterator<T, InputIterator>(begin), PtrToRefIterator<T, InputIterator>(end),
                                    eigenvalues, eigenvectors, centroid);
    }

}; // class PCA_N<T *>

// Principal component analysis of objects that map to single points in 2-space.
template <typename T, typename ScalarT>
class PCA_N<T, 2, ScalarT, typename boost::enable_if< IsPointN<T, 2> >::type>
{
  public:
    typedef VectorN<2, ScalarT> VectorT;

    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[2], VectorT eigenvectors[2],
                                                          VectorT * centroid = NULL)
    {
      typedef MatrixMN<2, 2, ScalarT> MatrixT;

      // Compute centroid
      VectorT ctr = CentroidN<T, 2>::compute(begin, end);

      // Construct covariance matrix
      MatrixT cov(0, 0, 0, 0);
      long n = 0;
      for (InputIterator i = begin; i != end; ++i, ++n)
      {
        VectorT p = PointTraitsN<T, 2, ScalarT>::getPosition(*i) - ctr;

        for (int r = 0; r < 2; ++r)
          for (int c = 0; c < 2; ++c)
            cov(r, c) += (p[r] * p[c]);
      }

      if (n < 2)
      {
        eigenvalues[0]  = eigenvalues[1]  = 0;
        eigenvectors[0] = eigenvectors[1] = VectorT::zero();
        return;
      }

      cov /= (n - 1);

      // Find eigenvalues of covariance matrix (TODO: Test this)
      ScalarT trc = cov(0, 0) + cov(1, 1);
      ScalarT det = cov(0, 0) * cov(1, 1) - cov(0, 1) * cov(1, 0);
      ScalarT dsc = 0.25f * trc * trc - det;

      // This should never happen: cov is a real symmetric matrix
      alwaysAssertM(dsc < 0, "PCA2: Covariance matrix has complex eigenvalues");

      ScalarT sqd = std::sqrt(dsc);
      eigenvalues[0] = 0.5f * trc + sqd;
      eigenvalues[1] = 0.5f * trc - sqd;

      if (cov(1, 0) > Math::eps<ScalarT>())
      {
        eigenvectors[0] = VectorT(eigenvalues[0] - cov(1, 1), cov(1, 0)).unit();
        eigenvectors[1] = VectorT(eigenvalues[2] - cov(1, 1), cov(1, 0)).unit();
      }
      else
      {
        eigenvectors[0] = VectorT(1, 0);
        eigenvectors[1] = VectorT(0, 1);
      }

      if (centroid)
        *centroid = ctr;
    }

}; // class PCA_N<Point2>

// Principal component analysis of objects that map to single points in 3-space.
template <typename T, typename ScalarT>
class PCA_N<T, 3, ScalarT, typename boost::enable_if< IsPointN<T, 3> >::type>
{
  public:
    typedef VectorN<3, ScalarT> VectorT;

    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[3], VectorT eigenvectors[3],
                                                          VectorT * centroid = NULL)
    {
      typedef MatrixMN<3, 3, ScalarT> MatrixT;

      // Compute centroid
      VectorT ctr = CentroidN<T, 3, ScalarT>::compute(begin, end);

      // Construct covariance matrix
      MatrixT cov = MatrixT::zero();
      long n = 0;
      for (InputIterator i = begin; i != end; ++i, ++n)
      {
        VectorT p = PointTraitsN<T, 3, ScalarT>::getPosition(*i) - ctr;

        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
            cov(r, c) += (p[r] * p[c]);
      }

      if (n < 2)
      {
        eigenvalues[0]  = eigenvalues[1]  = eigenvalues[2]  = 0;
        eigenvectors[0] = eigenvectors[1] = eigenvectors[2] = VectorT::zero();
        return;
      }

      cov /= (n - 1);

      // Find eigenvalues of covariance matrix
      ScalarT eval[3];
      VectorT evec[3];
      cov.eigenSolveSymmetric(eval, evec);

      // Sort from largest to smallest eigenvalue
      int i0, i1, i2;
      if (eval[0] > eval[1])
      {
        if (eval[0] > eval[2])  // 0, 1, 2 or 0, 2, 1
        {
          i0 = 0;
          i1 = (eval[1] > eval[2] ? 1 : 2);
          i2 = (i1 == 1 ? 2 : 1);
        }
        else  // 2, 0, 1
        {
          i0 = 2;
          i1 = 0;
          i2 = 1;
        }
      }
      else
      {
        if (eval[2] > eval[0])  // 2, 1, 0 or 1, 2, 0
        {
          i0 = (eval[2] > eval[1] ? 2 : 1);
          i1 = (i0 == 1 ? 2 : 1);
          i2 = 0;
        }
        else  // 1, 0, 2
        {
          i0 = 1;
          i1 = 0;
          i2 = 2;
        }
      }

      eigenvalues[0] = eval[i0]; eigenvectors[0] = evec[i0].unit();
      eigenvalues[1] = eval[i1]; eigenvectors[1] = evec[i1].unit();
      eigenvalues[2] = eval[i2]; eigenvectors[2] = evec[i2].unit();

      if (centroid)
        *centroid = ctr;
    }

}; // class PCA_N<Point3>

} // namespace Thea
} // namespace Algorithms

#endif
