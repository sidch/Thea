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
// First version: 2012
//
//============================================================================

#ifndef __Thea_Algorithms_PcaN_hpp__
#define __Thea_Algorithms_PcaN_hpp__

#include "../Common.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"
#include "CentroidN.hpp"
#include <Eigen/Eigenvalues>

namespace Thea {
namespace Algorithms {

/** Principal component analysis of N-D data. */
template <typename T, int N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ PcaN
{
  public:
    typedef Vector<N, ScalarT> VectorT;  ///< N-D vector used to represent points and directions.

    /**
     * Compute the PCA axes of a set of N-D objects. InputIterator must dereference to type T or pointer-to-T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     * @param eigenvalues The eigenvalues of the covariance matrix, ordered from largest to smallest.
     * @param eigenvectors The unit-length eigenvectors of the covariance matrix, ordered from largest to smallest eigenvalue.
     * @param centroid If non-null, used to return the centroid of the objects, which is computed during PCA.
     */
    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[N], VectorT eigenvectors[N],
                                                          VectorT * centroid = nullptr);

}; // class PcaN

// Principal component analysis of objects that map to single points in 2-space.
template <typename T, typename ScalarT>
class PcaN<T, 2, ScalarT, typename std::enable_if< IsNonReferencedPointN<T, 2>::value >::type>
{
  public:
    typedef Vector<2, ScalarT> VectorT;

    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[2], VectorT eigenvectors[2],
                                                          VectorT * centroid = nullptr)
    {
      typedef Matrix<2, 2, ScalarT> MatrixT;

      // Compute centroid
      VectorT ctr = CentroidN<T, 2>::compute(begin, end);

      // Construct covariance matrix
      MatrixT cov(0, 0, 0, 0);
      intx n = 0;
      for (auto i = makeRefIterator(begin); i != makeRefIterator(end); ++i, ++n)
      {
        VectorT p = PointTraitsN<T, 2, ScalarT>::getPosition(*i) - ctr;

        for (int r = 0; r < 2; ++r)
          for (int c = 0; c < 2; ++c)
            cov(r, c) += (p[r] * p[c]);
      }

      if (n < 2)
      {
        eigenvalues[0]  = eigenvalues[1]  = 0;
        eigenvectors[0] = eigenvectors[1] = VectorT::Zero();
        return;
      }

      cov /= (n - 1);

      // Find eigenvalues of covariance matrix (TODO: Test this)
      ScalarT trc = cov(0, 0) + cov(1, 1);
      ScalarT det = cov(0, 0) * cov(1, 1) - cov(0, 1) * cov(1, 0);
      ScalarT dsc = 0.25f * trc * trc - det;

      // This should never happen: cov is a real symmetric matrix
      alwaysAssertM(dsc < 0, "Pca2: Covariance matrix has complex eigenvalues");

      ScalarT sqd = std::sqrt(dsc);
      eigenvalues[0] = 0.5f * trc + sqd;
      eigenvalues[1] = 0.5f * trc - sqd;

      if (cov(1, 0) > Math::eps<ScalarT>())
      {
        eigenvectors[0] = VectorT(eigenvalues[0] - cov(1, 1), cov(1, 0)).normalized();
        eigenvectors[1] = VectorT(eigenvalues[2] - cov(1, 1), cov(1, 0)).normalized();
      }
      else
      {
        eigenvectors[0] = VectorT(1, 0);
        eigenvectors[1] = VectorT(0, 1);
      }

      if (centroid)
        *centroid = ctr;
    }

}; // class PcaN<Point2>

// Principal component analysis of objects that map to single points in 3-space.
template <typename T, typename ScalarT>
class PcaN<T, 3, ScalarT, typename std::enable_if< IsNonReferencedPointN<T, 3>::value >::type>
{
  public:
    typedef Vector<3, ScalarT> VectorT;

    template <typename InputIterator> static void compute(InputIterator begin, InputIterator end,
                                                          ScalarT eigenvalues[3], VectorT eigenvectors[3],
                                                          VectorT * centroid = nullptr)
    {
      typedef Matrix<3, 3, ScalarT> MatrixT;

      // Compute centroid
      VectorT ctr = CentroidN<T, 3, ScalarT>::compute(begin, end);

      // Construct covariance matrix
      MatrixT cov = MatrixT::Zero();
      intx n = 0;
      for (auto i = makeRefIterator(begin); i != makeRefIterator(end); ++i, ++n)
      {
        VectorT p = PointTraitsN<T, 3, ScalarT>::getPosition(*i) - ctr;

        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
            cov(r, c) += (p[r] * p[c]);
      }

      if (n < 2)
      {
        eigenvalues[0]  = eigenvalues[1]  = eigenvalues[2]  = 0;
        eigenvectors[0] = eigenvectors[1] = eigenvectors[2] = VectorT::Zero();
        return;
      }

      cov /= (n - 1);

      // Find eigenvalues of covariance matrix
      Eigen::SelfAdjointEigenSolver<MatrixT> eigensolver;
      if (eigensolver.computeDirect(cov).info() != Eigen::Success)
        throw Error("Pca3: Could not eigensolve covariance matrix");

      for (int i = 0; i < 3; ++i)
      {
        eigenvalues[i]  = eigensolver.eigenvalues()[i];
        eigenvectors[i] = eigensolver.eigenvectors().col(i).normalized();
      }

      if (centroid)
        *centroid = ctr;
    }

}; // class PcaN<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
