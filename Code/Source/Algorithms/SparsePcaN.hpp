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

#ifndef __Thea_Algorithms_SparsePcaN_hpp__
#define __Thea_Algorithms_SparsePcaN_hpp__

#include "../Common.hpp"
#include "CentroidN.hpp"
#include "IteratorModifiers.hpp"
#include "PcaN.hpp"
#include "PointTraitsN.hpp"
#include "../Math.hpp"
#include <algorithm>

namespace Thea {
namespace Algorithms {

/**
 * Sparse Principal Component Analysis -- biases basis vectors to have more zeros. See:
 *
 * C. Chennubhotla and A. Jepson, "Sparse PCA (SPCA): Extracting Multi-Scale Structure from Data", Proc. ICCV, 641--647, 2001.
 */
template <typename T, int N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ SparsePcaN
{
  public:
    typedef Vector<N, ScalarT> VectorT;  ///< N-dimensional vector class used for positions and directions.

    /**
     * Compute the sparse PCA axes of a set of N-dimensional objects. InputIterator must dereference to type T or pointer-to-T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     * @param variances The variances of the data along the sparse PCA axes, sorted from maxium variances to least.
     * @param axes The sparse PCA axes, sorted from maximum variance to least.
     * @param centroid If non-null, used to return the centroid of the objects, which is computed during PCA.
     * @param lambda The bias towards sparsity. A negative value selects a default bias.
     * @param eps When the cost changes less than this between successive iterations, the search is assumed to have converged.
     *   A negative value selects a default threshold.
     * @param max_iters Maximum number of iterations. A negative value selects a default number.
     *
     * @return True if the search converged, false otherwise.
     */
    template <typename InputIterator> static bool compute(InputIterator begin, InputIterator end, ScalarT variances[N],
                                                          VectorT axes[N], VectorT * centroid = nullptr, ScalarT lambda = -1,
                                                          ScalarT eps = -1, intx max_iters = -1);

}; // class SparsePcaN

// Principal component analysis of objects that map to single points in N-space.
template <int N, typename T, typename ScalarT>
class SparsePcaN<T, N, ScalarT, typename std::enable_if< IsNonReferencedPointN<T, N>::value >::type>
{
  public:
    typedef Vector<N, ScalarT> VectorT;

    template <typename InputIterator> static bool compute(InputIterator begin, InputIterator end, ScalarT variances[N],
                                                          VectorT axes[N], VectorT * centroid = nullptr, ScalarT lambda = -1,
                                                          ScalarT eps = -1, intx max_iters = -1)
    {
      if (lambda < 0)
        lambda = (ScalarT)1.0 / ((ScalarT)N * (ScalarT)std::log((double)N));

      if (eps < 0)
        eps = Math::eps<ScalarT>();

      if (max_iters < 0)
        max_iters = 100 * N * N;

      // Compute initial PCA basis
      VectorT ctr;
      if (!centroid)
        centroid = &ctr;

      PcaN<T, N>::compute(begin, end, variances, axes, centroid);

      // Normalize initial variances
      ScalarT normalized_variances[N];
      normalizeVariances<N>(variances, normalized_variances);

      // Compute initial cost
      ScalarT old_cost = cost<N>(normalized_variances, axes, lambda);

      // Compute the lambda for the 2-D sweep
      ScalarT lambda2 = lambda * (ScalarT)(std::log((double)N) / std::log(2.0));

      // Iterate to minimize the cost function
      bool updated = false;
      for (intx iter = 0; iter < max_iters; ++iter)
      {
        bool updated_in_iter = false;

        for (intx i = 0; i < N; ++i)
          for (intx j = i + 1; j < N; ++j)
          {
            VectorT selected_axes[2] = { axes[i], axes[j] };

            ScalarT vars[2] = { variances[i], variances[j] };
            ScalarT nvars[2]; normalizeVariances<2>(vars, nvars);
            ScalarT best_cost = cost<2>(nvars, selected_axes, lambda2);

            static intx const ANG_ITERS = 180;
            for (intx k = 1; k < ANG_ITERS; ++k)
            {
              double ang = Math::pi() * (k / (double)ANG_ITERS);
              double s = std::sin(ang);
              double c = std::cos(ang);
              VectorT rot_axes[2] = { c * selected_axes[0] - s * selected_axes[1],
                                      s * selected_axes[0] + c * selected_axes[1] };

              vars[0] = computeVariance(begin, end, rot_axes[0], *centroid);
              vars[1] = computeVariance(begin, end, rot_axes[1], *centroid);
              normalizeVariances<2>(vars, nvars);

              ScalarT rot_cost = cost<2>(nvars, rot_axes, lambda2);
              if (rot_cost < best_cost)
              {
                axes[i] = rot_axes[0]; axes[j] = rot_axes[1];
                variances[i] = vars[0]; variances[j] = vars[1];

                best_cost = rot_cost;
                updated = updated_in_iter = true;
              }
            }
          }

        if (updated_in_iter)
        {
          normalizeVariances<N>(variances, normalized_variances);
          ScalarT new_cost = cost<N>(normalized_variances, axes, lambda);

          // THEA_CONSOLE << "Cost after iteration " << iter << " = " << new_cost;

          if (std::abs(old_cost - new_cost) < eps)
          {
            sortAxes(variances, axes);
            return true;
          }

          old_cost = new_cost;
        }
        else
        {
          if (updated)
            sortAxes(variances, axes);

          return true;
        }
      }

      if (updated)
        sortAxes(variances, axes);

      THEA_WARNING << "SparsePCA" << N << ": No convergence after " << max_iters
                   << " iterations, returning best result found so far (cost = " << old_cost << ')';
      return false;
    }

  private:
    // Compute variance along a basis direction.
    template <typename InputIterator>
    static ScalarT computeVariance(InputIterator begin, InputIterator end, VectorT const & axis, VectorT const & centroid)
    {
      ScalarT var = 0;
      for (auto pi = makeRefIterator(begin); pi != makeRefIterator(end); ++pi)
        var += Math::square(axis.dot(PointTraitsN<T, N, ScalarT>::getPosition(*pi) - centroid));

      return var;
    }

    // Compute normalized variances.
    template <int M> static void normalizeVariances(ScalarT const * variances, ScalarT * normalized_variances)
    {
      ScalarT sum = 0;
      for (intx i = 0; i < M; ++i)
        sum += variances[i];

      if (sum > Math::eps<ScalarT>())
      {
        for (intx i = 0; i < M; ++i)
          normalized_variances[i] = variances[i] / sum;
      }
      else
      {
        for (intx i = 0; i < M; ++i)
          normalized_variances[i] = 0;
      }
    }

    // Cost function.
    template <int M> static ScalarT cost(ScalarT const * normalized_variances, VectorT axes[M], ScalarT lambda)
    {
      return cost1<M>(normalized_variances) + lambda * cost2<M>(axes);
    }

    // First term of cost function.
    template <int M> static ScalarT cost1(ScalarT const * normalized_variances)
    {
      // \sum_{i = 0}^{N - 1} -d_m \log d_m

      ScalarT c = 0;
      for (intx i = 0; i < M; ++i)
        if (normalized_variances[i] > Math::eps<ScalarT>())
          c += -normalized_variances[i] * std::log(normalized_variances[i]);

      return c;
    }

    // Second term of cost function.
    template <int M> static ScalarT cost2(VectorT axes[M])
    {
      // \sum_{i = 0}^{N - 1} \sum_{j = 0}^{N - 1} -axes_{i, j}^2 \log -axes_{i, j}^2

      ScalarT c = 0;
      for (intx i = 0; i < M; ++i)
        for (intx j = 0; j < N; ++j)
        {
          ScalarT sq = axes[i][j] * axes[i][j];
          if (sq > Math::eps<ScalarT>())
            c += -sq * std::log(sq);
        }

      return c;
    }

    // Sort axes from largest to smallest variance.
    static void sortAxes(ScalarT variances[N], VectorT axes[N])
    {
      for (intx i = 0; i < N; ++i)
        for (intx j = i + 1; j < N; ++j)
        {
          if (variances[i] < variances[j])
          {
            std::swap(variances[i], variances[j]);
            std::swap(axes[i], axes[j]);
          }
        }
    }

}; // class SparsePcaN<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
