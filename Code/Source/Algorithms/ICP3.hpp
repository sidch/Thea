//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_ICP3_hpp__
#define __Thea_Algorithms_ICP3_hpp__

#include "../Common.hpp"
#include "KDTreeN.hpp"
#include "MetricL2.hpp"
#include "PointTraitsN.hpp"
#include "SVD.hpp"
#include "../AffineTransformN.hpp"
#include "../Math.hpp"
#include "../MatrixMN.hpp"
#include "../Plane3.hpp"
#include "../VectorN.hpp"

namespace Thea {
namespace Algorithms {

/** Align two sets of points in 3D using the Iterative Closest Point (ICP) algorithm. */
template <typename ScalarT = Real>
class ICP3
{
  private:
    typedef VectorN<3, ScalarT>      VectorT;  ///< 3D vector.
    typedef MatrixMN<3, 3, ScalarT>  MatrixT;  ///< 3x3 matrix.

    /** The default weight per point. */
    template <typename T> struct DefaultWeightFunc
    {
      double getTranslationWeight(T const & t) const { return 1; }
      double getRotationWeight(T const & t) const { return 1; }
    };

  public:
    typedef AffineTransformN<3, ScalarT> AffineTransformT;  ///< Affine transform in 3 dimensions.

    /**
     * Constructor.
     *
     * @param fractional_error_threshold_ The maximum fractional change in error to determine convergence.
     * @param max_iterations_ The maximum number of iterations.
     */
    ICP3(ScalarT fractional_error_threshold_ = -1, long max_iterations_ = -1, bool verbose_ = false)
    : fractional_error_threshold(fractional_error_threshold_), max_iterations(max_iterations_), verbose(verbose_)
    {
      if (fractional_error_threshold < 0)
        fractional_error_threshold = 0.0001;

      if (max_iterations < 1)
        max_iterations = 10;
    }

    /** Find the transform that best aligns the point set \a from to the point set \a to. */
    template <typename FromT, typename ToT>
    AffineTransformT align(long from_num_pts, FromT const * from, long to_num_pts, ToT const * to, ScalarT * error = NULL) const
    {
      return align(from_num_pts, from, (DefaultWeightFunc<FromT> *)NULL, to_num_pts, to, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, with a per-point weight assigned to the
     * cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToT, typename FromWeightFuncT>
    AffineTransformT align(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           long to_num_pts, ToT const * to, ScalarT * error = NULL) const
    {
      if (from_num_pts <= 0 || to_num_pts <= 0)
        return AffineTransformT::identity();

      KDTreeN<ToT, 3> to_kdtree(to, to + to_num_pts);
      return align(from_num_pts, from, from_weight_func, to_kdtree, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT align(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           ToProximityQueryStructureT const & to, ScalarT * error = NULL) const
    {
      return align(from_num_pts, from, from_weight_func, NULL, to, NULL, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, assuming both sets have known symmetry
     * planes. The computed alignment will ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToT>
    AffineTransformT alignSymmetric(long from_num_pts, FromT const * from, Plane3 const & from_symmetry_plane,
                                    long to_num_pts, ToT const * to, Plane3 const & to_symmetry_plane, ScalarT * error = NULL)
                                    const
    {
      return alignSymmetric(from_num_pts, from, (DefaultWeightFunc<FromT> *)NULL, from_symmetry_plane,
                            to_num_pts, to, to_symmetry_plane, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, with a per-point weight assigned to the
     * cost of aligning each element of \a from, and assuming both sets have known symmetry planes. The computed alignment will
     * ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToT, typename FromWeightFuncT>
    AffineTransformT alignSymmetric(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                    Plane3 const & from_symmetry_plane,
                                    long to_num_pts, ToT const * to, Plane3 const & to_symmetry_plane, ScalarT * error = NULL)
                                    const
    {
      if (from_num_pts <= 0 || to_num_pts <= 0)
        return AffineTransformT::identity();

      KDTreeN<ToT, 3> to_kdtree(to, to + to_num_pts);
      return alignSymmetric(from_num_pts, from, from_weight_func, from_symmetry_plane, to_kdtree, to_symmetry_plane, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from, and assuming both \a from and \a to have known symmetry planes.
     * The computed alignment will ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT alignSymmetric(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                    Plane3 const & from_symmetry_plane,
                                    ToProximityQueryStructureT const & to, Plane3 const & to_symmetry_plane,
                                    ScalarT * error = NULL) const
    {
      return align(from_num_pts, from, from_weight_func, &from_symmetry_plane, to, &to_symmetry_plane, error);
    }

  private:
    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT align(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           Plane3 const * from_symmetry_plane, ToProximityQueryStructureT const & to,
                           Plane3 const * to_symmetry_plane, ScalarT * error = NULL) const
    {
      if (verbose)
        THEA_CONSOLE << "ICP3(fractional_error_change = " << fractional_error_threshold
                     <<    ", max_iterations = " << max_iterations << ')';

      if (from_num_pts <= 0)
      {
        if (error) *error = 0;
        return AffineTransformT::identity();
      }

      TheaArray<Vector3> from_points((array_size_t)from_num_pts);
      TheaArray<Vector3> to_points((array_size_t)from_num_pts);
      for (array_size_t i = 0; i < from_points.size(); ++i)
        from_points[i] = PointTraitsN<FromT, 3>::getPosition(from[i]);

      AffineTransformT tr = AffineTransformT::identity();
      ScalarT old_error = measureError(tr, from_num_pts, &from_points[0], from_weight_func, to, &to_points[0]);
      if (verbose) THEA_CONSOLE << "Initial error: " << old_error;

      if (old_error <= std::numeric_limits<ScalarT>::min())
      {
        if (error) *error = old_error;
        return tr;
      }

      ScalarT new_error = old_error;
      for (long i = 0; i < max_iterations; ++i)
      {
        // Align using the point mapping created by the last call to measureError()
        AffineTransformT inc_tr = alignOneStep(from_num_pts, &from_points[0], from_weight_func, from_symmetry_plane,
                                               &to_points[0], to_symmetry_plane);

        // Update the overall transform
        tr = inc_tr * tr;

        THEA_CONSOLE << "tr = " << tr.toString();
        THEA_CONSOLE << "inc_tr = " << inc_tr.toString();

        // Compute the new error and the new mapping between points
        if (i < max_iterations - 1 || error)
          new_error = measureError(inc_tr, from_num_pts, &from_points[0], from_weight_func, to, &to_points[0]);

        if (i < max_iterations - 1)
        {
          ScalarT frac_change = (old_error - new_error) / old_error;

          if (verbose)
            THEA_CONSOLE << "Iteration " << i << " error: " << new_error << " (fractional change: " << frac_change << ')';

          if (frac_change < fractional_error_threshold)
          {
            if (error) *error = new_error;
            return tr;
          }

          for (array_size_t j = 0; j < from_points.size(); ++j)
            from_points[j] = tr * PointTraitsN<FromT, 3>::getPosition(from[j]);  // transform original points to prevent drift
        }
        else
        {
          if (verbose)
            THEA_CONSOLE << "Iteration " << i << " error: " << new_error;
        }

        old_error = new_error;
      }

      if (error) *error = new_error;
      return tr;
    }

    /**
     * Align a point set to a proximity query structure, in one ICP step, after finding the nearest neighbor of each point in
     * the first set in the structure.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    static AffineTransformT alignOneStep(long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                         Plane3 const * from_sym_plane, ToProximityQueryStructureT const & to,
                                         Plane3 const * to_sym_plane, Vector3 * to_points)
    {
      for (long i = 0; i < from_num_pts; ++i)
      {
        long index = to.template closestElement<MetricL2>(PointTraitsN<FromT, 3>::getPosition(from[i]), -1, NULL,
                                                          &to_points[i]);
        if (index < 0)
          throw Error(format("ICP3: Couldn't get nearest neighbor of source point %ld", i));
      }

      return alignOneStep(from_num_pts, from, from_weight_func, from_sym_plane, to_points, to_sym_plane);
    }

    /** Align one point set to another, in one ICP step, with a known bijective mapping between the points. */
    template <typename FromT, typename FromWeightFuncT>
    static AffineTransformT alignOneStep(long num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                         Plane3 const * from_sym_plane, Vector3 const * to, Plane3 const * to_sym_plane)
    {
      // When both point sets have symmetry planes, we'll align the projections. This does *NOT* minimize the correct error
      // function, but provides a simple approximation to the true solution which is ok for now.
      bool use_symmetry = from_sym_plane && to_sym_plane;

      // Find centroids of the two sets
      VectorT p_mean = VectorT::zero();
      VectorT q_mean = VectorT::zero();
      ScalarT sum_weights = 0;
      for (long i = 0; i < num_pts; ++i)
      {
        Vector3 p = PointTraitsN<FromT, 3>::getPosition(from[i]);
        Vector3 q = to[i];
        if (use_symmetry)
        {
          p = from_sym_plane->closestPoint(p);
          q = to_sym_plane->closestPoint(q);
        }

        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getTranslationWeight(from[i]);
          p_mean += weight * VectorT(p);
          q_mean += weight * VectorT(q);
          sum_weights += weight;
        }
        else
        {
          p_mean += VectorT(p);
          q_mean += VectorT(q);
          sum_weights += 1;
        }
      }

      p_mean /= sum_weights;
      q_mean /= sum_weights;

      // Find covariance matrix
      MatrixT cov(0);
      for (long i = 0; i < num_pts; ++i)
      {
        Vector3 p = PointTraitsN<FromT, 3>::getPosition(from[i]);
        Vector3 q = to[i];
        if (use_symmetry)
        {
          p = from_sym_plane->closestPoint(p);
          q = to_sym_plane->closestPoint(q);
        }

        VectorT dp = VectorT(p) - p_mean;
        VectorT dq = VectorT(q) - q_mean;
        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getRotationWeight(from[i]);
          dp *= weight;
          dq *= weight;
        }

        cov(0, 0) += dp[0] * dq[0];  cov(0, 1) += dp[0] * dq[1];  cov(0, 2) += dp[0] * dq[2];
        cov(1, 0) += dp[1] * dq[0];  cov(1, 1) += dp[1] * dq[1];  cov(1, 2) += dp[1] * dq[2];
        cov(2, 0) += dp[2] * dq[0];  cov(2, 1) += dp[2] * dq[1];  cov(2, 2) += dp[2] * dq[2];
      }

      MatrixT u, v;
      TheaArray<ScalarT> d(3);
      if (!SVD::compute(cov, u, d, v))
        return AffineTransformT::identity();

      MatrixT rot = (u * v).transpose();
      if (use_symmetry && rot.determinant() < 0)  // matching two planar projections can cause flips in the symmetry plane
      {
        // Generate the transformation that flips in the (origin-centered) target symmetry plane
        Vector3 n = to_sym_plane->getNormal();
        MatrixT nn(n[0] * n[0], n[0] * n[1], n[0] * n[2],
                   n[1] * n[0], n[1] * n[1], n[1] * n[2],
                   n[2] * n[0], n[2] * n[1], n[2] * n[2]);
        MatrixT flip = MatrixT::identity() - static_cast<ScalarT>(2) * nn;
        rot = flip * rot;
      }

      VectorT trans = q_mean - rot * p_mean;

      return AffineTransformT(rot, trans);
    }

    /** Measure the alignment error of a mapping from each point of \a from to its nearest neighbor in \a to. */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    static ScalarT measureError(AffineTransformT const & tr,
                                long from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                ToProximityQueryStructureT const & to, Vector3 * to_points)
    {
      if (from_num_pts <= 0)
        return 0;

      for (long i = 0; i < from_num_pts; ++i)
      {
        long index = to.template closestElement<MetricL2>(PointTraitsN<FromT, 3>::getPosition(from[i]), -1, NULL,
                                                                                              &to_points[(array_size_t)i]);
        if (index < 0)
          throw Error(format("ICP3: Couldn't get nearest neighbor of source point %ld", i));
      }

      return measureError(tr, from_num_pts, from, from_weight_func, to_points);
    }

    /**
     * Measure the alignment error of a mapping from each point of \a from to each point of \a to, with a known bijective
     * mapping between them.
     */
    template <typename FromT, typename FromWeightFuncT>
    static ScalarT measureError(AffineTransformT const & tr,
                                long num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                Vector3 const * to)
    {
      if (num_pts <= 0)
        return 0;

      VectorT p_mean = VectorT::zero();
      VectorT q_mean = VectorT::zero();

      ScalarT sum_weights = 0;
      for (long i = 0; i < num_pts; ++i)
      {
        Vector3 p = PointTraitsN<FromT, 3>::getPosition(from[i]);
        Vector3 const & q = to[i];

        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getTranslationWeight(from[i]);
          p_mean += weight * VectorT(p);
          q_mean += weight * VectorT(q);
          sum_weights += weight;
        }
        else
        {
          p_mean += VectorT(p);
          q_mean += VectorT(q);
          sum_weights += 1;
        }
      }

      ScalarT err = (tr * p_mean - q_mean).squaredLength();

      p_mean /= sum_weights;
      q_mean /= sum_weights;

      for (long i = 0; i < num_pts; ++i)
      {
        Vector3 p = PointTraitsN<FromT, 3>::getPosition(from[i]);
        Vector3 q = to[i];

        VectorT dp = VectorT(p) - p_mean;
        VectorT dq = VectorT(q) - q_mean;
        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getRotationWeight(from[i]);
          err += weight * weight * (tr.getLinear() * dp - dq).squaredLength();
        }
        else
          err += (tr.getLinear() * dp - dq).squaredLength();
      }

      return err;
    }

    double fractional_error_threshold;  ///< Maximum fractional change in error to determine convergence.
    long max_iterations;  ///< Maximum number of iterations.
    bool verbose;  ///< Print lots of debugging information?

}; // class ICP3

} // namespace Thea
} // namespace Algorithms

#endif
