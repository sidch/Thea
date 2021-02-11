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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_Icp3_hpp__
#define __Thea_Algorithms_Icp3_hpp__

#include "../Common.hpp"
#include "KdTreeN.hpp"
#include "MetricL2.hpp"
#include "PointTraitsN.hpp"
#include "../AffineTransformN.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../HyperplaneN.hpp"
#include <Eigen/SVD>

namespace Thea {
namespace Algorithms {

/** Align two sets of points in 3D using the Iterative Closest Point (ICP) algorithm. */
template <typename ScalarT = Real>
class Icp3
{
  private:
    typedef Vector<3, ScalarT>       VectorT;  ///< 3D vector.
    typedef Matrix<3, 3, ScalarT>    MatrixT;  ///< 3x3 matrix.
    typedef HyperplaneN<3, ScalarT>  PlaneT;   ///< Plane in 3-space.

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
     * @param fractional_error_threshold_ The maximum fractional change in error to determine convergence (negative for
     *   default).
     * @param min_iterations_ The minimum number of iterations (negative for default).
     * @param max_iterations_ The maximum number of iterations (negative for default).
     * @param verbose_ If true, print extra debugging information.
     */
    Icp3(ScalarT fractional_error_threshold_ = -1, intx min_iterations_ = -1, intx max_iterations_ = -1, bool verbose_ = false)
    : fractional_error_threshold(fractional_error_threshold_), min_iterations(min_iterations_), max_iterations(max_iterations_),
      has_up(false), verbose(verbose_)
    {
      if (fractional_error_threshold < 0) fractional_error_threshold = 0.0001;
      if (min_iterations < 0) min_iterations = 3;
      if (max_iterations < 0) max_iterations = 10;
    }

    /** Set the up vector. Subsequent alignments will only rotate around the up vector. */
    void setUpVector(VectorT const & up_) { up = up_.normalized(); has_up = true; }

    /** Check if the up vector has been set. */
    bool hasUpVector() const { return has_up; }

    /**
     * Get the up vector, if it has been set.
     *
     * @see hasUpVector();
     */
    VectorT const & getUpVector() const { return up; }

    /** Clear the up vector. Subsequent alignments will be unconstrained. */
    void clearUpVector() { has_up = false; }

    /** Find the transform that best aligns the point set \a from to the point set \a to. */
    template <typename FromT, typename ToT>
    AffineTransformT align(intx from_num_pts, FromT const * from, intx to_num_pts, ToT const * to, ScalarT * error = nullptr)
                     const
    {
      return align(from_num_pts, from, (DefaultWeightFunc<FromT> *)nullptr, to_num_pts, to, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, with a per-point weight assigned to the
     * cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToT, typename FromWeightFuncT>
    AffineTransformT align(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           intx to_num_pts, ToT const * to, ScalarT * error = nullptr) const
    {
      if (from_num_pts <= 0 || to_num_pts <= 0)
        return AffineTransformT::identity();

      KdTreeN<ToT, 3, ScalarT> to_kdtree(to, to + to_num_pts);
      to_kdtree.enableNearestNeighborAcceleration();

      return align(from_num_pts, from, from_weight_func, to_kdtree, error);
    }

    /** Find the transform that best aligns the point set \a from to the proximity query structure \a to. */
    template <typename FromT, typename ToProximityQueryStructureT>
    AffineTransformT align(intx from_num_pts, FromT const * from, ToProximityQueryStructureT const & to,
                           ScalarT * error = nullptr) const
    {
      return align(from_num_pts, from, (DefaultWeightFunc<FromT> *)nullptr, to, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT align(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           ToProximityQueryStructureT const & to, ScalarT * error = nullptr) const
    {
      return align(from_num_pts, from, from_weight_func, nullptr, to, nullptr, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, assuming both sets have known symmetry
     * planes. The computed alignment will ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToT>
    AffineTransformT alignSymmetric(intx from_num_pts, FromT const * from, PlaneT const & from_symmetry_plane,
                                    intx to_num_pts, ToT const * to, PlaneT const & to_symmetry_plane,
                                    ScalarT * error = nullptr) const
    {
      return alignSymmetric(from_num_pts, from, (DefaultWeightFunc<FromT> *)nullptr, from_symmetry_plane,
                            to_num_pts, to, to_symmetry_plane, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the point set \a to, with a per-point weight assigned to the
     * cost of aligning each element of \a from, and assuming both sets have known symmetry planes. The computed alignment will
     * ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToT, typename FromWeightFuncT>
    AffineTransformT alignSymmetric(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                    PlaneT const & from_symmetry_plane,
                                    intx to_num_pts, ToT const * to, PlaneT const & to_symmetry_plane,
                                    ScalarT * error = nullptr) const
    {
      if (from_num_pts <= 0 || to_num_pts <= 0)
        return AffineTransformT::identity();

      KdTreeN<ToT, 3, ScalarT> to_kdtree(to, to + to_num_pts);
      to_kdtree.enableNearestNeighborAcceleration();

      return alignSymmetric(from_num_pts, from, from_weight_func, from_symmetry_plane, to_kdtree, to_symmetry_plane, error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, assuming both \a from
     * and \a to have known symmetry planes. The computed alignment will ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToProximityQueryStructureT>
    AffineTransformT alignSymmetric(intx from_num_pts, FromT const * from, PlaneT const & from_symmetry_plane,
                                    ToProximityQueryStructureT const & to, PlaneT const & to_symmetry_plane,
                                    ScalarT * error = nullptr) const
    {
      return alignSymmetric(from_num_pts, from, (DefaultWeightFunc<FromT> *)nullptr, from_symmetry_plane, to, to_symmetry_plane,
                            error);
    }

    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from, and assuming both \a from and \a to have known symmetry planes.
     * The computed alignment will ensure the symmetry planes coincide.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT alignSymmetric(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                    PlaneT const & from_symmetry_plane,
                                    ToProximityQueryStructureT const & to, PlaneT const & to_symmetry_plane,
                                    ScalarT * error = nullptr) const
    {
      return align(from_num_pts, from, from_weight_func, &from_symmetry_plane, to, &to_symmetry_plane, error);
    }

  private:
    /**
     * Find the transform that best aligns the point set \a from to the proximity query structure \a to, with a per-point weight
     * assigned to the cost of aligning each element of \a from.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT align(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                           PlaneT const * from_symmetry_plane, ToProximityQueryStructureT const & to,
                           PlaneT const * to_symmetry_plane, ScalarT * error = nullptr) const
    {
      if (verbose)
        THEA_CONSOLE << "Icp3(fractional_error_change = " << fractional_error_threshold
                     <<    ", max_iterations = " << max_iterations << ')';

      if (from_num_pts <= 0)
      {
        if (error) *error = 0;
        return AffineTransformT::identity();
      }

      Array<VectorT> from_points((size_t)from_num_pts);
      Array<VectorT> to_points((size_t)from_num_pts);
      for (size_t i = 0; i < from_points.size(); ++i)
        from_points[i] = PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]);

      AffineTransformT old_tr = AffineTransformT::identity();
      ScalarT old_error = measureError(old_tr, from_num_pts, &from_points[0], from_weight_func, to, &to_points[0]);
      if (verbose) THEA_CONSOLE << "[Icp3] Initial error: " << old_error;

      if (old_error <= std::numeric_limits<ScalarT>::min())
      {
        if (error) *error = old_error;
        return old_tr;
      }

      ScalarT new_error = old_error;
      AffineTransformT new_tr = old_tr;
      for (intx i = 0; i < max_iterations; ++i)
      {
        // Align using the point mapping created by the last call to measureError()
        AffineTransformT inc_tr = alignOneStep(from_num_pts, &from_points[0], from_weight_func, from_symmetry_plane,
                                               &to_points[0], to_symmetry_plane);

        // Update the overall transform
        new_tr = inc_tr * old_tr;

        // Compute the new error and the new mapping between points
        if (i < max_iterations - 1 || error)
          new_error = measureError(inc_tr, from_num_pts, &from_points[0], from_weight_func, to, &to_points[0]);

        if (i < max_iterations - 1)
        {
          ScalarT frac_change = (old_error - new_error) / old_error;

          if (verbose)
            THEA_CONSOLE << "[Icp3] Iteration " << i << " error: " << new_error << " (fractional change: " << frac_change << ')';

          if (i >= min_iterations && frac_change < fractional_error_threshold)
          {
            if (frac_change > 0)  // we improved slightly
            {
              if (error) *error = new_error;
              return new_tr;
            }
            else  // the previous alignment was better
            {
              if (error) *error = old_error;
              return old_tr;
            }
          }

          for (size_t j = 0; j < from_points.size(); ++j)  // transform original points to prevent drift
            from_points[j] = new_tr * PointTraitsN<FromT, 3, ScalarT>::getPosition(from[j]);
        }
        else
        {
          if (verbose)
            THEA_CONSOLE << "[Icp3] Iteration " << i << " error: " << new_error;
        }

        old_error = new_error;
        old_tr = new_tr;
      }

      if (error) *error = new_error;
      return new_tr;
    }

    /**
     * Align a point set to a proximity query structure, in one ICP step, after finding the nearest neighbor of each point in
     * the first set in the structure.
     */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    AffineTransformT alignOneStep(intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                  PlaneT const * from_sym_plane, ToProximityQueryStructureT const & to,
                                  PlaneT const * to_sym_plane, VectorT * to_points) const
    {
      for (intx i = 0; i < from_num_pts; ++i)
      {
        intx index = to.template closestElement<MetricL2>(PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]), -1, nullptr,
                                                          &to_points[i]);
        if (index < 0)
          throw Error(format("Icp3: Couldn't get nearest neighbor of source point %ld", i));
      }

      return alignOneStep(from_num_pts, from, from_weight_func, from_sym_plane, to_points, to_sym_plane);
    }

    /** Align one point set to another, in one ICP step, with a known bijective mapping between the points. */
    template <typename FromT, typename FromWeightFuncT>
    AffineTransformT alignOneStep(intx num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                  PlaneT const * from_sym_plane, VectorT const * to, PlaneT const * to_sym_plane) const
    {
      // When both point sets have symmetry planes, we'll align the projections. This does *NOT* minimize the correct error
      // function, but provides a simple approximation to the true solution which is ok for now.
      bool use_symmetry = from_sym_plane && to_sym_plane;

      // Find centroids of the two sets
      VectorT p_mean = VectorT::Zero();
      VectorT q_mean = VectorT::Zero();
      ScalarT sum_weights = 0;
      for (intx i = 0; i < num_pts; ++i)
      {
        VectorT p = PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]);
        VectorT q = to[i];
        if (use_symmetry)
        {
          p = from_sym_plane->closestPoint(p);
          q = to_sym_plane->closestPoint(q);
        }

        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getTranslationWeight(from[i]);
          p_mean += weight * p;
          q_mean += weight * q;
          sum_weights += weight;
        }
        else
        {
          p_mean += p;
          q_mean += q;
          sum_weights += 1;
        }
      }

      p_mean /= sum_weights;
      q_mean /= sum_weights;

      // Find covariance matrix
      MatrixT cov; cov.setZero();
      for (intx i = 0; i < num_pts; ++i)
      {
        VectorT p = PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]);
        VectorT q = to[i];
        if (use_symmetry)
        {
          p = from_sym_plane->closestPoint(p);
          q = to_sym_plane->closestPoint(q);
        }

        VectorT dp = p - p_mean;
        VectorT dq = q - q_mean;

        if (has_up)  // remove the component in the up direction
        {
          dp = dp - (dp.dot(up) * up);
          dq = dq - (dq.dot(up) * up);
        }

        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getRotationWeight(from[i]);
          dp *= weight;
          dq *= weight;
        }

        cov += (dp * dq.transpose());  // outer product
      }

      Eigen::JacobiSVD<MatrixT> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
      MatrixT rot = svd.matrixU() * svd.matrixV().transpose();

      if (use_symmetry && rot.determinant() < 0)  // matching two planar projections can cause flips in the symmetry plane
      {
        // Generate the transformation that flips in the (origin-centered) target symmetry plane
        VectorT n = to_sym_plane->getNormal();
        MatrixT nn = n * n.transpose();  // outer product
        MatrixT flip = MatrixT::Identity() - static_cast<ScalarT>(2) * nn;
        rot = flip * rot;
      }

      VectorT trans = q_mean - rot * p_mean;

      return AffineTransformT(rot, trans);
    }

    /** Measure the alignment error of a mapping from each point of \a from to its nearest neighbor in \a to. */
    template <typename FromT, typename ToProximityQueryStructureT, typename FromWeightFuncT>
    static ScalarT measureError(AffineTransformT const & tr,
                                intx from_num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                ToProximityQueryStructureT const & to, VectorT * to_points)
    {
      if (from_num_pts <= 0)
        return 0;

      for (intx i = 0; i < from_num_pts; ++i)
      {
        intx index = to.template closestElement<MetricL2>(PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]), -1, nullptr,
                                                                                                       &to_points[i]);
        if (index < 0)
          throw Error(format("Icp3: Couldn't get nearest neighbor of source point %ld", i));
      }

      return measureError(tr, from_num_pts, from, from_weight_func, to_points);
    }

    /**
     * Measure the alignment error of a mapping from each point of \a from to each point of \a to, with a known bijective
     * mapping between them.
     */
    template <typename FromT, typename FromWeightFuncT>
    static ScalarT measureError(AffineTransformT const & tr,
                                intx num_pts, FromT const * from, FromWeightFuncT * from_weight_func,
                                VectorT const * to)
    {
      if (num_pts <= 0)
        return 0;

      VectorT p_mean = VectorT::Zero();
      VectorT q_mean = VectorT::Zero();

      ScalarT sum_weights = 0;
      for (intx i = 0; i < num_pts; ++i)
      {
        VectorT p = PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]);
        VectorT const & q = to[i];

        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getTranslationWeight(from[i]);
          p_mean += weight * p;
          q_mean += weight * q;
          sum_weights += weight;
        }
        else
        {
          p_mean += p;
          q_mean += q;
          sum_weights += 1;
        }
      }

      p_mean /= sum_weights;
      q_mean /= sum_weights;

      ScalarT err = (tr * p_mean - q_mean).squaredNorm();

      for (intx i = 0; i < num_pts; ++i)
      {
        VectorT p = PointTraitsN<FromT, 3, ScalarT>::getPosition(from[i]);
        VectorT q = to[i];

        VectorT dp = p - p_mean;
        VectorT dq = q - q_mean;
        if (from_weight_func)
        {
          ScalarT weight = (ScalarT)from_weight_func->getRotationWeight(from[i]);
          err += weight * weight * (tr.getLinear() * dp - dq).squaredNorm();
        }
        else
          err += (tr.getLinear() * dp - dq).squaredNorm();
      }

      return err;
    }

    double fractional_error_threshold;  ///< Maximum fractional change in error to determine convergence.
    intx min_iterations;  ///< Minimum number of iterations.
    intx max_iterations;  ///< Maximum number of iterations.
    bool has_up;
    VectorT up;
    bool verbose;  ///< Print lots of debugging information?

}; // class Icp3

} // namespace Algorithms
} // namespace Thea

#endif
