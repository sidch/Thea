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

#ifndef __Thea_Algorithms_Symmetry3_hpp__
#define __Thea_Algorithms_Symmetry3_hpp__

#include "../Common.hpp"
#include "IteratorModifiers.hpp"
#include "KDTreeN.hpp"
#include "MetricL2.hpp"
#include "PCA_N.hpp"
#include "PointTraitsN.hpp"
#include <boost/utility/enable_if.hpp>
#include <cmath>

namespace Thea {
namespace Algorithms {

/** Symmetry detection in 3D. */
template <typename T, typename Enable = void>
class /* THEA_API */ Symmetry3
{
  public:
    /**
     * Find a plane of (global) reflective symmetry of a set of objects. Passing a negative value for any of the numeric
     * parameters indicates that a default value should be chosen. A kd-tree will be constructed to answer proximity queries on
     * the set of objects. To use a precomputed query structure instead, use the other form of this function.
     *
     * @param begin The first object.
     * @param end One position beyond the last object.
     * @param plane Used to store the resulting symmetry plane.
     * @param precomputed_centroid The precomputed centroid of the objects. If NULL, the centroid will be computed from the
     *   input data.
     * @param max_rms_error The maximum allowed error, measured as the root of the average squared distance (RMS) of the
     *   reflections of the objects in the candidate plane to their respective closest neighbors.
     */
    template <typename InputIterator>
    static bool findPlane(InputIterator begin, InputIterator end, Plane3 & plane,
                          Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1);

    /**
     * Find a plane of (global) reflective symmetry of a set of objects. Passing a negative value for any of the numeric
     * parameters indicates that a default value should be chosen.
     *
     * @param begin The first objects.
     * @param end One position beyond the last object.
     * @param plane Used to store the resulting symmetry plane.
     * @param prox_query_struct A structure that quickly answers proximity queries on the input set.
     * @param precomputed_centroid The precomputed centroid of the objects. If NULL, the centroid will be computed from the
     *   input data.
     * @param max_rms_error The maximum allowed error, measured as the root of the average squared distance (RMS) of the
     *   reflections of the objects in the candidate plane to their respective closest neighbors.
     */
    template <typename InputIterator, typename ProximityQueryStructureT>
    static bool findPlane(InputIterator begin, InputIterator end, ProximityQueryStructureT const & prox_query_struct,
                          Plane3 & plane, Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1);

}; // class Symmetry3

// Symmetries of sets of objects passed as pointers.
template <typename T>
class /* THEA_API */ Symmetry3<T *>
{
  public:
    template <typename InputIterator>
    static bool findPlane(InputIterator begin, InputIterator end, Plane3 & plane,
                          Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1)
    {
      return Symmetry3<T>::findPlane(PtrToRefIterator<T, InputIterator>(begin), PtrToRefIterator<T, InputIterator>(end),
                                     plane, precomputed_centroid, max_rms_error);
    }

    template <typename InputIterator, typename ProximityQueryStructureT>
    static bool findPlane(InputIterator begin, InputIterator end, ProximityQueryStructureT const & prox_query_struct,
                          Plane3 & plane, Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1)
    {
      return Symmetry3<T>::findPlane(PtrToRefIterator<T, InputIterator>(begin), PtrToRefIterator<T, InputIterator>(end),
                                     prox_query_struct, plane, precomputed_centroid, max_rms_error);
    }

}; // class Symmetry3<T *>

// Symmetries of point clouds.
template <typename T>
class /* THEA_API */ Symmetry3<T, typename boost::enable_if< IsPointN<T, 3> >::type>
{
  public:
    template <typename InputIterator>
    static bool findPlane(InputIterator begin, InputIterator end, Plane3 & plane,
                          Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1)
    {
      KDTreeN<T, 3> kdtree(begin, end);
      return findPlane(begin, end, kdtree, plane, precomputed_centroid, max_rms_error);
    }

    template <typename InputIterator, typename ProximityQueryStructureT>
    static bool findPlane(InputIterator begin, InputIterator end, ProximityQueryStructureT const & prox_query_struct,
                          Plane3 & plane, Vector3 const * precomputed_centroid = NULL, double max_rms_error = -1)
    {
      if (begin == end)  // no points, early exit
        return false;

      // Choose a default RMS error if none was specified.
      double max_mean_squared_error;
      if (max_rms_error < 0)
      {
        // The threshold should be proportional to the extent of the point set
        AxisAlignedBox3 aab;
        long num_points = 0;
        for (InputIterator pi = begin; pi != end; ++pi, ++num_points)
          aab.merge(PointTraitsN<T, 3>::getPosition(*pi));

        double scale = num_points > 1000 ? 10.0 / num_points : 0.01;
        max_mean_squared_error = scale * scale * aab.getExtent().squaredLength();
      }
      else
        max_mean_squared_error = max_rms_error * max_rms_error;

      // Error attributed to an unmatched point
      double unmatched_error = 100 * max_mean_squared_error;  // allow 10 times the maximum mean separation

      // Compute PCA axes of the point set
      Real eigenvalues[3];
      Vector3 eigenvectors[3];
      Vector3 centroid;
      PCA_N<T, 3>::compute(begin, end, eigenvalues, eigenvectors, &centroid);

      double best_error = max_mean_squared_error;
      bool found = false;
      for (int i = 0; i < 3; ++i)
      {
        Vector3 test_normal = eigenvectors[i].fastUnit();
        Plane3 test_plane = Plane3::fromPointAndNormal(centroid, test_normal);

        double err = computeMeanSquaredError(test_plane, begin, end, prox_query_struct, unmatched_error);
        if (err < best_error)
        {
          plane = test_plane;
          best_error = err;
          found = true;
        }
      }

      return found;
    }

  private:
    // Compute the mean squared matching error between the reflection of the point set and the unreflected data.
    template <typename InputIterator, typename ProximityQueryStructureT>
    static double computeMeanSquaredError(Plane3 const & plane, InputIterator begin, InputIterator end,
                                          ProximityQueryStructureT const & prox_query_struct, double unmatched_error)
    {
      Vector3 normal = plane.getNormal();
      Vector3 cp;
      double sqrt_unmatched_error = std::sqrt(unmatched_error);
      double total_error = 0;
      long total_count = 0;
      for (InputIterator pi = begin; pi != end; ++pi, ++total_count)
      {
        Vector3 reflected = plane.reflect(PointTraitsN<T, 3>::getPosition(*pi));

        // Without the ".template", C++ thinks the '<' is "less-than". See Vandevoorde and Josuttis, "C++ Templates: The
        // Complete Guide", p.44.
        long index = prox_query_struct.template closestElement<MetricL2>(reflected, sqrt_unmatched_error, NULL, &cp);
        if (index >= 0)
          total_error += (reflected - cp).squaredLength();
        else
          total_error += unmatched_error;
      }

      return total_count <= 0 ? 0 : total_error / total_count;
    }

}; // class Symmetry3<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
