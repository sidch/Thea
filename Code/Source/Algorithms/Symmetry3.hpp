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
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE),
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#ifndef __Thea_Algorithms_Symmetry3_hpp__
#define __Thea_Algorithms_Symmetry3_hpp__

#include "../Common.hpp"
#include "CentroidN.hpp"
#include "GeodesicSphere3.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"
#include "../Math.hpp"
#include "../UnorderedSet.hpp"
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
     * Find a plane of (global) reflective symmetry of a set of objects.
     *
     * @param begin The first object.
     * @param end One position beyond the last object.
     * @param plane Used to store the resulting symmetry plane.
     * @param precomputed_centroid The precomputed centroid of the objects. If NULL, the centroid will be computed from the
     *   input data.
     *
     * @return The error, in the range 0 (best) to 1 (worst), of the symmetry relation.
     */
    template <typename InputIterator>
    static double findPlane(InputIterator begin, InputIterator end, Plane3 & plane,
                            Vector3 const * precomputed_centroid = NULL);

}; // class Symmetry3

// Symmetries of sets of objects passed as pointers.
template <typename T>
class /* THEA_API */ Symmetry3<T *>
{
  public:
    template <typename InputIterator>
    static double findPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 const * precomputed_centroid = NULL)
    {
      return Symmetry3<T>::findPlane(PtrToRefIterator<T, InputIterator>(begin), PtrToRefIterator<T, InputIterator>(end),
                                     plane, precomputed_centroid);
    }

}; // class Symmetry3<T *>

// Symmetries of point clouds.
template <typename T>
class /* THEA_API */ Symmetry3<T, typename boost::enable_if< IsNonReferencedPointN<T, 3> >::type>
{
  public:
    template <typename InputIterator>
    static double findPlane(InputIterator begin, InputIterator end, Plane3 & plane, Vector3 const * precomputed_centroid = NULL,
                            long num_rounds = -1)
    {
      if (begin == end)  // no points, early exit
        return false;

      if (num_rounds <= 0)
        num_rounds = 5;

      Vector3 centroid = (precomputed_centroid ? *precomputed_centroid : CentroidN<T, 3>::compute(begin, end));
      Real radius = 0;
      for (InputIterator pi = begin; pi != end; ++pi)
        radius = std::max(radius, (PointTraitsN<T, 3>::getPosition(*pi) - centroid).squaredLength());

      radius = std::sqrt(radius);

      TheaArray<Vector3> vertices;
      TheaArray<long> triangles;
      GeodesicSphere3::compute(2, vertices, &triangles);

      long best_dir = -1;
      double best_error = -1;
      size_t first_vertex_of_round = 0;
      for (long round = 0; round < num_rounds; ++round)
      {
        for (size_t i = first_vertex_of_round; i < vertices.size(); ++i)
        {
          Plane3 candidate = Plane3::fromPointAndNormal(centroid, vertices[i]);
          double err = symmetryError(begin, end, candidate, centroid, radius);
          if (best_dir < 0 || err < best_error)
          {
            best_dir = (long)i;
            best_error = err;
            plane = candidate;
          }
        }

        THEA_DEBUG << "GeodesicSphere3: Round " << round << ", best plane = " << plane.toString() << ", error = " << best_error;

        // Further localize the search to a smaller set of directions more tightly clustered around the best one
        if (round < num_rounds - 1)
        {
          first_vertex_of_round = vertices.size();
          TheaArray<long> tris_to_subdivide;

          // Collect triangles incident on the best direction vertex
          TheaUnorderedSet<long> nbr_verts;
          for (size_t i = 0; i < triangles.size(); i += 3)
          {
            if ((long)triangles[i] == best_dir || (long)triangles[i + 1] == best_dir || (long)triangles[i + 2] == best_dir)
            {
              tris_to_subdivide.push_back(triangles[i    ]);
              tris_to_subdivide.push_back(triangles[i + 1]);
              tris_to_subdivide.push_back(triangles[i + 2]);

              if (triangles[i    ] != best_dir) nbr_verts.insert(triangles[i    ]);
              if (triangles[i + 1] != best_dir) nbr_verts.insert(triangles[i + 1]);
              if (triangles[i + 2] != best_dir) nbr_verts.insert(triangles[i + 2]);
            }
          }

          // Collect triangles incident on the one-hop neighbor vertices
          for (size_t i = 0; i < triangles.size(); i += 3)
          {
            if ((long)triangles[i] == best_dir || (long)triangles[i + 1] == best_dir || (long)triangles[i + 2] == best_dir)
              continue;  // already added

            if (nbr_verts.find((long)triangles[i    ]) != nbr_verts.end()
             || nbr_verts.find((long)triangles[i + 1]) != nbr_verts.end()
             || nbr_verts.find((long)triangles[i + 2]) != nbr_verts.end())
            {
              tris_to_subdivide.push_back(triangles[i    ]);
              tris_to_subdivide.push_back(triangles[i + 1]);
              tris_to_subdivide.push_back(triangles[i + 2]);
            }
          }

          triangles.clear();
          GeodesicSphere3::compute(1, vertices, tris_to_subdivide, &triangles);
        }
      }

      return best_error;
    }

  private:
    /** Measure the quality of a candidate symmetry plane. Returns a number between 0 (best) and 1 (worst). */
    template <typename InputIterator>
    static double symmetryError(InputIterator begin, InputIterator end, Plane3 const & plane, Vector3 const & centroid,
                                Real radius)
    {
      static size_t const NUM_BINS = 10;

      double bins[2][NUM_BINS][NUM_BINS];  // 0: negative side, 1: positive side
      long count[2][NUM_BINS][NUM_BINS];

      for (size_t i = 0; i < NUM_BINS; ++i)
        for (size_t j = 0; j < NUM_BINS; ++j)
        {
          bins[0][i][j] = bins[1][i][j] = 0.0;
          count[0][i][j] = count[1][i][j] = 0;
        }

      Vector3 nrm = plane.getNormal();
      Vector3 u, v; nrm.createOrthonormalBasis(u, v);

      long total_count[2] = { 0, 0 };
      for (InputIterator pi = begin; pi != end; ++pi)
      {
        Vector3 p = PointTraitsN<T, 3>::getPosition(*pi) - centroid;
        Real pu = 0.5f * (p.dot(u) / radius + 1);
        Real pv = 0.5f * (p.dot(v) / radius + 1);

        int bin_i = Math::clamp((int)(NUM_BINS * pu), 0, NUM_BINS - 1);
        int bin_j = Math::clamp((int)(NUM_BINS * pv), 0, NUM_BINS - 1);

        Real pw = p.dot(nrm);
        if (pw >= 0)
        {
          bins[1][bin_i][bin_j] += pw;
          count[1][bin_i][bin_j]++;
          total_count[1]++;
        }

        if (pw <= 0)  // double-count points exactly on the plane, if any
        {
          bins[0][bin_i][bin_j] += std::fabs(pw);
          count[0][bin_i][bin_j]++;
          total_count[0]++;
        }
      }

      double dist_error = 0, count_error = 0;
      for (size_t i = 0; i < NUM_BINS; ++i)
        for (size_t j = 0; j < NUM_BINS; ++j)
        {
          long c0 = count[0][i][j];
          long c1 = count[1][i][j];
          if (c0 == 0 && c1 == 0)
            continue;

          double avg0 = c0 > 0 ? bins[0][i][j] / c0 : 0.0;
          double avg1 = c1 > 0 ? bins[1][i][j] / c1 : 0.0;

          double weight = c0 + c1;
          dist_error += (std::fabs(avg0 - avg1) / radius) * weight;
          count_error += (c0 < c1 ? 1.0 - c0 / (double)c1 : 1.0 - c1 / (double)c0) * weight;
        }

      double sum_weights = total_count[0] + total_count[1];
      dist_error /= sum_weights;
      count_error /= sum_weights;

      return 0.5f * (dist_error + count_error);
    }

}; // class Symmetry3<Point3>

} // namespace Algorithms
} // namespace Thea

#endif
