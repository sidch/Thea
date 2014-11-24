//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_MeshSampler_hpp__
#define __Thea_Algorithms_MeshSampler_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "MeshTriangles.hpp"
#include "SampleGraph.hpp"
#include "ShortestPaths.hpp"
#include "../Random.hpp"
#include "../Vector3.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace Thea {
namespace Algorithms {

/** Sample points from a mesh. */
template <typename MeshT>
class MeshSampler
{
  public:
    typedef MeshT                         Mesh;       ///< The mesh class.
    typedef MeshTriangles<Mesh>           Triangles;  ///< A set of mesh triangles.
    typedef typename Triangles::Triangle  Triangle;   ///< A triangle belonging to a face of the mesh.

    /** Whether to return exactly or approximately the desired number of sampled points. */
    struct CountMode
    {
      /** Supported values. */
      enum Value
      {
        EXACT,       ///< Return exactly the desired number of samples.
        APPROXIMATE  ///< Return approximately the desired number of samples (likely to be faster than EXACT).
      };

      THEA_ENUM_CLASS_BODY(CountMode)
    };

    /**
     * Constructs the object to compute sample points on a given mesh. The mesh must persist as long as this object does.
     * Initializes internal data structures that do not need to be recomputed for successive calls to functions that generate
     * samples on this mesh.
     */
    MeshSampler(Mesh const & mesh) : num_external_tris(0), external_tris(NULL)
    {
      tris.add(const_cast<Mesh &>(mesh));
    }

    /**
     * Constructs the object to compute sample points on a given mesh group. The mesh group must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to functions that
     * generate samples on this mesh.
     */
    MeshSampler(Graphics::MeshGroup<Mesh> const & mesh_group) : num_external_tris(0), external_tris(NULL)
    {
      tris.add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));
    }

    /**
     * Constructs the object to compute sample points on a mesh with a precomputed triangulation. The triangles must persist as
     * long as this object does.
     */
    MeshSampler(long num_tris, Triangle const * tris) : num_external_tris(num_tris), external_tris(tris)
    {
      alwaysAssertM(num_tris <= 0 || tris, "MeshSampler: Triangle list cannot be null");
    }

    /** Destructor. */
    ~MeshSampler() {}

    /**
     * Samples points from the mesh evenly by area.
     *
     * @param desired_num_samples The number of samples desired. This need <em>not</em> be exactly the number actually computed
     *   (and returned by the function), unless \a count_mode == CountMode::EXACT.
     * @param positions Used to return the sample positions.
     * @param face_normals If not null, used to return the normals of the parent faces of the samples.
     * @param triangles If not null, used to return the mesh triangles from which the samples were selected.
     * @param count_mode If CountMode::EXACT, exactly \a desired_num_samples samples are returned. Else, approximately
     *   \a desired_num_samples samples are returned. The latter is likely to be faster.
     * @param verbose Prints progress messages.
     *
     * @return The number of samples computed.
     */
    long sampleEvenlyByArea(long desired_num_samples,
                            TheaArray<Vector3> & positions,
                            TheaArray<Vector3> * face_normals = NULL,
                            TheaArray<Triangle const *> * triangles = NULL,
                            CountMode count_mode = CountMode::EXACT,
                            bool verbose = false) const
    {
      positions.clear();
      if (face_normals) face_normals->clear();
      if (triangles) triangles->clear();

      if (desired_num_samples <= 0)
      {
        THEA_CONSOLE << "MeshSampler: Selected 0 sample(s) uniformly by area";
        return 0;
      }

      Triangle const * tarray;
      long tcount;
      if (num_external_tris > 0)
      {
        tarray = external_tris;
        tcount = num_external_tris;
      }
      else
      {
        if (tris.isEmpty())
        {
          THEA_CONSOLE << "MeshSampler: Selected 0 sample(s) uniformly by area";
          return 0;
        }
        else
        {
          tarray = &tris.getTriangles()[0];
          tcount = tris.numTriangles();
        }
      }

      if (verbose)
        THEA_CONSOLE << "MeshSampler: Obtained triangle list";

      if (tcount <= 0)
      {
        THEA_CONSOLE << "MeshSampler: Selected 0 sample(s) uniformly by area";
        return 0;
      }

      double total_area = 0;
      for (long i = 0; i < tcount; ++i)
        total_area += tarray[i].getArea();

      if (total_area <= 0)
      {
        THEA_CONSOLE << "MeshSampler: Selected 0 sample(s) uniformly by area";
        return 0;
      }

      if (verbose)
      {
        THEA_CONSOLE << "MeshSampler: Computed triangle areas";
        std::cout << "MeshSampler: Selecting samples: " << std::flush;
      }

      int prev_percent = 0;
      if (count_mode == CountMode::EXACT)
      {
        TheaArray<double> cum_areas((array_size_t)tcount);
        cum_areas[0] = tarray[0].getArea();
        for (array_size_t i = 1; i < cum_areas.size(); ++i)
          cum_areas[i] = cum_areas[i - 1] + tarray[i].getArea();

        positions.resize((array_size_t)desired_num_samples);
        if (face_normals) face_normals->resize(positions.size());
        if (triangles) triangles->resize(positions.size());

        // Sample N triangles with replacement, weighted by area
        double const * first = &cum_areas[0];
        double const * last = first + cum_areas.size();
        for (array_size_t i = 0; i < positions.size(); ++i)
        {
          double t = Random::common().uniform01() * total_area;
          long tri_index = std::min((long)(std::lower_bound(first, last, t) - first), tcount - 1);

          positions[i] = tarray[tri_index].randomPoint();
          if (face_normals) (*face_normals)[i] = tarray[tri_index].getNormal();
          if (triangles) (*triangles)[i] = &tarray[tri_index];

          if (verbose)
          {
            int curr_percent = (int)std::floor(100 * (i / (float)desired_num_samples));
            if (curr_percent >= prev_percent + 2)
            {
              for (prev_percent += 2; prev_percent <= curr_percent; prev_percent += 2)
              {
                if (prev_percent % 10 == 0)
                  std::cout << prev_percent << '%' << std::flush;
                else
                  std::cout << '.' << std::flush;
              }

              prev_percent = curr_percent;
            }
          }
        }
      }
      else
      {
        int prev_percent = 0;
        for (long i = 0; i < tcount; ++i)
        {
          double frac_samples = (desired_num_samples * tarray[i].getArea()) / total_area;
          while (frac_samples > 0)
          {
            // Last (possible) sample?
            if (frac_samples < 1 && Random::common().uniform01() > frac_samples)
              break;

            positions.push_back(tarray[i].randomPoint());
            if (face_normals) face_normals->push_back(tarray[i].getNormal());
            if (triangles) triangles->push_back(&tarray[i]);

            frac_samples -= 1.0;

            if (verbose)
            {
              int curr_percent = (int)std::floor(100 * (positions.size() / (float)desired_num_samples));
              if (curr_percent >= prev_percent + 2 && curr_percent < 100)
              {
                for (prev_percent += 2; prev_percent <= curr_percent; prev_percent += 2)
                {
                  if (prev_percent % 10 == 0)
                    std::cout << prev_percent << '%' << std::flush;
                  else
                    std::cout << '.' << std::flush;
                }

                prev_percent = curr_percent;
              }
            }
          }
        }
      }

      if (verbose)
      {
        std::cout << "done" << std::endl;
        THEA_CONSOLE << "MeshSampler: Selected " << positions.size() << " sample(s) uniformly by area";
      }

      return (long)positions.size();
    }

    /**
     * Samples approximately uniformly separately points from the mesh. The function returns an array of samples ordered so that
     * for any K, the first K samples form an approximately uniformly separated subsampling. The function initially computes an
     * oversampling of the mesh by area.
     *
     * @param desired_num_samples The number of samples desired. This need <em>not</em> be exactly the number actually computed
     *   (and returned by the function), unless \a count_mode == CountMode::EXACT.
     * @param positions Used to return the sample positions.
     * @param face_normals If not null, used to return the normals of the parent faces of the samples.
     * @param triangles If not null, used to return the mesh triangles from which the samples were selected.
     * @param count_mode If CountMode::EXACT, exactly \a desired_num_samples samples are returned. Else, approximately
     *   \a desired_num_samples samples are returned. The latter is likely to be faster.
     * @param oversampling_factor The ratio of the number of initially selected points to the number of finally returned points.
     *   The default, selected with a negative value, is 3.
     * @param verbose Prints progress messages.
     *
     * @return The number of samples computed.
     */
    long sampleEvenlyBySeparation(long desired_num_samples,
                                  TheaArray<Vector3> & positions,
                                  TheaArray<Vector3> * face_normals = NULL,
                                  TheaArray<Triangle const *> * triangles = NULL,
                                  CountMode count_mode = CountMode::EXACT,
                                  Real oversampling_factor = -1,
                                  bool verbose = false) const
    {
      positions.clear();
      if (face_normals) face_normals->clear();
      if (triangles) triangles->clear();

      if (desired_num_samples <= 0)
      {
        THEA_CONSOLE << "MeshSampler: Selected 0 sample(s) uniformly by separation";
        return 0;
      }

      // Compute oversampling by area
      TheaArray<Vector3> orig_positions;
      TheaArray<Vector3> orig_face_normals;
      TheaArray<Triangle const *> orig_triangles;

      static Real const DEFAULT_OVERSAMPLING_FACTOR = 3;
      if (oversampling_factor < 0)
        oversampling_factor = DEFAULT_OVERSAMPLING_FACTOR;

      long num_oversampling = (long)std::ceil(oversampling_factor * desired_num_samples);
      long orig_num_samples = sampleEvenlyByArea(num_oversampling, orig_positions,
                                                 face_normals ? &orig_face_normals : NULL,
                                                 triangles ? &orig_triangles : NULL, CountMode::EXACT);
      if (orig_num_samples < num_oversampling)
      {
        THEA_ERROR << "MeshSampler: Could not compute oversampling";
        return 0;
      }

      if (verbose)
      {
        THEA_CONSOLE << "MeshSampler: Computed " << oversampling_factor << "x oversampling with " << orig_num_samples
                     << " sample(s)";
      }

      // Compute proximity graph
      SampleGraph graph;
      graph.setSamples(orig_num_samples, &orig_positions[0]);
      graph.init();

      if (verbose)
      {
        THEA_CONSOLE << "MeshSampler: Computed proximity graph";
        std::cout << "MeshSampler: Selecting samples: " << std::flush;
      }

      // Repeatedly add the furthest sample from the selected set to the selected set
      ShortestPaths<SampleGraph> shortest_paths;
      TheaUnorderedMap<SampleGraph::VertexHandle, double> src_region;
      int prev_percent = 0;
      for (long i = 0; i < desired_num_samples; ++i)
      {
        SampleGraph::VertexHandle furthest_sample = NULL;
        if (src_region.empty())
        {
          // Just pick the first sample
          furthest_sample = const_cast<SampleGraph::VertexHandle>(&graph.getSample(0));
        }
        else
        {
          DijkstraCallback callback;
          shortest_paths.dijkstraWithCallback(graph, NULL, &callback, -1, &src_region, /* include_unreachable = */ true);
          if (!callback.furthest_sample || src_region.find(callback.furthest_sample) != src_region.end())
          {
            THEA_WARNING << "MeshSampler: Could not return enough uniformly separated samples";
            break;
          }

          furthest_sample = callback.furthest_sample;
        }

        if (furthest_sample)
        {
          array_size_t index = (array_size_t)furthest_sample->getIndex();

          positions.push_back(orig_positions[index]);
          if (face_normals) face_normals->push_back(orig_face_normals[index]);
          if (triangles) triangles->push_back(orig_triangles[index]);

          src_region[furthest_sample] = 0;
        }
        else
          break;

        if (verbose)
        {
          int curr_percent = (int)std::floor(100 * (i / (float)desired_num_samples));
          if (curr_percent >= prev_percent + 2)
          {
            for (prev_percent += 2; prev_percent <= curr_percent; prev_percent += 2)
            {
              if (prev_percent % 10 == 0)
                std::cout << prev_percent << '%' << std::flush;
              else
                std::cout << '.' << std::flush;
            }

            prev_percent = curr_percent;
          }
        }
      }

      if (verbose)
      {
        std::cout << "done" << std::endl;
        THEA_CONSOLE << "MeshSampler: Selected " << positions.size() << " sample(s) uniformly by separation";
      }

      return (long)positions.size();
    }

  private:
    /** Called during Dijkstra search. */
    struct DijkstraCallback
    {
      DijkstraCallback() : furthest_sample(NULL) {}

      bool operator()(SampleGraph::VertexHandle vertex, double distance, bool has_pred, SampleGraph::VertexHandle pred)
      {
        furthest_sample = vertex;  // the callback is always called in order of increasing distance
        return false;
      }

      SampleGraph::VertexHandle furthest_sample;
    };

    Triangles tris;  ///< An internally computed triangulation of the mesh.
    long num_external_tris;  ///< Number of externally supplied, precomputed mesh triangles.
    Triangle const * external_tris;  ///< Array of externally supplied, precomputed mesh triangles.

}; // class MeshSampler

} // namespace Algorithms
} // namespace Thea

#endif
