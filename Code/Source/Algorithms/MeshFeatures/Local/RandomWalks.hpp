//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2016, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_MeshFeatures_Local_RandomWalks_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_RandomWalks_hpp__

#include "../../../Common.hpp"
#include "../SampledSurface.hpp"
#include "../../Histogram.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../PointTraitsN.hpp"
#include "../../SampleGraph.hpp"
#include "../../ShortestPaths.hpp"
#include "../../../Ball3.hpp"
#include "../../../Random.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/** Compute the average offset, from the query position, after each step of an n-step random walk. */
template < typename ExternalSampleKDTreeT = KDTreeN<Vector3, 3> >
class RandomWalks : public SampledSurface<ExternalSampleKDTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKDTreeT> BaseT;  ///< Base class.
    static long const DEFAULT_NUM_SAMPLES = 5000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute random walk patterns on a given mesh. The mesh need not persist beyond the execution of
     * the constructor.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     */
    template <typename MeshT>
    RandomWalks(MeshT const & mesh, long num_samples = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples))
    {}

    /**
     * Constructs the object to compute random walk patterns on a given mesh group. The mesh group need not persist beyond the
     * execution of the constructor.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     */
    template <typename MeshT>
    RandomWalks(Graphics::MeshGroup<MeshT> const & mesh_group, long num_samples = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples))
    {}

    /**
     * Constructs the object to compute random walk patterns on a shape with a precomputed kd-tree on surface samples. The
     * kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ The kd-tree representing the shape.
     */
    RandomWalks(ExternalSampleKDTreeT const * sample_kdtree_)
    : BaseT(sample_kdtree_)
    {}

    /**
     * Compute the the average offset, from the query position, after each step of an n-step random walk on the shape's sample
     * graph.
     *
     * @param position The position of the query point.
     * @param num_steps The number of random steps to take from the query point.
     * @param features Used to return the point features. Should be preallocated to \a num_steps * 3 entries.
     * @param num_walks The number of random walks over which to take averages.
     */
    void compute(Vector3 const & position, long num_steps, double * features, long num_walks = -1) const
    {
      static long const DEFAULT_NUM_WALKS = 1000;

      if (num_walks < 0)
        num_walks = DEFAULT_NUM_WALKS;

      // Zero out all features initially
      for (long i = 0; i < num_steps; ++i)
      {
        features[3 * i    ] = 0.0;
        features[3 * i + 1] = 0.0;
        features[3 * i + 2] = 0.0;
      }

      // Find the sample closest to the query position and use it as the source for all distance calculations
      long seed_index = -1;
      if (this->hasExternalKDTree())
        seed_index = this->getMutableExternalKDTree()->template closestElement<MetricL2>(position);
      else
        seed_index = this->getMutableInternalKDTree()->template closestElement<MetricL2>(position);

      alwaysAssertM(seed_index >= 0, "RandomWalks: Seed sample for random walks not found");

      TheaArray<long> counts((size_t)num_steps, 0);
      for (long i = 0; i < num_walks; ++i)
      {
        long walk_steps = walk(seed_index, num_steps, features);
        for (long j = 0; j < walk_steps; ++j)
          counts[(size_t)j]++;
      }

      for (long i = 0; i < num_steps; ++i)
      {
        long count = counts[(size_t)i];
        features[3 * i    ] = (features[3 * i    ] / count - position[0]);
        features[3 * i + 1] = (features[3 * i + 1] / count - position[1]);
        features[3 * i + 2] = (features[3 * i + 2] / count - position[2]);
      }
    }

  private:
    /**
     * Do a random walk upto \a num_steps steps, and return the number of steps actually taken (= \a num_steps except in corner
     * cases).
     */
    long walk(long seed_index, long num_steps, double * features) const
    {
      typedef SampleGraph::SurfaceSample::NeighborSet NeighborSet;

      // Assume the graph and the kd-tree have samples in the same sequence
      SampleGraph * graph = const_cast<SampleGraph *>(this->getSampleGraph());
      alwaysAssertM(graph, "RandomWalks: Non-null sample graph required for random walks");

      SampleGraph::SurfaceSample * sample = const_cast<SampleGraph::SurfaceSample *>(&graph->getSample(seed_index));

      long base_index = 0;
      for (long i = 0; i < num_steps; ++i, base_index += 3)
      {
        NeighborSet const & nbrs = sample->getNeighbors();
        if (nbrs.isEmpty())
          return i;

        long next_index = Random::common().integer(0, nbrs.size() - 1);
        sample = nbrs[next_index].getSample();

        Vector3 const & p = sample->getPosition();
        features[base_index    ] += p[0];
        features[base_index + 1] += p[1];
        features[base_index + 2] += p[2];
      }

      return num_steps;
    }

}; // class RandomWalks

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
