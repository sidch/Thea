//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2018, Siddhartha Chaudhuri
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

#include "../Common.hpp"
#include "FurthestPointSampling.hpp"
#include "SampleGraph.hpp"
#include "ShortestPaths.hpp"
#include "../Vector3.hpp"
#include "../UnorderedMap.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>

namespace Thea {
namespace Algorithms {

namespace FurthestPointSamplingInternal {

// Called during Dijkstra search.
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

} // namespace FurthestPointSamplingInternal

long
FurthestPointSampling::subsample(long num_orig_points, Vector3 const * orig_points, long num_desired_points,
                                 long * selected_indices, DistanceType dist_type, bool verbose)
{
  alwaysAssertM(num_desired_points >= 0, "FurthestPointSampling: Can't sample a negative number of points");
  alwaysAssertM(num_orig_points >= num_desired_points,
                format("FurthestPointSampling: Can't subsample %ld point(s) from %ld point(s)",
                       num_desired_points, num_orig_points));
  alwaysAssertM(dist_type == DistanceType::GEODESIC, "FurthestPointSampling: Only geodesic distances currently supported");

  if (num_desired_points == 0)
    return 0;

  // Compute proximity graph
  SampleGraph graph;
  graph.setSamples(num_orig_points, orig_points);
  graph.init();

  if (verbose)
  {
    THEA_CONSOLE << "FurthestPointSampling: Computed proximity graph";
    std::cout << "FurthestPointSampling: Selecting samples: " << std::flush;
  }

  // Repeatedly add the furthest sample from the selected set to the selected set
  ShortestPaths<SampleGraph> shortest_paths;
  TheaUnorderedMap<SampleGraph::VertexHandle, double> src_region;
  int prev_percent = 0;
  for (long i = 0; i < num_desired_points; ++i)
  {
    SampleGraph::VertexHandle furthest_sample = NULL;
    if (src_region.empty())
    {
      // Just pick the first sample
      furthest_sample = const_cast<SampleGraph::VertexHandle>(&graph.getSample(0));
    }
    else
    {
      FurthestPointSamplingInternal::DijkstraCallback callback;
      shortest_paths.dijkstraWithCallback(graph, NULL, &callback, -1, &src_region, /* include_unreachable = */ true);
      if (!callback.furthest_sample || src_region.find(callback.furthest_sample) != src_region.end())
      {
        THEA_ERROR << "FurthestPointSampling: Could not return enough uniformly separated points";
        return i;
      }

      furthest_sample = callback.furthest_sample;
    }

    alwaysAssertM(furthest_sample, "FurthestPointSampling: Furthest point is null");

    selected_indices[i] = furthest_sample->getIndex();
    src_region[furthest_sample] = 0;

    if (verbose)
    {
      int curr_percent = (int)std::floor(100 * (i / (float)num_desired_points));
      if (curr_percent >= prev_percent + 2)
      {
        for (prev_percent += 2; prev_percent <= curr_percent; prev_percent += 2)
        {
          if (prev_percent % 10 == 0) std::cout << prev_percent << '%' << std::flush;
          else                        std::cout << '.' << std::flush;
        }

        prev_percent = curr_percent;
      }
    }
  }

  if (verbose)
  {
    std::cout << "done" << std::endl;
    THEA_CONSOLE << "FurthestPointSampling: Selected " << num_desired_points << " point(s) uniformly by separation";
  }

  return num_desired_points;
}

} // namespace Algorithms
} // namespace Thea
