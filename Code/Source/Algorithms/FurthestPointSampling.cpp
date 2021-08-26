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
// First version: 2018
//
//============================================================================

#include "../Common.hpp"
#include "FurthestPointSampling.hpp"
#include "PointCloud3.hpp"
#include "ShortestPaths.hpp"
#include "../Noncopyable.hpp"
#include "../UnorderedMap.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>

namespace Thea {
namespace Algorithms {

namespace FurthestPointSamplingInternal {

// Called during Dijkstra search.
struct DijkstraCallback : public Noncopyable
{
  DijkstraCallback() : furthest_sample(nullptr) {}

  bool operator()(PointCloud3::SampleGraph::VertexHandle vertex, double distance, bool has_pred,
                  PointCloud3::SampleGraph::VertexHandle pred)
  {
    furthest_sample = vertex;  // the callback is always called in order of increasing distance
    return false;
  }

  PointCloud3::SampleGraph::VertexHandle furthest_sample;
};

} // namespace FurthestPointSamplingInternal

intx
FurthestPointSampling::subsample(intx num_orig_points, Vector3 const * orig_points, intx num_desired_points,
                                 intx * selected_indices, DistanceType dist_type, bool verbose)
{
  alwaysAssertM(num_desired_points >= 0, "FurthestPointSampling: Can't sample a negative number of points");
  alwaysAssertM(num_orig_points >= num_desired_points,
                format("FurthestPointSampling: Can't subsample %ld point(s) from %ld point(s)",
                       num_desired_points, num_orig_points));
  alwaysAssertM(dist_type == DistanceType::GEODESIC, "FurthestPointSampling: Only geodesic distances currently supported");

  if (num_desired_points == 0)
    return 0;

  // Compute proximity graph
  PointCloud3 surf;
  surf.addSamples(num_orig_points, orig_points);

  if (verbose)
  {
    THEA_CONSOLE << "FurthestPointSampling: Constructed point-sampled surface";
    std::cout << "FurthestPointSampling: Selecting samples: " << std::flush;
  }

  // Repeatedly add the furthest sample from the selected set to the selected set
  ShortestPaths<PointCloud3::SampleGraph> shortest_paths;
  UnorderedMap<PointCloud3::SampleGraph::VertexHandle, double> src_region;
  int prev_percent = 0;
  for (intx i = 0; i < num_desired_points; ++i)
  {
    PointCloud3::SampleGraph::VertexHandle furthest_sample = nullptr;
    if (src_region.empty())
    {
      // Just pick the first sample
      furthest_sample = const_cast<PointCloud3::SampleGraph::VertexHandle>(&surf.getSample(0));
    }
    else
    {
      FurthestPointSamplingInternal::DijkstraCallback callback;
      shortest_paths.dijkstraWithCallback(const_cast<PointCloud3::SampleGraph &>(surf.getGraph()), nullptr, std::ref(callback),
                                          -1, &src_region, /* include_unreachable = */ true);
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
