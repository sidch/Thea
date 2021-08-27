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
// First version: 2016
//
//============================================================================

#include "RandomWalks.hpp"
#include "../../MetricL2.hpp"
#include "../../../Random.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

RandomWalks::RandomWalks(PointSet3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "RandomWalks: Cannot construct with a null surface");
}

namespace RandomWalksInternal {

// Do a random walk upto \a num_steps steps, and return the number of steps actually taken (= \a num_steps except in corner
// cases).
intx
walk(PointSet3 const & surf, intx seed_index, intx num_steps, double * features)
{
  SamplePoint3 const * sample = &surf.getSample(seed_index);
  intx base_index = 0;
  for (intx i = 0; i < num_steps; ++i, base_index += 3)
  {
    auto const & nbrs = sample->getNeighbors();
    if (nbrs.isEmpty())
      return i;

    intx next_index = Random::common().integer(0, nbrs.size() - 1);
    sample = nbrs[next_index].getSample();

    Vector3 const & p = sample->getPosition();
    features[base_index    ] += p[0];
    features[base_index + 1] += p[1];
    features[base_index + 2] += p[2];
  }

  return num_steps;
}

} // namespace RandomWalksInternal

/**
 * Compute the average offset, from the query position, after each step of an n-step random walk on the shape's sample
 * graph.
 *
 * @param position The position of the query point.
 * @param num_steps The number of random steps to take from the query point.
 * @param features Used to return the point features. Should be preallocated to \a num_steps * 3 entries.
 * @param num_walks The number of random walks over which to take averages.
 */
void
RandomWalks::compute(Vector3 const & position, intx num_steps, double * features, intx num_walks) const
{
  static intx const DEFAULT_NUM_WALKS = 1000;

  if (num_walks < 0)
    num_walks = DEFAULT_NUM_WALKS;

  // Zero out all features initially
  for (intx i = 0; i < num_steps; ++i)
  {
    features[3 * i    ] = 0.0;
    features[3 * i + 1] = 0.0;
    features[3 * i + 2] = 0.0;
  }

  // Find the sample closest to the query position and use it as the source for walks
  intx seed_index = surf->getKdTree().closestElement<MetricL2>(position);
  alwaysAssertM(seed_index >= 0, "RandomWalks: Seed sample for random walks not found");

  // Make sure the proximity graph exists, since we don't actually call getGraph()
  surf->updateGraph();

  // Do the walks
  Array<intx> counts((size_t)num_steps, 0);
  for (intx i = 0; i < num_walks; ++i)
  {
    intx walk_steps = RandomWalksInternal::walk(*surf, seed_index, num_steps, features);
    for (intx j = 0; j < walk_steps; ++j)
      counts[(size_t)j]++;
  }

  // Compute average offsets
  for (intx i = 0; i < num_steps; ++i)
  {
    intx count = counts[(size_t)i];
    features[3 * i    ] = (features[3 * i    ] / count - position[0]);
    features[3 * i + 1] = (features[3 * i + 1] / count - position[1]);
    features[3 * i + 2] = (features[3 * i + 2] / count - position[2]);
  }
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
