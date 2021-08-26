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
// First version: 2015
//
//============================================================================

#include "DistanceHistogram.hpp"
#include "../../../Math.hpp"
#include "../../../Random.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Global {

DistanceHistogram::DistanceHistogram(PointCloud3 const * surf_)
: ldh(surf_)
{}

void
DistanceHistogram::compute(Histogram & histogram, DistanceType dist_type, Real max_distance, Real pair_reduction_ratio) const
{
  intx num_samples = getSurface()->numSamples();
  intx num_distinct_unordered = (num_samples - 1) * num_samples;

  if (pair_reduction_ratio < 0)
    pair_reduction_ratio = (Real)std::min(0.5, 1000000.0 / num_distinct_unordered);  // don't count (x, x)

  histogram.setZero();  // don't bother setting the range

  Real local_reduction_ratio = std::sqrt(pair_reduction_ratio);
  intx num_queries = Math::clamp((intx)std::ceil(local_reduction_ratio * num_samples), 0, num_samples - 1);
  if (num_queries <= 0)
    return;

  Array<int32> query_indices((size_t)num_queries);
  Random::common().sortedIntegers(0, (int32)num_samples - 1, (int32)num_queries, &query_indices[0]);

  Histogram local_histogram(histogram.numBins());
  auto const & samples = getSurface()->getSamples();
  for (size_t i = 0; i < query_indices.size(); ++i)
  {
    Vector3 const & p = samples[(size_t)query_indices[i]].getPosition();
    ldh.compute(p, local_histogram, dist_type, max_distance, local_reduction_ratio);

    // Remove the zero distance from the query point to itself
    local_histogram.remove(0.0);

    if (i == 0)
      histogram.setRange(local_histogram.minValue(), local_histogram.maxValue());

    histogram.insert(local_histogram);
  }
}

} // namespace Global
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
