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

#include "LocalDistanceHistogram.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../PointTraitsN.hpp"
#include "../../ShortestPaths.hpp"
#include "../../../Ball3.hpp"
#include "../../../Random.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

LocalDistanceHistogram::LocalDistanceHistogram(PointCloud3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "LocalDistanceHistogram: Cannot construct with a null surface");
}

namespace LocalDistanceHistogramInternal {

// Called for each point in the euclidean neighborhood.
struct EuclideanCallback
{
  EuclideanCallback(Vector3 const & position_, Histogram & histogram_, Real acceptance_probability_)
  : position(position_), histogram(histogram_), acceptance_probability(acceptance_probability_)
  {}

  template <typename SampleT> bool operator()(intx index, SampleT const & t)
  {
    if (acceptance_probability < 1 && Random::common().uniform01() > acceptance_probability)
      return false;

    Real d = (PointTraitsN<SampleT, 3>::getPosition(t) - position).norm();
    histogram.insert(d);

    return false;
  }

  Vector3 position;
  Histogram & histogram;
  Real acceptance_probability;

}; // struct EuclideanCallback

void
computeEuclidean(PointCloud3 const & surf, Vector3 const & position, Histogram & histogram, Real max_distance,
                 Real sample_reduction_ratio)
{
  if (sample_reduction_ratio < 0)
    sample_reduction_ratio = 1.1;  // play safe

  bool process_all = (max_distance < 0);
  if (process_all)
    max_distance = surf.getScale();

  histogram.setRange(0, std::max((double)max_distance, 1.0e-30));
  histogram.setZero();

  EuclideanCallback callback(position, histogram, sample_reduction_ratio);

  if (process_all)
  {
    intx num_samples = surf.numSamples();
    for (intx i = 0; i < num_samples; ++i)
      callback(i, surf.getSample(i).getPosition());
  }
  else
  {
    Ball3 ball(position, max_distance);
    const_cast<PointCloud3::SampleKdTree &>(surf.getKdTree()).processRangeUntil<IntersectionTester>(ball, callback);
  }
}

// Called for each point in the geodesic neighborhood.
struct GeodesicCallback
{
  GeodesicCallback(Histogram & histogram_, Real acceptance_probability_)
  : histogram(histogram_), acceptance_probability(acceptance_probability_)
  {}

  bool operator()(PointCloud3::SampleGraph::VertexHandle vertex, double distance, bool has_pred,
                  PointCloud3::SampleGraph::VertexHandle pred)
  {
    if (acceptance_probability < 1 && Random::common().uniform01() > acceptance_probability)
      return false;

    histogram.insert(distance);

    return false;
  }

  Histogram & histogram;
  Real acceptance_probability;

}; // struct GeodesicCallback

void computeGeodesic(PointCloud3 const & surf, Vector3 const & position, Histogram & histogram, Real max_distance,
                     Real sample_reduction_ratio)
{
  if (sample_reduction_ratio < 0)
    sample_reduction_ratio = 1.1;  // play safe

  bool process_all = (max_distance < 0);
  if (process_all)
    max_distance = surf.getScale();

  histogram.setRange(0, std::max((double)max_distance, 1.0e-30));
  histogram.setZero();

  // Find the sample closest to the query position and use it as the source for all distance calculations
  intx seed_index =  surf.getKdTree().closestElement<MetricL2>(position);
  alwaysAssertM(seed_index >= 0, "LocalDistanceHistogram: Seed sample for geodesic distances not found");

  // Assume the graph and the kd-tree have samples in the same sequence
  PointCloud3::SampleGraph::VertexHandle seed_sample = const_cast<SamplePoint3 *>(&surf.getSample(seed_index));

  ShortestPaths<PointCloud3::SampleGraph> shortest_paths;
  GeodesicCallback callback(histogram, sample_reduction_ratio);
  shortest_paths.dijkstraWithCallback(const_cast<PointCloud3::SampleGraph &>(surf.getGraph()), seed_sample, callback,
                                      (process_all ? -1 : max_distance));
}

} // namespace LocalDistanceHistogramInternal

void
LocalDistanceHistogram::compute(Vector3 const & position, Histogram & histogram, DistanceType dist_type, Real max_distance,
                                Real sample_reduction_ratio) const
{
  switch (dist_type)
  {
    case DistanceType::EUCLIDEAN:
      LocalDistanceHistogramInternal::computeEuclidean(*surf, position, histogram, max_distance, sample_reduction_ratio);
      return;

    case DistanceType::GEODESIC:
      LocalDistanceHistogramInternal::computeGeodesic(*surf, position, histogram, max_distance, sample_reduction_ratio);
      return;

    default: throw Error("Unknown distance metric");
  }
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
