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

#include "AverageDistance.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../PointTraitsN.hpp"
#include "../../ShortestPaths.hpp"
#include "../../../Ball3.hpp"
#include "../../../Noncopyable.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

AverageDistance::AverageDistance(PointSet3 const * surf_)
: surf(surf_)
{
  alwaysAssertM(surf_, "AverageDistance: Cannot construct with a null surface");
}

namespace AverageDistanceInternal {

// Called for each point in the euclidean neighborhood.
struct EuclideanCallback : private Noncopyable
{
  EuclideanCallback(Vector3 const & position_) : position(position_), sum_distances(0), num_points(0) {}

  template <typename SampleT> bool operator()(intx index, SampleT const & t)
  {
    Real d = (PointTraitsN<SampleT, 3>::getPosition(t) - position).norm();
    sum_distances += d;
    num_points++;

    return false;
  }

  double getAverageDistance() const
  {
    return num_points <= 0 ? 0.0 : sum_distances / num_points;
  }

  Vector3 position;
  double sum_distances;
  intx num_points;

}; // struct EuclideanCallback

double
computeEuclidean(PointSet3 const & surf, Vector3 const & position, Real max_distance)
{
  bool process_all = (max_distance < 0);
  if (process_all)
    max_distance = surf.getScale();

  EuclideanCallback callback(position);

  if (process_all)
  {
    intx num_samples = surf.numSamples();
    for (intx i = 0; i < num_samples; ++i)
      callback(i, surf.getSample(i).getPosition());
  }
  else
  {
    Ball3 ball(position, max_distance);
    const_cast<PointSet3::SampleKdTree &>(surf.getKdTree()).processRangeUntil<IntersectionTester>(ball, std::ref(callback));
  }

  return callback.getAverageDistance() / max_distance;
}

// Called for each point in the geodesic neighborhood.
struct GeodesicCallback : private Noncopyable
{
  GeodesicCallback() : sum_distances(0), num_points(0) {}

  bool operator()(PointSet3::SampleGraph::VertexHandle vertex, double distance, bool has_pred,
                  PointSet3::SampleGraph::VertexHandle pred)
  {
    sum_distances += distance;
    num_points++;

    return false;
  }

  double getAverageDistance() const
  {
    return num_points <= 0 ? 0.0 : sum_distances / num_points;
  }

  double sum_distances;
  intx num_points;

}; // struct GeodesicCallback

double
computeGeodesic(PointSet3 const & surf, Vector3 const & position, Real max_distance)
{
  bool process_all = (max_distance < 0);
  if (process_all)
    max_distance = surf.getScale();

  // Find the sample closest to the query position and use it as the source for all distance calculations
  intx seed_index =  surf.getKdTree().closestElement<MetricL2>(position);
  alwaysAssertM(seed_index >= 0, "AverageDistance: Seed sample for geodesic distances not found");

  // Assume the graph and the kd-tree have samples in the same sequence
  PointSet3::SampleGraph::VertexHandle seed_sample = const_cast<SamplePoint3 *>(&surf.getSample(seed_index));

  ShortestPaths<PointSet3::SampleGraph> shortest_paths;
  GeodesicCallback callback;
  shortest_paths.dijkstraWithCallback(const_cast<PointSet3::SampleGraph &>(surf.getGraph()), seed_sample, std::ref(callback),
                                      (process_all ? -1 : max_distance));

  return callback.getAverageDistance() / max_distance;
}

} // namespace AverageDistanceInternal

double
AverageDistance::compute(Vector3 const & position, DistanceType dist_type, Real max_distance) const
{
  switch (dist_type)
  {
    case DistanceType::EUCLIDEAN: return AverageDistanceInternal::computeEuclidean(*surf, position, max_distance);
    case DistanceType::GEODESIC:  return AverageDistanceInternal::computeGeodesic(*surf, position, max_distance);
    default: throw Error("Unknown distance metric");
  }
}

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea
