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

#ifndef __Thea_Algorithms_MeshFeatures_Local_AverageDistance_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_AverageDistance_hpp__

#include "../../../Common.hpp"
#include "../SampledSurface.hpp"
#include "../../IntersectionTester.hpp"
#include "../../MetricL2.hpp"
#include "../../PointTraitsN.hpp"
#include "../../SampleGraph.hpp"
#include "../../ShortestPaths.hpp"
#include "../../../Ball3.hpp"
#include "../../../Noncopyable.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/** Compute the average distance from a query point to other points on a shape. */
template < typename ExternalSampleKdTreeT = KdTreeN<Vector3, 3> >
class AverageDistance : public SampledSurface<ExternalSampleKdTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKdTreeT> BaseT;  ///< Base class.
    static intx const DEFAULT_NUM_SAMPLES = 5000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute the histogram of distances to sample points on a given mesh. The mesh need not persist
     * beyond the execution of the constructor.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used as the default way to normalize the average distance (an actual
     *   distance of \a normalization_scale will be mapped to 1). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    AverageDistance(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances to sample points on a given mesh group. The mesh group need
     * not persist beyond the execution of the constructor.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used as the default way to normalize the average distance (an actual
     *   distance of \a normalization_scale will be mapped to 1). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    AverageDistance(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances to sample points of a shape with a precomputed kd-tree on
     * these points. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ The kd-tree representing the shape.
     * @param normalization_scale The scale of the shape, used as the default way to normalize the average distance (an actual
     *   distance of \a normalization_scale will be mapped to 1). If <= 0, the bounding sphere diameter will be used.
     */
    AverageDistance(ExternalSampleKdTreeT const * sample_kdtree_, Real normalization_scale = -1)
    : BaseT(sample_kdtree_, normalization_scale)
    {}

    /**
     * Compute the average distance from a query point to sample points on the shape. The returned distance is normalized by
     * dividing by the normalization scale, which is \a max_distance (if non-negative), else as specified in the constructor.
     *
     * @param position The position of the query point.
     * @param dist_type The type of distance metric.
     * @param max_distance The distance to the furthest point to consider. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor).
     */
    double compute(Vector3 const & position, DistanceType dist_type = DistanceType::EUCLIDEAN, Real max_distance = -1) const
    {
      switch (dist_type)
      {
        case DistanceType::EUCLIDEAN: return computeEuclidean(position, max_distance);
        case DistanceType::GEODESIC:  return computeGeodesic(position, max_distance);
        default: throw Error("Unknown distance metric");
      }
    }

  private:
    /**
     * Compute the average euclidean distance from a query point to sample points on the shape. The returned distance is
     * normalized by dividing by the normalization scale, which is \a max_distance (if non-negative), else as specified in the
     * constructor.
     *
     * @param position The position of the query point.
     * @param max_distance The distance to the furthest point to consider. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor).
     */
    double computeEuclidean(Vector3 const & position, Real max_distance) const
    {
      bool process_all = (max_distance < 0);
      if (process_all)
        max_distance = this->getNormalizationScale();

      EuclideanCallback callback(position);

      if (process_all)
      {
        intx num_samples = this->numSamples();
        for (intx i = 0; i < num_samples; ++i)
          callback(i, this->getSamplePosition(i));
      }
      else
      {
        Ball3 ball(position, max_distance);

        if (this->hasExternalKdTree())
          this->getMutableExternalKdTree()->template processRangeUntil<IntersectionTester>(ball, std::ref(callback));
        else
          this->getMutableInternalKdTree()->template processRangeUntil<IntersectionTester>(ball, std::ref(callback));
      }

      return callback.getAverageDistance() / max_distance;
    }

    /**
     * Compute the average geodesic distance from a query point to sample points on the shape. The returned distance is
     * normalized by dividing by the normalization scale, which is \a max_distance (if non-negative), else as specified in the
     * constructor.
     *
     * @param position The position of the query point.
     * @param max_distance The distance to the furthest point to consider. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor).
     */
    double computeGeodesic(Vector3 const & position, Real max_distance) const
    {
      bool process_all = (max_distance < 0);
      if (process_all)
        max_distance = this->getNormalizationScale();

      SampleGraph * graph = const_cast<SampleGraph *>(this->getSampleGraph());
      alwaysAssertM(graph, "AverageDistance: Non-null sample graph required to compute geodesic distances");

      // Find the sample closest to the query position and use it as the source for all distance calculations
      intx seed_index = -1;
      if (this->hasExternalKdTree())
        seed_index = this->getMutableExternalKdTree()->template closestElement<MetricL2>(position);
      else
        seed_index = this->getMutableInternalKdTree()->template closestElement<MetricL2>(position);

      alwaysAssertM(seed_index >= 0, "AverageDistance: Seed sample for geodesic distances not found");

      // Assume the graph and the kd-tree have samples in the same sequence
      SampleGraph::SurfaceSample * seed_sample = const_cast<SampleGraph::SurfaceSample *>(&graph->getSample(seed_index));

      ShortestPaths<SampleGraph> shortest_paths;
      GeodesicCallback callback;
      shortest_paths.dijkstraWithCallback(*graph, seed_sample, std::ref(callback), (process_all ? -1 : max_distance));

      return callback.getAverageDistance() / max_distance;
    }

    /** Called for each point in the euclidean neighborhood. */
    struct EuclideanCallback : public Noncopyable
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

    /** Called for each point in the geodesic neighborhood. */
    struct GeodesicCallback : public Noncopyable
    {
      GeodesicCallback() : sum_distances(0), num_points(0) {}

      bool operator()(SampleGraph::VertexHandle vertex, double distance, bool has_pred, SampleGraph::VertexHandle pred)
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

}; // class AverageDistance

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
