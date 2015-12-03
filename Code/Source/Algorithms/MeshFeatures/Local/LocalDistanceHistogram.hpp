//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_MeshFeatures_Local_LocalDistanceHistogram_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_LocalDistanceHistogram_hpp__

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

/** Compute the histogram of distances from a query point to other points on a shape. */
template < typename ExternalSampleKDTreeT = KDTreeN<Vector3, 3> >
class LocalDistanceHistogram : public SampledSurface<ExternalSampleKDTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKDTreeT> BaseT;  ///< Base class.
    static long const DEFAULT_NUM_SAMPLES = 5000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute the histogram of distances to sample points on a given mesh. The mesh need not persist
     * beyond the execution of the constructor.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    LocalDistanceHistogram(MeshT const & mesh, long num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances to sample points on a given mesh group. The mesh group need
     * not persist beyond the execution of the constructor.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    template <typename MeshT>
    LocalDistanceHistogram(Graphics::MeshGroup<MeshT> const & mesh_group, long num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute the histogram of distances to sample points of a shape with a precomputed kd-tree on
     * these points. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ The kd-tree representing the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    LocalDistanceHistogram(ExternalSampleKDTreeT const * sample_kdtree_, Real normalization_scale = -1)
    : BaseT(sample_kdtree_, normalization_scale)
    {}

    /**
     * Compute the histogram of distances from a query point to sample points on the shape. The histogram bins uniformly
     * subdivide the range of distances from zero to \a max_distance. If \a max_distance is negative, the shape scale specified
     * in the constructor will be used.
     *
     * @param position The position of the query point.
     * @param histogram The histogram to be computed.
     * @param dist_type The type of distance metric.
     * @param max_distance The distance to the furthest point to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor). The histogram range is set appropriately.
     * @param sample_reduction_ratio The fraction of the available set of samples -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. This may be useful for getting a more evenly
     *   sampled set of pairwise distances when calling this function with multiple query points (and an extra-large initial set
     *   of points). A negative value, or a value of 1, indicates all sample points will be used.
     */
    void compute(Vector3 const & position, Histogram & histogram, DistanceType dist_type = DistanceType::EUCLIDEAN,
                 Real max_distance = -1, Real sample_reduction_ratio = -1) const
    {
      switch (dist_type)
      {
        case DistanceType::EUCLIDEAN: computeEuclidean(position, histogram, max_distance, sample_reduction_ratio); return;
        case DistanceType::GEODESIC:  computeGeodesic(position, histogram, max_distance, sample_reduction_ratio); return;
        default: throw Error("Unknown distance metric");
      }
    }

  private:
    /**
     * Compute the histogram of euclidean distances from a query point to sample points on the shape. The histogram bins
     * uniformly subdivide the range of distances from zero to \a max_distance. If \a max_distance is negative, the shape scale
     * specified in the constructor will be used.
     *
     * @param position The position of the query point.
     * @param histogram The histogram to be computed.
     * @param max_distance The distance to the furthest point to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor). The histogram range is set appropriately.
     * @param sample_reduction_ratio The fraction of the available set of samples -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. This may be useful for getting a more evenly
     *   sampled set of pairwise distances when calling this function with multiple query points (and an extra-large initial set
     *   of points). A negative value, or a value of 1, indicates all sample points will be used.
     */
    void computeEuclidean(Vector3 const & position, Histogram & histogram, Real max_distance, Real sample_reduction_ratio) const
    {
      if (sample_reduction_ratio < 0)
        sample_reduction_ratio = 1.1;  // play safe

      bool process_all = (max_distance < 0);
      if (process_all)
        max_distance = this->getNormalizationScale();

      histogram.setRange(0, std::max((double)max_distance, 1.0e-30));
      histogram.setZero();

      EuclideanCallback callback(position, histogram, sample_reduction_ratio);

      if (process_all)
      {
        long num_samples = this->numSamples();
        for (long i = 0; i < num_samples; ++i)
          callback(i, this->getSamplePosition(i));
      }
      else
      {
        Ball3 ball(position, max_distance);

        if (this->hasExternalKDTree())
          this->getMutableExternalKDTree()->template processRangeUntil<IntersectionTester>(ball, &callback);
        else
          this->getMutableInternalKDTree()->template processRangeUntil<IntersectionTester>(ball, &callback);
      }
    }

    /**
     * Compute the histogram of euclidean distances from a query point to sample points on the shape. The histogram bins
     * uniformly subdivide the range of distances from zero to \a max_distance. If \a max_distance is negative, the shape scale
     * specified in the constructor will be used.
     *
     * @param position The position of the query point.
     * @param histogram The histogram to be computed.
     * @param max_distance The distance to the furthest point to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor). The histogram range is set appropriately.
     * @param sample_reduction_ratio The fraction of the available set of samples -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. This may be useful for getting a more evenly
     *   sampled set of pairwise distances when calling this function with multiple query points (and an extra-large initial set
     *   of points). A negative value, or a value of 1, indicates all sample points will be used.
     */
    void computeGeodesic(Vector3 const & position, Histogram & histogram, Real max_distance, Real sample_reduction_ratio) const
    {
      if (sample_reduction_ratio < 0)
        sample_reduction_ratio = 1.1;  // play safe

      bool process_all = (max_distance < 0);
      if (process_all)
        max_distance = this->getNormalizationScale();

      histogram.setRange(0, std::max((double)max_distance, 1.0e-30));
      histogram.setZero();

      SampleGraph * graph = const_cast<SampleGraph *>(this->getSampleGraph());
      alwaysAssertM(graph, "LocalDistanceHistogram: Non-null sample graph required to compute geodesic distances");

      // Find the sample closest to the query position and use it as the source for all distance calculations
      long seed_index = -1;
      if (this->hasExternalKDTree())
        seed_index = this->getMutableExternalKDTree()->template closestElement<MetricL2>(position);
      else
        seed_index = this->getMutableInternalKDTree()->template closestElement<MetricL2>(position);

      alwaysAssertM(seed_index >= 0, "LocalDistanceHistogram: Seed sample for geodesic distances not found");

      // Assume the graph and the kd-tree have samples in the same sequence
      SampleGraph::SurfaceSample * seed_sample = const_cast<SampleGraph::SurfaceSample *>(&graph->getSample(seed_index));

      ShortestPaths<SampleGraph> shortest_paths;
      GeodesicCallback callback(histogram, sample_reduction_ratio);
      shortest_paths.dijkstraWithCallback(*graph, seed_sample, &callback, (process_all ? -1 : max_distance));
    }

    /** Called for each point in the euclidean neighborhood. */
    struct EuclideanCallback
    {
      EuclideanCallback(Vector3 const & position_, Histogram & histogram_, Real acceptance_probability_)
      : position(position_), histogram(histogram_), acceptance_probability(acceptance_probability_)
      {}

      template <typename SampleT> bool operator()(long index, SampleT const & t)
      {
        if (acceptance_probability < 1 && Random::common().uniform01() > acceptance_probability)
          return false;

        Real d = (PointTraitsN<SampleT, 3>::getPosition(t) - position).length();
        histogram.insert(d);

        return false;
      }

      Vector3 position;
      Histogram & histogram;
      Real acceptance_probability;

    }; // struct EuclideanCallback

    /** Called for each point in the geodesic neighborhood. */
    struct GeodesicCallback
    {
      GeodesicCallback(Histogram & histogram_, Real acceptance_probability_)
      : histogram(histogram_), acceptance_probability(acceptance_probability_)
      {}

      bool operator()(SampleGraph::VertexHandle vertex, double distance, bool has_pred, SampleGraph::VertexHandle pred)
      {
        if (acceptance_probability < 1 && Random::common().uniform01() > acceptance_probability)
          return false;

        histogram.insert(distance);

        return false;
      }

      Histogram & histogram;
      Real acceptance_probability;

    }; // struct GeodesicCallback

}; // class LocalDistanceHistogram

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
