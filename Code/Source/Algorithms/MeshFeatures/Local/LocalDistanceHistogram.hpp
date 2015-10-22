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
#include "../../../Graphics/MeshGroup.hpp"
#include "../../BestFitSphere3.hpp"
#include "../../IntersectionTester.hpp"
#include "../../KDTreeN.hpp"
#include "../../MeshSampler.hpp"
#include "../../PointCollectorN.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Ball3.hpp"
#include "../../../Math.hpp"
#include "../../../Random.hpp"
#include "../../../Vector3.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/** Compute the histogram of distances from a query point to other points on a shape. */
template < typename ExternalSampleKDTreeT = KDTreeN<Vector3, 3> >
class LocalDistanceHistogram
{
  public:
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    typedef KDTreeN<Vector3, 3> SampleKDTree;  ///< A kd-tree on mesh samples.

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
    : sample_kdtree(NULL), precomp_kdtree(NULL), scale(normalization_scale)
    {
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      MeshSampler<MeshT> sampler(mesh);
      sampler.sampleEvenlyByArea(num_samples, samples);

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh);
        scale = bsphere.getDiameter();
      }
    }

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
    : sample_kdtree(NULL), precomp_kdtree(NULL), scale(normalization_scale)
    {
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      MeshSampler<MeshT> sampler(mesh_group);
      sampler.sampleEvenlyByArea(num_samples, samples);

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh_group);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Constructs the object to compute the histogram of distances to sample points of a shape with a precomputed kd-tree on
     * these points. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ The kd-tree representing the shape.
     * @param normalization_scale The scale of the shape, used to define the size of histogram bins if the latter is not
     *   explicitly specified when calling compute(). If <= 0, the bounding sphere diameter will be used.
     */
    LocalDistanceHistogram(ExternalSampleKDTree const * sample_kdtree_, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(sample_kdtree_), scale(normalization_scale)
    {
      alwaysAssertM(precomp_kdtree, "LocalDistanceHistogram: Precomputed KD-tree cannot be null");

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(precomp_kdtree->getElements(),
                                                               precomp_kdtree->getElements() + precomp_kdtree->numElements());
        scale = bsphere.getDiameter();
      }
    }

    /** Destructor. */
    ~LocalDistanceHistogram()
    {
      delete sample_kdtree;
    }

    /** Get the number of surface samples. */
    long numSamples() const
    {
      if (precomp_kdtree)
        return precomp_kdtree->numElements();
      else
        return (long)samples.size();
    }

    /** Get the position of the surface sample with index \a index. */
    Vector3 getSamplePosition(long index) const
    {
      debugAssertM(index >= 0 && index < numSamples(), format("LocalDistanceHistogram: Sample index %ld out of bounds", index));

      if (precomp_kdtree)
      {
        typedef typename ExternalSampleKDTree::Element ExternalSample;
        return PointTraitsN<ExternalSample, 3>::getPosition(precomp_kdtree->getElements()[(array_size_t)index]);
      }
      else
        return samples[(array_size_t)index];
    }

    /**
     * Compute the histogram of distances from a query point to sample points on the shape. The histogram bins uniformly
     * subdivide the range of distances from zero to \a max_distance.
     *
     * @param position The position of the query point.
     * @param num_bins The number of bins in the output histogram.
     * @param histogram An array preallocated to \a num_bins entries, each entry corresponding to a histogram bin.
     * @param max_distance The distance to the furthest point to consider for the histogram. A negative value indicates the
     *   entire shape is to be considered (in which case \a max_distance is set to the shape scale specified in the
     *   constructor).
     * @param sample_reduction_ratio The fraction of the available set of samples -- expressed as a number between 0 and 1 --
     *   that will be randomly selected and used to actually build the histogram. This may be useful for getting a more evenly
     *   sampled set of pairwise distances when calling this function with multiple query points (and an extra-large initial set
     *   of points). A negative value, or a value of 1, indicates all sample points will be used.
     *
     * @return The sum of the values in the histogram bins. Can be used to normalize the histogram, especially if the exact
     *   number of pairs used to build the histogram is unknown.
     */
    double compute(Vector3 const & position, long num_bins, double * histogram, Real max_distance = -1,
                   Real sample_reduction_ratio = -1) const
    {
      if (sample_reduction_ratio < 0)
        sample_reduction_ratio = 1.1;  // play safe

      double hist_sum = 0;

      if (max_distance < 0)
      {
        max_distance = scale;
        LocalDistanceHistogramFunctor func(position, num_bins, histogram, max_distance, sample_reduction_ratio);

        if (precomp_kdtree)
        {
          typename ExternalSampleKDTree::value_type const * elems = precomp_kdtree->getElements();
          long num_elems = precomp_kdtree->numElements();

          for (long i = 0; i < num_elems; ++i)
            func(i, elems[i]);
        }
        else
        {
          for (array_size_t i = 0; i < samples.size(); ++i)
            func((long)i, samples[i]);
        }

        hist_sum = func.sum_values;
      }
      else
      {
        if (!precomp_kdtree && !sample_kdtree)
          sample_kdtree = new SampleKDTree(samples.begin(), samples.end());

        LocalDistanceHistogramFunctor func(position, num_bins, histogram, max_distance, sample_reduction_ratio);
        Ball3 ball(position, max_distance);

        if (precomp_kdtree)
          const_cast<ExternalSampleKDTree *>(precomp_kdtree)->template processRangeUntil<IntersectionTester>(ball, &func);
        else
          sample_kdtree->template processRangeUntil<IntersectionTester>(ball, &func);

        hist_sum = func.sum_values;
      }

      return hist_sum;
    }

  private:
    /** Called for each point in the neighborhood. */
    struct LocalDistanceHistogramFunctor
    {
      LocalDistanceHistogramFunctor(Vector3 const & position_, long num_bins_, double * histogram_, Real max_distance_,
                                    Real acceptance_probability_)
      : position(position_), num_bins(num_bins_), histogram(histogram_), bin_scale(0),
        acceptance_probability(acceptance_probability_), sum_values(0)
      {
        alwaysAssertM(num_bins_ > 0, "LocalDistanceHistogram: Number of bins must be positive");
        alwaysAssertM(max_distance_ >= 0, "LocalDistanceHistogram: Maximum distance must be non-negative");

        bin_scale = std::max(max_distance_ / num_bins_, 1.0e-30f);

        for (long i = 0; i < num_bins_; ++i)
          histogram_[i] = 0.0;
      }

      template <typename SampleT> bool operator()(long index, SampleT & t)
      {
        if (acceptance_probability < 1 && Random::common().uniform01() > acceptance_probability)
          return false;

        Real d = (PointTraitsN<SampleT, 3>::getPosition(t) - position).length();
        long bin = Math::clamp((long)std::floor(d / bin_scale), 0, num_bins - 1);
        histogram[bin] += 1.0;
        sum_values += 1.0;

        return false;
      }

      Vector3 position;
      long num_bins;
      double * histogram;
      Real bin_scale;
      Real acceptance_probability;
      double sum_values;

    }; // struct LocalDistanceHistogramFunctor

    TheaArray<Vector3> samples;  ///< MeshT sample points computed by this object.
    mutable SampleKDTree * sample_kdtree;  ///< KD-tree on mesh samples.
    ExternalSampleKDTree const * precomp_kdtree;  ///< Precomputed KD-tree on mesh samples.
    Real scale;  ///< The normalization length.

}; // class LocalDistanceHistogram

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
