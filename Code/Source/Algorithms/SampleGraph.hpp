//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Cornell University
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

#ifndef __Thea_Algorithms_SampleGraph_hpp__
#define __Thea_Algorithms_SampleGraph_hpp__

#include "../Common.hpp"
#include "Filter.hpp"
#include "KDTreeN.hpp"
#include "IntersectionTester.hpp"
#include "MetricL2.hpp"
#include "PointTraitsN.hpp"
#include "RayIntersectionTester.hpp"
#include "RayQueryStructureN.hpp"
#include "../BoundedSortedArray.hpp"
#include "../Noncopyable.hpp"
#include "../Vector3.hpp"
#include <boost/type_traits/has_trivial_assign.hpp>

namespace Thea {
namespace Algorithms {

namespace SampleGraphInternal {

// Forward declarations
class SurfaceSample;

/** A sample neighboring another sample. */
class NeighboringSample
{
  public:
    /** Constructor. */
    NeighboringSample(SurfaceSample * sample_ = NULL, Real separation_ = 0) : sample(sample_), separation(separation_) {}

    /** Get a reference to this sample. */
    SurfaceSample * getSample() const { return sample; }

    /** Set the reference to the neighboring sample. */
    void setSample(SurfaceSample * sample_) { sample = sample_; }

    /** Get the separation of this sample from the source sample. */
    Real getSeparation() const { return separation; }

    /** Set the separation of this sample from the source sample. */
    void setSeparation(Real separation_) { separation = separation_; }

    /** Compare neighbors by their separation from the source. */
    bool operator<(NeighboringSample const & rhs) const
    {
      return separation < rhs.separation || (separation == rhs.separation && std::less<SurfaceSample *>()(sample, rhs.sample));
    }

  private:
    SurfaceSample * sample;  ///< Reference to this sample.
    Real separation;         ///< Separation from the source sample.

}; // class NeighboringSample

} // namespace SampleGraphInternal
} // namespace Algorithms
} // namespace Thea

namespace boost {

template <> struct has_trivial_assign<Thea::Algorithms::SampleGraphInternal::NeighboringSample> : public boost::true_type {};

} // namespace boost

namespace Thea {
namespace Algorithms {
namespace SampleGraphInternal {

/**
 * A class used to store surface samples with adjacency information. Encapsulates a sample position, normal and list of
 * neighboring samples.
 */
class SurfaceSample
{
  public:
    typedef NeighboringSample Neighbor;                ///< A neighboring sample.
    typedef BoundedSortedArray<Neighbor> NeighborSet;  ///< Set of neighboring samples, sorted by distance from this sample.

    /** Default constructor. */
    SurfaceSample() {}

    /** Construct with an index and maximum number of neighbors. */
    SurfaceSample(long index_, int max_nbrs) : index(index_), nbrs(max_nbrs) {}

    /** Get the sample index. */
    long getIndex() const { return index; }

    /** Set the sample index. */
    void setIndex(long index_) { index = index_; }

    /** Get the position of the sample. */
    Vector3 const & getPosition() const { return p; }

    /** Set the position of the sample. */
    void setPosition(Vector3 const & p_) { p = p_; }

    /** Get the surface normal at the sample position. */
    Vector3 const & getNormal() const { return n; }

    /** Set the surface normal at the sample position. */
    void setNormal(Vector3 const & n_) { n = n_; }

    /** Get the number of neighboring samples. */
    int numNeighbors() const { return nbrs.size(); }

    /** Get the set of neighboring samples. */
    NeighborSet const & getNeighbors() const { return nbrs; }

    /** Get a non-const reference to the set of neighboring samples. */
    NeighborSet & getNeighbors() { return nbrs; }

  private:
    long index;        ///< Index of the sample.
    Vector3 p;         ///< Sample position.
    Vector3 n;         ///< Sample normal.
    NeighborSet nbrs;  ///< Set of neighboring samples.

}; // class SurfaceSample

} // namespace SampleGraphInternal

template <>
class IsPointN<SampleGraphInternal::SurfaceSample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<SampleGraphInternal::SurfaceSample, 3>::getPosition(SampleGraphInternal::SurfaceSample const & sample)
{
  return sample.getPosition();
}

/**
 * Proximity graph on surface samples. Satisfies the IsAdjacencyGraph concept. Typical usage is to first specify the sample
 * positions (and optionally normals) via setSamples(). Then, specify an optional dense oversampling via setOversampling() --
 * the oversampling makes it easier to verify if two samples are actually adjacent on the surface. Finally, call init(),
 * optionally with a representation of the underlying surface, to compute the sample adjacencies (graph edges).
 *
 * @note The graph is <b>not created</b> until init() is called.
 */
class SampleGraph : private Noncopyable
{
  private:
    /** Dummy surface representation used to disambiguate specialization of init(). */
    struct DummyRayQueryStructure3 : public RayQueryStructureN<3, Real>
    {
      typedef RayQueryStructureN<3, Real> BaseT;

      template <typename RayIntersectionTesterT> bool rayIntersects(RayT const & ray, Real max_time = -1) const
      { return false; }

      template <typename RayIntersectionTesterT> Real rayIntersectionTime(RayT const & ray, Real max_time = -1) const
      { return -1; }

      template <typename RayIntersectionTesterT>
      typename BaseT::RayStructureIntersectionT rayStructureIntersection(RayT const & ray, Real max_time = -1) const
      { return typename BaseT::RayStructureIntersectionT(); }

    }; // struct DummyRayQueryStructure3

  public:
    /** %Options for constructing the sample graph. */
    class Options
    {
      public:
        /** Set the maximum degree of the graph. */
        Options & setMaxDegree(int max_degree_) { max_degree = max_degree_; return *this; }

        /** Construct with default values. */
        Options() : max_degree(8) {}

        /** Get a set of options with default values. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        int max_degree;  ///< Maximum degree of the graph.

        friend class SampleGraph;

    }; // class Options

    typedef SampleGraphInternal::SurfaceSample  SurfaceSample;          ///< Point sampled from a surface.
    typedef TheaArray<SurfaceSample>            SampleArray;            ///< Array of samples.

  private:
    typedef KDTreeN<SurfaceSample *, 3>         SampleKDTree;           ///< KD-tree on surface samples.

  public:
    //=========================================================================================================================
    //
    // Types and functions required to make this a bona fide graph satisfying the IsAdjacencyGraph concept.
    //
    //=========================================================================================================================

    typedef SurfaceSample *                     VertexHandle;           ///< Handle to a graph vertex (a sample).
    typedef SurfaceSample const *               VertexConstHandle;      ///< Const handle to a graph vertex (a sample).
    typedef SampleArray::iterator               VertexIterator;         ///< Iterator over vertices.
    typedef SampleArray::const_iterator         VertexConstIterator;    ///< Const iterator over vertices.
    typedef SurfaceSample::Neighbor *           NeighborIterator;       ///< Iterator over neighbors of a vertex.
    typedef SurfaceSample::Neighbor const *     NeighborConstIterator;  ///< Const iterator over neighbors of a vertex.

    /** Get the number of vertices (samples) in the graph. */
    long numVertices() const { return (long)samples.size(); }

    /** Get an iterator to the first vertex. */
    VertexIterator verticesBegin() { return samples.begin(); }

    /** Get a const iterator to the first vertex. */
    VertexConstIterator verticesBegin() const { return samples.begin(); }

    /** Get an iterator to one position beyond the last vertex. */
    VertexIterator verticesEnd() { return samples.end(); }

    /** Get a const iterator to one position beyond the last vertex. */
    VertexConstIterator verticesEnd() const { return samples.end(); }

    /** Get a handle to the vertex referenced by an iterator. */
    VertexHandle getVertex(VertexIterator vi) { return &(*vi); }

    /** Get a handle to the vertex referenced by a const iterator. */
    VertexConstHandle getVertex(VertexConstIterator vi) const { return &(*vi); }

    /** Get the number of neighbors of a vertex. */
    long numNeighbors(VertexConstHandle vertex) const { return vertex->getNeighbors().size(); }

    /** Get an iterator to the first neighbor of a vertex. */
    NeighborIterator neighborsBegin(VertexHandle vertex)
    {
      return vertex->getNeighbors().isEmpty() ? NULL : const_cast<SurfaceSample::Neighbor *>(&vertex->getNeighbors()[0]);
    }

    /** Get a const iterator to the first neighbor of a vertex. */
    NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
    {
      return vertex->getNeighbors().isEmpty() ? NULL : &vertex->getNeighbors()[0];
    }

    /** Get an iterator to the one position beyond the last neighbor of a vertex. */
    NeighborIterator neighborsEnd(VertexHandle vertex) { return neighborsBegin(vertex) + numNeighbors(vertex); }

    /** Get a const iterator to the one position beyond the last neighbor of a vertex. */
    NeighborConstIterator neighborsEnd(VertexConstHandle vertex) const { return neighborsBegin(vertex) + numNeighbors(vertex); }

    /** Get a handle to the neighboring vertex referenced by an iterator. */
    VertexHandle getVertex(NeighborIterator ni) { return ni->getSample(); }

    /** Get a handle to the neighboring vertex referenced by a const iterator. */
    VertexConstHandle getVertex(NeighborConstIterator ni) const { return ni->getSample(); }

    /** Get the distance between a vertex and its neighbor. */
    double distance(VertexConstHandle v, NeighborConstIterator ni) const { return ni->getSeparation(); }

    //=========================================================================================================================
    //
    // Functions to create the graph from a set of samples and an optional base shape.
    //
    //=========================================================================================================================

    /** Construct with a set of options for creating the graph. */
    SampleGraph(Options const & options_ = Options::defaults())
    : options(options_), has_normals(false), avg_separation(0), initialized(false)
    {
      alwaysAssertM(options.max_degree > 0, "SampleGraph: Maximum degree must be positive");
    }

    /** Copy constructor. */
    SampleGraph(SampleGraph const & src);

    /** Assignment operator. */
    SampleGraph & operator=(SampleGraph const & src);

    /** Get the number of samples in the graph. */
    long numSamples() const { return (long)samples.size(); }

    /** Get the set of samples in the graph. */
    SampleArray const & getSamples() const { return samples; }

    /** Get a sample by its index. */
    SurfaceSample const & getSample(long index) const
    {
      debugAssertM(index >= 0 && index < (long)samples.size(), "SampleGraph: Sample index out of bounds");
      return samples[(size_t)index];
    }

    /** Set the sample positions and (optionally) normals. All prior samples will be cleared. */
    void setSamples(long num_samples, Vector3 const * positions, Vector3 const * normals = NULL)
    {
      alwaysAssertM(num_samples >= 0, "SampleGraph: Cannot specify a negative number of samples");
      alwaysAssertM(num_samples == 0 || positions, "SampleGraph: Sample positions must be specified");

      samples.resize((size_t)num_samples);
      for (size_t i = 0; i < samples.size(); ++i)
      {
        samples[i].setIndex((long)i);
        samples[i].setPosition(positions[i]);
        samples[i].getNeighbors().setCapacity(options.max_degree);
      }

      if (normals && num_samples > 0)
      {
        for (size_t i = 0; i < samples.size(); ++i)
          samples[i].setNormal(normals[i]);

        has_normals = true;
      }
      else
        has_normals = false;

      initialized = false;
    }

    /**
     * Specify a previously computed oversampling of the surface. These additional samples will not constitute vertices of the
     * final graph, but will be used to more accurately compute adjacencies. If normals are not specified, they must also not
     * have been specified for the main samples.
     */
    void setOversampling(long num_samples, Vector3 const * dense_positions, Vector3 const * dense_normals = NULL)
    {
      alwaysAssertM(num_samples >= 0, "SampleGraph: Cannot specify a negative number of dense samples");
      alwaysAssertM(!has_normals || dense_normals,
                    "SampleGraph: Main samples have normals, so oversampling must also have normals");

      dense_samples.resize((size_t)num_samples);
      for (size_t i = 0; i < dense_samples.size(); ++i)
      {
        dense_samples[i].setPosition(dense_positions[i]);
        dense_samples[i].getNeighbors().setCapacity(options.max_degree);
      }

      if (dense_normals && has_normals)
      {
        for (size_t i = 0; i < dense_samples.size(); ++i)
          dense_samples[i].setNormal(dense_normals[i]);
      }

      initialized = false;
    }

    /** Construct the sample graph from a set of samples, without tests that require access to the underlying surface. */
    void init()
    {
      init<DummyRayQueryStructure3>(NULL);
    }

    /**
     * Construct from a set of sample positions and normals, and a representation of the underlying surface (such as a kd-tree
     * on mesh triangles). All possible adjacency tests, including those that test if the two samples are on opposite sides of a
     * thin surface, will be performed by this constructor.
     */
    template <typename RayQueryStructureT> void init(RayQueryStructureT const * surface)
    {
      if (initialized)  // don't initialize twice
        return;

      // Aggregate samples
      TheaArray<SurfaceSample *> sample_ptrs(samples.size() + dense_samples.size());

      for (size_t i = 0; i < samples.size(); ++i)
        sample_ptrs[i] = &samples[i];

      for (size_t i = 0; i < dense_samples.size(); ++i)
      {
        dense_samples[i].setIndex((long)(samples.size() + i));
        sample_ptrs[samples.size() + i] = &dense_samples[i];
      }

      // Clear any prior adjacency data
      for (size_t i = 0; i < sample_ptrs.size(); ++i)
        sample_ptrs[i]->getNeighbors().clear();

      if (sample_ptrs.size() <= 1)  // nothing to do
        return;

      // Create kd-tree on samples
      SampleKDTree sample_kdtree(sample_ptrs.begin(), sample_ptrs.end());

      // Get a measure of the average pairwise separation of samples
      avg_separation = 0;
      long num_trials = std::min(100L, (long)sample_ptrs.size());
      for (long i = 0; i < num_trials; ++i)
      {
        size_t index = (size_t)Random::common().integer(0, (int32)sample_ptrs.size() - 1);
        FilterSelf filter(sample_ptrs[index]);
        sample_kdtree.pushFilter(&filter);
          long nn_index = sample_kdtree.closestElement<MetricL2>(sample_ptrs[index]->getPosition());
          alwaysAssertM(nn_index >= 0, "SampleGraph: Nearest neighbor of sample not found");
          avg_separation += (sample_ptrs[(size_t)nn_index]->getPosition()
                           - sample_ptrs[index]->getPosition()).squaredLength();
        sample_kdtree.popFilter();
      }

      avg_separation = std::sqrt(avg_separation / num_trials);  // RMS average

      // Find the neighbors of each sample
      Real sep_scale = std::sqrt((Real)options.max_degree);  // assume uniform distribution on 2D surface
      for (size_t i = 0; i < sample_ptrs.size(); ++i)
        findSampleNeighbors(sample_ptrs[i], sample_kdtree, sep_scale * avg_separation, surface);

      // Extract adjacencies between original set of samples
      if (!dense_samples.empty())
        extractOriginalAdjacencies(sample_ptrs);

      // Update the average separation to the correct value
      updateAverageSeparation();

      initialized = true;
    }

    /**
     * Get the average separation of samples from their neighbors, that is, the average edge length of the graph. This may be
     * computed by sampling and hence should not be assumed to be a true average, just a representative value. Has a meaningful
     * value <b>only after graph initialization</b>, i.e. calling init(). It should not be called before then.
     */
    Real getAverageSeparation() const { return avg_separation; }

    /** Clear the graph. */
    void clear();

    /** Load the graph and samples from disk files. */
    bool load(std::string const & graph_path, std::string const & samples_path);

    /**
     * Save the graph (and optionally the samples) to disk.
     *
     * @param graph_path Output file for graph.
     * @param samples_path Output file for samples.
     * @param write_distances If true, the distance of each neighbor (which may be different from the Euclidean distance if
     *   the graph was computed via an oversampling) is also written to the file.
     */
    bool save(std::string const & graph_path, std::string const & samples_path = "", bool write_distances = false) const;

  private:
    /** Allows every sample except one. */
    struct FilterSelf : public Filter<SurfaceSample *>
    {
      FilterSelf(SurfaceSample * self_) : self(self_) {}
      bool allows(SurfaceSample * const & s) const { return s != self; }

      SurfaceSample * self;

    }; // struct FilterSelf

    /** Functor that validates and adds neighbors to a given sample. */
    template <typename RayQueryStructureT> class NeighborFunctor
    {
      public:
        /** Constructor. */
        NeighborFunctor(SurfaceSample * sample_, bool has_normals_, RayQueryStructureT const * surface_)
        : sample(sample_), has_normals(has_normals_), surface(surface_)
        {
          debugAssertM(sample, "SampleGraph: Can't create neighbor functor without valid source sample");
        }

        /** Called for every candidate neighbor. */
        bool operator()(long index, SurfaceSample * nbr)
        {
          Real sep = 0;
          if (isValidNeighbor(nbr, sep))
            sample->getNeighbors().insert(SurfaceSample::Neighbor(nbr, sep));

          return false;
        }

      private:
        /**
         * Check if a candidate neighboring sample is a valid neighbor, that is, if it is near the original sample as measured
         * on the surface.
         */
        bool isValidNeighbor(SurfaceSample const * nbr, Real & sep)
        {
          // Tests are (hopefully) in decreasing order of speed

          // Identity test
          if (sample == nbr)
            return false;

          // Angle test
          if (has_normals)
          {
            static Real const MIN_DOT = -0.5f;
            if (sample->getNormal().dot(nbr->getNormal()) < MIN_DOT)
              return false;
          }

          // Duplication test which compares only the samples themselves (since numerical error can cause multiple attempts to
          // insert the same sample with slightly different separations)
          for (int i = 0; i < sample->numNeighbors(); ++i)
            if (nbr == sample->getNeighbors()[i].getSample())
              return false;

          // Compute the separation of the samples
          Vector3 diff = nbr->getPosition() - sample->getPosition();
          Real sqsep = diff.squaredLength();
          sep = Math::fastSqrt(sqsep);

          // Reachability test: lift points off the surface by an amount proportional to their separation, and see if the line
          // connecting them is blocked by the surface
          if (surface && has_normals)
          {
            static Real const LIFT_FACTOR = 5;
            Vector3 lift_dir = (sample->getNormal() + nbr->getNormal()).fastUnit();
            typename RayQueryStructureT::RayT ray(sample->getPosition() + LIFT_FACTOR * sep * lift_dir, diff);
            if (surface->template rayIntersects<RayIntersectionTester>(ray, 1))
              return false;
          }

          return true;
        }

        SurfaceSample * sample;
        bool has_normals;
        RayQueryStructureT const * surface;

    }; // struct NeighborFunctor

    /** Find the samples adjacent to a given sample. */
    template <typename RayQueryStructureT>
    void findSampleNeighbors(SurfaceSample * sample, SampleKDTree & sample_kdtree, Real init_radius,
                             RayQueryStructureT const * surface)
    {
      static int const MAX_ITERS = 3;
      static Real const RADIUS_EXPANSION_FACTOR = 2.0f;
      int min_degree = Math::clamp((int)(0.25 * options.max_degree), 4, options.max_degree);

      Real radius = init_radius;
      for (int i = 0; i < MAX_ITERS; ++i)
      {
        sample->getNeighbors().clear();
        NeighborFunctor<RayQueryStructureT> func(sample, has_normals, surface);

        Ball3 nbd(sample->getPosition(), radius);
        sample_kdtree.processRangeUntil<IntersectionTester>(nbd, &func);

        if (sample->getNeighbors().size() >= min_degree)
          break;
        else
          radius *= RADIUS_EXPANSION_FACTOR;
      }
    }

    /**
     * Extract adjacencies between the original set of samples, by following geodesic paths in the graph of oversampled points.
     */
    void extractOriginalAdjacencies(TheaArray<SurfaceSample *> & sample_ptrs);

    /** Update the average separation to the correct value. We'll compute this by brute force instead of sampling for now. */
    void updateAverageSeparation()
    {
      avg_separation = 0;
      long num_edges = 0;  // double-counts, but we'll ignore that for now
      for (size_t i = 0; i < samples.size(); ++i)
      {
        SurfaceSample::NeighborSet const & nbrs = samples[i].getNeighbors();
        for (int j = 0; j < nbrs.size(); ++j)
          avg_separation += nbrs[j].getSeparation();

        num_edges += (long)nbrs.size();
      }

      if (num_edges > 0)
        avg_separation /= num_edges;
    }

    Options options;            ///< Options for creating the graph.
    bool has_normals;           ///< Do the samples have associated surface normals?
    SampleArray samples;        ///< Array of samples constituting the vertex set of the graph.
    SampleArray dense_samples;  ///< Array of samples constituting an oversampling of the surface.
    Real avg_separation;        ///< Average separation of samples.
    bool initialized;           ///< Has the graph been initialized?

}; // class SampleGraph

} // namespace Algorithms
} // namespace Thea

#endif
