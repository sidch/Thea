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

#ifndef __Thea_Algorithms_PointSet3_hpp__
#define __Thea_Algorithms_PointSet3_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "BestFitSphere3.hpp"
#include "BvhN.hpp"
#include "Filter.hpp"
#include "Iterators.hpp"
#include "IntersectionTester.hpp"
#include "MeshSampler.hpp"
#include "MetricL2.hpp"
#include "NormalTraitsN.hpp"
#include "PointCollectorN.hpp"
#include "PointTraitsN.hpp"
#include "RayIntersectionTester.hpp"
#include "RayQueryStructureN.hpp"
#include "SamplePoint3.hpp"
#include "../MatVec.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

/**
 * A set of 3D points, with a proximity graph and BVH on the samples. Typical usage is to first specify the sample positions
 * (and optionally normals) via setSamples(). Then, specify an optional dense oversampling via setOversampling() -- the
 * oversampling makes it easier to verify if two samples are actually adjacent on the underlying manifold (if there is one) when
 * computing its proximity graph. This graph is lazily computed by getGraph(): it can be improved if a representation of the
 * underlying surface is provided via an explicit call to <tt>updateGraph(RayQueryStructureT const *)</tt>.
 */
class PointSet3
{
  public:
    /** %Options for constructing and processing the point set. */
    class Options
    {
      public:
        /** Set the maximum degree of the proximity graph. */
        Options & setMaxDegree(int max_degree_) { max_degree = max_degree_; return *this; }

        /** Construct with default values. */
        Options() : max_degree(8) {}

        /** Get a set of options with default values. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        int max_degree;  ///< Maximum degree of the proximity graph.

        friend class PointSet3;

    }; // class Options

    typedef SamplePoint3 Sample;          ///< A single sample point.
    typedef Array<Sample> SampleArray;    ///< Array of samples.
    typedef BvhN<Sample *, 3> SampleBvh;  ///< BVH on surface samples.

    /** Proximity graph on surface samples, satisfying the IsAdjacencyGraph concept. */
    class SampleGraph
    {
      public:
        typedef Sample *                     VertexHandle;           ///< Handle to a graph vertex (a sample).
        typedef Sample const *               VertexConstHandle;      ///< Const handle to a graph vertex (a sample).
        typedef SampleArray::iterator        VertexIterator;         ///< Iterator over vertices.
        typedef SampleArray::const_iterator  VertexConstIterator;    ///< Const iterator over vertices.
        typedef Sample::Neighbor *           NeighborIterator;       ///< Iterator over neighbors of a vertex.
        typedef Sample::Neighbor const *     NeighborConstIterator;  ///< Const iterator over neighbors of a vertex.

        /** Get the number of vertices (samples) in the graph. */
        intx numVertices() const { return (intx)samples.size(); }

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
        intx numNeighbors(VertexConstHandle vertex) const { return vertex->numNeighbors(); }

        /** Get an iterator to the first neighbor of a vertex. */
        NeighborIterator neighborsBegin(VertexHandle vertex)
        {
          return vertex->numNeighbors() <= 0 ? nullptr : const_cast<Sample::Neighbor *>(&vertex->getNeighbor(0));
        }

        /** Get a const iterator to the first neighbor of a vertex. */
        NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
        {
          return vertex->numNeighbors() <= 0 ? nullptr : &vertex->getNeighbor(0);
        }

        /** Get an iterator to one position beyond the last neighbor of a vertex. */
        NeighborIterator neighborsEnd(VertexHandle vertex) { return neighborsBegin(vertex) + numNeighbors(vertex); }

        /** Get a const iterator to one position beyond the last neighbor of a vertex. */
        NeighborConstIterator neighborsEnd(VertexConstHandle vertex) const
        { return neighborsBegin(vertex) + numNeighbors(vertex); }

        /** Get a handle to the neighboring vertex referenced by an iterator. */
        VertexHandle getVertex(NeighborIterator ni) { return ni->getSample(); }

        /** Get a handle to the neighboring vertex referenced by a const iterator. */
        VertexConstHandle getVertex(NeighborConstIterator ni) const { return ni->getSample(); }

        /** Get the distance between a vertex and its neighbor. */
        double distance(VertexConstHandle v, NeighborConstIterator ni) const { return ni->getSeparation(); }

      private:
        /** Construct the graph to reference an external set of samples. */
        SampleGraph(SampleArray & samples_) : samples(samples_) {}

        SampleArray & samples;  ///< The referenced set of samples which are the vertices of the graph.

        friend class PointSet3;

    }; // class SampleGraph

    //=========================================================================================================================
    //
    // Functions to create the graph from a set of samples and an optional base shape.
    //
    //=========================================================================================================================

    /** Construct with a set of options for creating the graph. */
    PointSet3(Options const & options_ = Options::defaults())
    : options(options_), has_normals(false), valid_graph(true), graph(samples), avg_separation(0), valid_scale(true), scale(0),
      valid_bvh(true)
    {
      alwaysAssertM(options.max_degree >= 0, "PointSet3: Maximum degree must be non-negative");
    }

    /** Copy constructor. */
    PointSet3(PointSet3 const & src);

    /** Assignment operator. */
    PointSet3 & operator=(PointSet3 const & src);

    /** Check if the samples have normals. The return value is arbitrary if there are no samples (or oversamples) at all. */
    bool hasNormals() const { return has_normals; }

    /** Get the number of samples. */
    intx numSamples() const { return (intx)samples.size(); }

    /** Get the set of samples. */
    SampleArray const & getSamples() const { return samples; }

    /** Get a sample by its index. */
    Sample const & getSample(intx index) const
    {
      theaAssertM(index >= 0 && index < (intx)samples.size(), "PointSet3: Sample index out of bounds");
      return samples[(size_t)index];
    }

    /**
     * Add new sample positions and (optionally) normals. If the existing samples have normals, the new ones must have normals
     * too.
     */
    void addSamples(intx num_samples, Vector3 const * positions, Vector3 const * normals = nullptr)
    {
      addSamples(num_samples, positions, normals, samples);
    }

    /**
     * Add new samples from an input sequence of objects that satisfy IsPointN and (optionally) HasNormalN. The traits classes
     * PointTraitsN and NormalTraitsN will be used to extract the position and normal respectively. If the existing samples have
     * normals, the new ones must have normals too.
     */
    template <typename SampleIteratorT>
    void addSamples(SampleIteratorT samples_begin, SampleIteratorT samples_end)
    {
      addSamples(samples_begin, samples_end, samples);
    }

    /**
     * Adds points sampled from a surface mesh. The new samples will have normals unless any previously added samples lack them.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape. If negative, a default number of samples are generated.
     * @param smooth If true, sample normals will be computed by averaging vertex normals rather than by using piecewise
     *   constant triangle normals.
     */
    template <typename MeshT> void addSamples(MeshT const & mesh, intx num_samples = -1, bool smooth = false)
    {
      addSamples(MeshSampler<MeshT>(mesh), num_samples, smooth, samples);
    }

    /**
     * Adds points sampled from a surface mesh group. The new samples will have normals unless any previously added samples lack
     * them.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape. If negative, a default number of samples are generated.
     * @param smooth If true, sample normals will be computed by averaging vertex normals rather than by using piecewise
     *   constant triangle normals.
     */
    template <typename MeshT>
    void addSamples(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, bool smooth = false)
    {
      addSamples(MeshSampler<MeshT>(mesh_group), num_samples, smooth, samples);
    }

    /**
     * Add points to an oversampling of the surface. These additional samples will not constitute vertices of the final graph,
     * but will be used to more accurately compute adjacencies. If the main samples have normals, then the additional samples
     * must also have normals specified.
     */
    void addOversampling(intx num_samples, Vector3 const * dense_positions, Vector3 const * dense_normals = nullptr)
    {
      addSamples(num_samples, dense_positions, dense_normals, dense_samples);
    }

    /** Get the number of oversamples (points in the dense oversampling). */
    intx numOversamples() const { return (intx)dense_samples.size(); }

    /**
     * Add points to an oversampling of the surface. These additional samples will not constitute vertices of the final graph,
     * but will be used to more accurately compute adjacencies. The input sequence should consist of objects that satisfy
     * IsPointN and (optionally) HasNormalN. The traits classes PointTraitsN and NormalTraitsN will be used to extract the
     * position and normal respectively. Normals of the oversampled points will be used if and only if the main samples have
     * normals (or if there are no main samples, in which case normals of the oversampling will be saved in case we need them
     * later).
     */
    template <typename SampleIteratorT>
    void addOversampling(SampleIteratorT samples_begin, SampleIteratorT samples_end)
    {
      addSamples(samples_begin, samples_end,
                 (samples.empty() || has_normals),  // if no main samples, cache dense normals anyway in case we need them
                 dense_samples);
    }

    /**
     * Add points sampled from a surface mesh to an oversampling of the surface. The new samples will have normals if and only
     * if the main samples have normals.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of additional samples to compute on the shape. If negative, a default number of samples are
     *   generated.
     * @param smooth If true, sample normals will be computed by averaging vertex normals rather than by using piecewise
     *   constant triangle normals.
     */
    template <typename MeshT> void addOversampling(MeshT const & mesh, intx num_samples = -1, bool smooth = false)
    {
      static intx const DEFAULT_NUM_DENSE_SAMPLES = 10000;
      if (num_samples < 0) num_samples = DEFAULT_NUM_DENSE_SAMPLES;

      addSamples(MeshSampler<MeshT>(mesh), num_samples, smooth, dense_samples);
    }

    /**
     * Add points sampled from a surface mesh group to an oversampling of the surface. The new samples will have normals if and
     * only if the main samples have normals.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of additional samples to compute on the shape. If negative, a default number of samples are
     *   generated.
     * @param smooth If true, sample normals will be computed by averaging vertex normals rather than by using piecewise
     *   constant triangle normals.
     */
    template <typename MeshT>
    void addOversampling(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, bool smooth = false)
    {
      static intx const DEFAULT_NUM_DENSE_SAMPLES = 10000;
      if (num_samples < 0) num_samples = DEFAULT_NUM_DENSE_SAMPLES;

      addSamples(MeshSampler<MeshT>(mesh_group), num_samples, smooth, dense_samples);
    }

    /**
     * Get the proximity graph on surface samples. The graph will be lazily computed by this function if it is currently
     * invalid.
     */
    SampleGraph const & getGraph() const { updateGraph(); return graph; }

    /** Invalidate the current proximity graph. Vertex information will be retained but edges will be considered invalid. */
    void invalidateGraph() { valid_graph = false; }

    /**
     * Construct the sample proximity graph from the points, without tests that require access to the underlying surface. This
     * function is automatically called by getGraph() to update the graph if it is invalid. It can also be called explicitly to
     * ensure the proximity graph exists.
     *
     * @see updateGraph(RayQueryStructureT const *)
     */
    void updateGraph() const
    {
      updateGraph<DummyRayQueryStructure3>(nullptr);
    }

    /**
     * Construct the sample proximity graph with the aid of a known representation of the underlying surface (such as a BVH
     * on mesh triangles). Additional adjacency tests, including those that test if the two samples are on opposite sides of a
     * thin surface, will be performed by this explicit graph update.
     *
     * @see updateGraph()
     */
    template <typename RayQueryStructureT> void updateGraph(RayQueryStructureT const * surface) const
    {
      if (valid_graph)  // don't initialize twice
        return;

      SampleBvh * h = nullptr;
      SampleBvh agg_bvh;
      if (!dense_samples.empty())
      {
        // Aggregate samples
        Array<Sample *> sample_ptrs(samples.size() + dense_samples.size());

        for (size_t i = 0; i < samples.size(); ++i)
          sample_ptrs[i] = const_cast<Sample *>(&samples[i]);

        for (size_t i = 0; i < dense_samples.size(); ++i)
        {
          auto & ds = const_cast<Sample &>(dense_samples[i]);
          ds.setIndex((intx)(samples.size() + i));  // indices of dense samples don't matter and can be changed internally
          sample_ptrs[samples.size() + i] = &ds;
        }

        // Clear any prior adjacency data
        for (size_t i = 0; i < sample_ptrs.size(); ++i)
          sample_ptrs[i]->clearNeighbors();

        if (sample_ptrs.size() <= 1)  // nothing to do
          return;

        // Create BVH on aggregated samples
        agg_bvh.init(sample_ptrs.begin(), sample_ptrs.end());
        h = &agg_bvh;
      }
      else
        h = const_cast<SampleBvh *>(&getBvh());

      // Get a measure of the average pairwise separation of samples
      auto h_num_elems = h->numElements();
      auto h_elems = h->getElements();
      avg_separation = 0;
      intx num_trials = std::min(100L, h_num_elems);
      for (intx i = 0; i < num_trials; ++i)
      {
        auto index = Random::common().integer(0, (int32)h_num_elems - 1);
        FilterSelf filter(h_elems[index]);
        h->pushFilter(&filter);
          intx nn_index = h->closestElement<MetricL2>(h_elems[index]->getPosition());
          alwaysAssertM(nn_index >= 0, "PointSet3: Nearest neighbor of sample not found");
          avg_separation += (h_elems[nn_index]->getPosition() - h_elems[index]->getPosition()).squaredNorm();
        h->popFilter();
      }

      avg_separation = std::sqrt(avg_separation / num_trials);  // RMS average

      // Find the neighbors of each sample
      Real sep_scale = std::sqrt((Real)options.max_degree);  // assume uniform distribution on 2D surface
      for (intx i = 0; i < h_num_elems; ++i)
        findSampleNeighbors(h_elems[i], *h, sep_scale * avg_separation, surface);

      // Extract adjacencies between original set of samples
      if (!dense_samples.empty())
        extractOriginalAdjacencies(h_num_elems, h_elems);

      // Update the average separation to the correct value
      updateAverageSeparation();

      valid_graph = true;
    }

    /**
     * Get the average separation of samples from their neighbors, that is, the average edge length of the graph. This may be
     * computed by sampling and hence should not be assumed to be a true average, just a representative value. This function,
     * like getGraph(), automatically computes the graph if it is invalid.
     */
    Real getAverageSeparation() const { updateGraph(); return avg_separation; }

    /**
     * Get the scale of the shape, which can be used to normalize distances etc. Currently, the shape is computed as the
     * diameter of the best-fit bounding sphere of the samples, unless a value has been manually specified via setScale().
     */
    Real getScale() const { updateScale(); return scale; }

    /**
     * Explicitly set the scale of the shape, which can be used to normalize distances etc. If a negative value is supplied, the
     * scale is reset to the default value computed by getScale(). Note that functions that change the surface, such as
     * addSamples(), will also cause the scale to be recomputed with the default calculation and you will need to manually
     * specify the scale again as needed.
     */
    void setScale(Real s)
    {
      scale = s;
      valid_scale = (s >= 0);
    }

    /** Get the BVH on the sample points. */
    SampleBvh const & getBvh() const { updateBvh(); return bvh; }

    /** Clear the graph. */
    void clear();

    /**
     * Load the samples and proximity graph from disk files. If \a graph_path is an empty string, the graph is not loaded, and
     * must be computed with an explicit call to updateGraph().
     */
    bool load(std::string const & samples_path, std::string const & graph_path = "");

    /**
     * Save the samples and proximity graph to disk. Either path may be an empty string, in which case the corresponding data
     * are not written.
     *
     * @param samples_path Output file for samples.
     * @param graph_path Output file for graph.
     * @param write_distances If true, the distance of each neighbor (which may be different from the Euclidean distance if
     *   the graph was computed via an oversampling) is also written to the file.
     */
    bool save(std::string const & samples_path, std::string const & graph_path = "", bool write_distances = false) const;

  private:
    /** Dummy surface representation used to disambiguate specialization of updateGraph(). */
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

    /** Invalidate all lazily computed/on-demand data. */
    void invalidateAll() { valid_graph = valid_scale = valid_bvh = false; }

    /**
     * Check if there is at least one sample with a normal vector in either the main samples array or in an additional
     * user-specified array (which may be the same as the main samples array).
     */
    bool hasAtLeastOneSampleNormal(SampleArray const & extra_arr) const
    { return has_normals && !(samples.empty() && extra_arr.empty()); }

    /**
     * If arr is the main samples array and it is empty, then set the has_normals flag to true if the new samples have normals
     * AND all existing dense samples have normals, else set has_normals to false. Else, if arr is the dense array, update
     * has_normals if and only if both main and dense arrays are empty.
     */
    void updateHasNormals(SampleArray const & arr, bool new_samples_have_normals)
    {
      if (samples.empty() && arr.empty())
        has_normals = new_samples_have_normals && (dense_samples.empty() || has_normals);
    }

    /** Helper function for adding new sample positions and (optionally) normals. */
    void addSamples(intx num_samples, Vector3 const * positions, Vector3 const * normals, SampleArray & arr)
    {
      alwaysAssertM(num_samples >= 0, "PointSet3: Cannot specify a negative number of samples");
      alwaysAssertM(num_samples == 0 || positions, "PointSet3: Sample positions must be specified");
      alwaysAssertM(!(hasAtLeastOneSampleNormal(arr) && !normals),
                    "PointSet3: Existing samples have normals, so new ones must have normals too");

      updateHasNormals(arr, (bool)normals);

      for (intx i = 0; i < num_samples; ++i)
      {
        intx index = (intx)arr.size();
        arr.push_back(Sample(index, options.max_degree));
        arr.back().setPosition(positions[i]);
        if (has_normals) { arr.back().setNormal(normals[i]); }
      }

      invalidateAll();
    }

    /**
     * Helper function for adding new samples from an input sequence of objects that satisfy IsPointN and (optionally)
     * HasNormalN. The traits classes PointTraitsN and NormalTraitsN will be used to extract the position and normal
     * respectively.
     */
    template <typename SampleIteratorT>
    void addSamples(SampleIteratorT samples_begin, SampleIteratorT samples_end, SampleArray & arr)
    {
      typedef typename std::iterator_traits<SampleIteratorT>::value_type SampleT;

      alwaysAssertM(!(hasAtLeastOneSampleNormal(arr) && !HasNormalN<SampleT, 3>::value),
                    "PointSet3: Existing samples have normals, so new ones must have normals too");

      updateHasNormals(arr, HasNormalN<SampleT, 3>::value);

      for (auto si = samples_begin; si != samples_end; ++si)
      {
        intx index = (intx)arr.size();
        arr.push_back(Sample(index, options.max_degree));
        arr.back().setPosition(PointTraitsN<SampleT, 3>::getPosition(*si));
        if (has_normals) { arr.back().setNormal(NormalTraitsN<SampleT, 3>::getNormal(*si)); }
      }

      invalidateAll();
    }

    /** Get the smoothed normal at a point on a triangle with vertex normals. */
    template <typename TriangleT> Vector3 smoothNormal(TriangleT const & tri, Vector3 const & p)
    {
      Vector3 b = tri.barycentricCoordinates(p);

      Vector3 n0 = tri.getVertices().getVertexNormal(0);
      Vector3 n1 = tri.getVertices().getVertexNormal(1);
      Vector3 n2 = tri.getVertices().getVertexNormal(2);

      return b[0] * n0 + b[1] * n1 + b[2] * n2;
    }

    /** Helper function for computing and adding sample points with normals on a mesh. */
    template <typename MeshT>
    void addSamples(MeshSampler<MeshT> const & sampler, intx num_samples, bool smooth, SampleArray & arr)
    {
      static intx const DEFAULT_NUM_SAMPLES = 5000;
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      updateHasNormals(arr, true);

      Array<Vector3> positions;
      Array< typename MeshSampler<MeshT>::Triangle const * > tris;
      sampler.sampleEvenlyByArea(num_samples, positions, nullptr, &tris);

      for (size_t i = 0; i < positions.size(); ++i)
      {
        intx index = (intx)arr.size();
        arr.push_back(Sample(index, options.max_degree));
        arr.back().setPosition(positions[i]);
        if (has_normals) { arr.back().setNormal(smooth ? smoothNormal(*tris[i], positions[i]) : tris[i]->getNormal()); }
      }

      invalidateAll();
    }

    /** Compute the scale of the shape as the diameter of the bounding sphere of the samples. */
    void updateScale() const
    {
      if (valid_scale) return;

      BestFitSphere3 bsphere;
      PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(samples.begin(), samples.end());
      scale = bsphere.getDiameter();

      valid_scale = true;
    }

    /** Allows every sample except one. */
    struct FilterSelf : public Filter<Sample *>
    {
      FilterSelf(Sample * self_) : self(self_) {}
      bool allows(Sample * const & s) const { return s != self; }

      Sample * self;

    }; // struct FilterSelf

    /** Functor that validates and adds neighbors to a given sample. */
    template <typename RayQueryStructureT> class NeighborFunctor
    {
      public:
        /** Constructor. */
        NeighborFunctor(Sample * sample_, bool has_normals_, RayQueryStructureT const * surface_)
        : sample(sample_), has_normals(has_normals_), surface(surface_)
        {
          theaAssertM(sample, "PointSet3: Can't create neighbor functor without valid source sample");
        }

        /** Called for every candidate neighbor. */
        bool operator()(intx index, Sample * nbr)
        {
          Real sep = 0;
          if (isValidNeighbor(nbr, sep))
            sample->addNeighbor(Sample::Neighbor(nbr, sep));

          return false;
        }

      private:
        /**
         * Check if a candidate neighboring sample is a valid neighbor, that is, if it is near the original sample as measured
         * on the surface. If the function returns true, it puts the computed separation between the original and candidate
         * samples in \a sep.
         */
        bool isValidNeighbor(Sample const * nbr, Real & sep)
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
            if (nbr == sample->getNeighbor(i).getSample())
              return false;

          // Compute the separation of the samples
          Vector3 diff = nbr->getPosition() - sample->getPosition();
          Real sqsep = diff.squaredNorm();
          sep = Math::fastSqrt(sqsep);

          // Reachability test: lift points off the surface by an amount proportional to their separation, and see if the line
          // connecting them is blocked by the surface
          if (surface && has_normals)
          {
            static Real const LIFT_FACTOR = 5;
            Vector3 lift_dir = (sample->getNormal() + nbr->getNormal()).normalized();
            typename RayQueryStructureT::RayT ray(sample->getPosition() + LIFT_FACTOR * sep * lift_dir, diff);
            if (surface->template rayIntersects<RayIntersectionTester>(ray, 1))
              return false;
          }

          return true;
        }

        Sample * sample;
        bool has_normals;
        RayQueryStructureT const * surface;

    }; // struct NeighborFunctor

    /** Find the samples adjacent to a given sample. */
    template <typename RayQueryStructureT>
    void findSampleNeighbors(Sample * sample, SampleBvh & h, Real init_radius, RayQueryStructureT const * surface) const
    {
      static int const MAX_ITERS = 3;
      static Real const RADIUS_EXPANSION_FACTOR = 2.0f;
      int min_degree = Math::clamp((int)(0.25 * options.max_degree), 4, options.max_degree);

      Real radius = init_radius;
      for (int i = 0; i < MAX_ITERS; ++i)
      {
        sample->clearNeighbors();
        NeighborFunctor<RayQueryStructureT> func(sample, has_normals, surface);

        Ball3 nbd(sample->getPosition(), radius);
        h.processRange<IntersectionTester>(nbd, func);

        if (sample->numNeighbors() >= min_degree)
          break;
        else
          radius *= RADIUS_EXPANSION_FACTOR;
      }
    }

    /**
     * Extract adjacencies between the original set of samples, by following geodesic paths in the graph of oversampled points.
     */
    void extractOriginalAdjacencies(intx num_samples, Sample ** sample_ptrs) const;

    /** Update the average separation to the correct value. We'll compute this by brute force instead of sampling for now. */
    void updateAverageSeparation() const
    {
      avg_separation = 0;
      intx num_edges = 0;  // double-counts, but we'll ignore that for now
      for (size_t i = 0; i < samples.size(); ++i)
      {
        for (auto const & n : samples[i].getNeighbors())
          avg_separation += n.getSeparation();

        num_edges += samples[i].numNeighbors();
      }

      if (num_edges > 0)
        avg_separation /= num_edges;
    }

    /**
     * Get a non-const reference to the internally created BVH, or null if no such BVH exists. <b>Use with caution</b>
     * -- you should never actually change the tree even if you have a non-const handle to it!
     */
    void updateBvh() const
    {
      if (valid_bvh) return;

      if (samples.empty())
        bvh.clear();
      else
      {
        // Need to use raw pointers else makePtrIterator can't create a read-write iterator
        bvh.init(makePtrIterator(const_cast<Sample *>(samples.data())),
                 makePtrIterator(const_cast<Sample *>(samples.data() + samples.size())));
      }

      valid_bvh = true;
    }

    Options options;              ///< Options for creating the graph.
    bool has_normals;             ///< Do the samples have associated surface normals?
    SampleArray samples;          ///< Array of samples constituting the vertex set of the graph.
    SampleArray dense_samples;    ///< Array of samples constituting an oversampling of the surface.

    mutable bool valid_graph;     ///< Has the proximity graph been initialized?
    mutable SampleGraph graph;    ///< Proximity graph on surface samples.
    mutable Real avg_separation;  ///< Average separation between neighboring samples.
    mutable bool valid_scale;     ///< Has the normalization length been computed?
    mutable Real scale;           ///< The normalization length.
    mutable bool valid_bvh;       ///< Has the sample BVH been initialized?
    mutable SampleBvh bvh;        ///< BVH on surface samples.

}; // class PointSet3

} // namespace Algorithms
} // namespace Thea

#endif
