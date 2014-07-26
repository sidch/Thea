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
#include "../BoundedSortedArray.hpp"
#include "../Vector3.hpp"
#include "PointTraitsN.hpp"

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

    /** Get the separation of this sample from the source sample. */
    Real getSeparation() const { return separation; }

    /** Compare neighbors by their separation from the source. */
    bool operator<(NeighboringSample const & rhs) const
    {
      return separation < rhs.separation || (separation == rhs.separation && std::less<SurfaceSample *>()(sample, rhs.sample));
    }

  private:
    friend class SampleGraph;

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
 * A class used to internal store surface samples with adjacency information. Encapsulates a sample position, normal and list of
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

  private:
    /** Get a non-const reference to the set of neighboring samples. */
    NeighborSet & getNeighbors() { return nbrs; }

    friend class SampleGraph;

    long index;        ///< Index of the sample.
    Vector3 p;         ///< Sample position.
    Vector3 n;         ///< Sample normal.
    NeighborSet nbrs;  ///< Set of neighboring samples.

}; // class SurfaceSample

} // namespace SampleGraphInternal

/** Proximity graph on surface samples. */
class SampleGraph
{
  private:
    /** Dummy surface representation used to disambiguate specialization of init(). */
    struct DummyRayQueryStructure3 : public RayQueryStructureN<3, Real>
    {
      typename RayQueryStructureN<3, Real> BaseT;

      template <typename RayIntersectionTesterT> bool rayIntersects(RayT const & ray, Real max_time = -1) const
      { return false; }

      template <typename RayIntersectionTesterT> Real rayIntersectionTime(RayT const & ray, Real max_time = -1) const
      { return -1; }

      template <typename RayIntersectionTesterT>
      BaseT::RayStructureIntersectionT rayStructureIntersection(RayT const & ray, Real max_time = -1) const
      { return BaseT::RayStructureIntersectionT; }

    }; // struct DummyRayQueryStructure3

  public:
    /** Options for constructing the sample graph. */
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

    //=========================================================================================================================
    //
    // Types and functions required to make this a bona fide graph satisfying the IsAdjacencyGraph concept.
    //
    //=========================================================================================================================

    typedef SurfaceSample *                     VertexHandle;           ///< Handle to a graph vertex (a sample).
    typedef SurfaceSample const *               VertexConstHandle;      ///< Const handle to a graph vertex (a sample).
    typedef SampleArray::iterator               VertexIterator;         ///< Iterator over vertices.
    typedef SampleArray::const_iterator         VertexConstIterator;    ///< Const iterator over vertices.
    typedef Neighbor *                          NeighborIterator;       ///< Iterator over neighbors of a vertex.
    typedef Neighbor const *                    NeighborConstIterator;  ///< Const iterator over neighbors of a vertex.

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
    long numNeighbors(VertexConstHandle vertex) const { return vertex->nbrs.size(); }

    /** Get an iterator to the first neighbor of a vertex. */
    NeighborIterator neighborsBegin(VertexHandle vertex)
    {
      return vertex->nbrs.isEmpty() ? NULL : const_cast<Neighbor *>(&vertex->nbrs[0]);
    }

    /** Get a const iterator to the first neighbor of a vertex. */
    NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
    {
      return vertex->nbrs.isEmpty() ? NULL : &vertex->nbrs[0];
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
    SampleGraph(Options const & options_ = Options::default())
    : options(options_), has_normals(false)
    {
      alwaysAssertM(options.max_degree > 0, "SampleGraph: Maximum degree must be positive");
    }

    /** Get the number of samples in the graph. */
    long numSamples() const { return (long)samples.size(); }

    /** Set the sample positions and (optionally) normals. All prior samples will be cleared. */
    void setPositions(TheaArray<Vector3> const & positions, TheaArray<Vector3> const * normals = NULL)
    {
      alwaysAssertM(!normals || normals->size() == positions.size(),
                    "SampleGraph: Number of normals does not match number of samples");

      samples.resize(positions.size());
      for (array_size_t i = 0; i < samples.size(); ++i)
      {
        samples[i].index = (long)i;
        samples[i].p = positions[i];
        samples[i].nbrs.setCapacity(options.max_degree);
      }

      if (normals)
      {
        for (array_size_t i = 0; i < samples.size(); ++i)
          samples[i].n = (*normals)[i];

        has_normals = true;
      }
      else
        has_normals = false;
    }

    /**
     * Specify a previously computed oversampling of the surface. These additional samples will not constitute vertices of the
     * final graph, but will be used to more accurately compute adjacencies. If normals are not specified, they must also not
     * have been specified for the main samples.
     */
    void setOversampling(TheaArray<Vector3> const & dense_positions, TheaArray<Vector3> const * dense_normals = NULL)
    {
      alwaysAssertM(!has_normals || dense_normals != NULL,
                    "SampleGraph: Main samples have normals, so oversampling must also have normals");
      alwaysAssertM(!dense_normals || dense_normals->size() == dense_positions.size(),
                    "SampleGraph: Number of normals does not match number of samples");

      dense_samples.resize(dense_positions.size());
      for (array_size_t i = 0; i < dense_samples.size(); ++i)
      {
        dense_samples[i].p = dense_positions[i];
        dense_samples[i].nbrs.setCapacity(options.max_degree);
      }

      if (dense_normals)
      {
        for (array_size_t i = 0; i < dense_samples.size(); ++i)
          dense_samples[i].n = (*dense_normals)[i];
      }
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
      TheaArray<SurfaceSample *> sample_ptrs(samples.size() + dense_samples.size());

      for (array_size_t i = 0; i < samples.size(); ++i)
        sample_ptrs[i] = &samples[i];

      for (array_size_t i = 0; i < dense_samples.size(); ++i)
      {
        dense_samples[i].index = (long)(samples.size() + i);
        sample_ptrs[samples.size() + i] = &dense_samples[i];
      }

      typedef KDTreeN<3, SurfaceSample *> SampleKDTree;
      SampleKDTree sample_kdtree(sample_ptrs.begin(), sample_ptrs.end());
    }

    Options options;      ///< Options for creating the graph.
    bool has_normals;     ///< Do the samples have associated surface normals?
    SampleArray samples;  ///< Array of samples constituting the vertex set of the graph.

}; // class SampleGraph

} // namespace Algorithms
} // namespace Thea

#endif
