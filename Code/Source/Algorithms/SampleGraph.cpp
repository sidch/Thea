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

#include "SampleGraph.hpp"
#include "ShortestPaths.hpp"

namespace Thea {
namespace Algorithms {

SampleGraph::SampleGraph(SampleGraph const & src)
{
  *this = src;
}

namespace SampleGraphInternal {

void
updateNeighborPointers(SampleGraph::SampleArray & samples, SampleGraph::SampleArray const & src_samples)
{
  alwaysAssertM(samples.size() == src_samples.size(),
                "SampleGraph: Can't update sample neighbor pointers from source array of different size");

  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    SurfaceSample::NeighborSet & nbrs = samples[i].getNeighbors();
    for (int j = 0; j < nbrs.size(); ++j)
    {
      array_size_t index = nbrs[j].getSample() - &src_samples[0];  // take advantage of array storage (this is NOT getIndex())
      debugAssertM(index >= 0 && index < samples.size(), "SampleGraph: Can't get array index of neighboring sample");

      // Again, because of array storage, this should not break the relative ordering of neighbors with equal separation, since
      // pointer less-than is preserved
      const_cast<SurfaceSample::Neighbor &>(nbrs[j]).setSample(&samples[index]);
    }
  }
}

} // namespace SampleGraphInternal

SampleGraph &
SampleGraph::operator=(SampleGraph const & src)
{
  options = src.options;
  has_normals = src.has_normals;
  samples = src.samples;
  dense_samples = src.dense_samples;
  avg_separation = src.avg_separation;
  initialized = src.initialized;

  // Update neighbor pointers
  updateNeighborPointers(samples, src.samples);
  updateNeighborPointers(dense_samples, src.dense_samples);

  return *this;
}

namespace SampleGraphInternal {

// A graph on samples specified as pointers, using the adjacency information already in the samples.
class SamplePointerGraph
{
  public:
    typedef TheaArray<SurfaceSample *> NodeArray;
    typedef SurfaceSample * VertexHandle;
    typedef SurfaceSample const * VertexConstHandle;
    typedef NodeArray::iterator VertexIterator;
    typedef NodeArray::const_iterator VertexConstIterator;
    typedef SurfaceSample::Neighbor * NeighborIterator;
    typedef SurfaceSample::Neighbor const * NeighborConstIterator;

    SamplePointerGraph(NodeArray * nodes_) : nodes(nodes_) {}
    long numVertices() const { return (long)nodes->size(); }

    VertexIterator verticesBegin() { return nodes->begin(); }
    VertexConstIterator verticesBegin() const { return nodes->begin(); }
    VertexIterator verticesEnd() { return nodes->end(); }
    VertexConstIterator verticesEnd() const { return nodes->end(); }
    VertexHandle getVertex(VertexIterator vi) { return *vi; }
    VertexConstHandle getVertex(VertexConstIterator vi) const { return *vi; }

    long numNeighbors(VertexConstHandle vertex) const { return vertex->getNeighbors().size(); }
    NeighborIterator neighborsBegin(VertexHandle vertex)
    {
      return vertex->getNeighbors().isEmpty() ? NULL : const_cast<SurfaceSample::Neighbor *>(&vertex->getNeighbors()[0]);
    }
    NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
    {
      return vertex->getNeighbors().isEmpty() ? NULL : &vertex->getNeighbors()[0];
    }
    NeighborIterator neighborsEnd(VertexHandle vertex) { return neighborsBegin(vertex) + numNeighbors(vertex); }
    NeighborConstIterator neighborsEnd(VertexConstHandle vertex) const { return neighborsBegin(vertex) + numNeighbors(vertex); }
    VertexHandle getVertex(NeighborIterator ni) { return ni->getSample(); }
    VertexConstHandle getVertex(NeighborConstIterator ni) const { return ni->getSample(); }
    double distance(VertexConstHandle v, NeighborConstIterator ni) const { return ni->getSeparation(); }

  private:
    NodeArray * nodes;

}; // class SamplePointerGraph

// Callback for shortest paths algorithm.
struct DijkstraCallback
{
  DijkstraCallback(SurfaceSample * sample_, long num_orig_samples_, int max_nbrs_)
  : sample(sample_), num_orig_samples(num_orig_samples_), max_nbrs(max_nbrs_)
  {
    sample->getNeighbors().clear();
  }

  bool operator()(SamplePointerGraph::VertexHandle vertex, double distance, bool has_pred,
                  SamplePointerGraph::VertexHandle pred)
  {
    if (vertex->getIndex() != sample->getIndex() && vertex->getIndex() < num_orig_samples)
      sample->getNeighbors().insert(SurfaceSample::Neighbor(vertex, (Real)distance));

    return sample->getNeighbors().size() >= max_nbrs;
  }

  SurfaceSample * sample;
  long num_orig_samples;
  int max_nbrs;
};

} // namespace SampleGraphInternal

void
SampleGraph::extractOriginalAdjacencies(TheaArray<SurfaceSample *> & sample_ptrs)
{
  using namespace SampleGraphInternal;

  SamplePointerGraph graph(&sample_ptrs);
  ShortestPaths<SamplePointerGraph> shortest_paths;

  TheaArray<SurfaceSample> samples_with_new_nbrs(samples.size());
  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    samples_with_new_nbrs[i] = samples[i];
    DijkstraCallback callback(&samples_with_new_nbrs[i], (long)samples.size(), options.max_degree);
    shortest_paths.dijkstraWithCallback(graph, &samples[i], &callback);
  }

  for (array_size_t i = 0; i < samples.size(); ++i)
    samples[i].getNeighbors() = samples_with_new_nbrs[i].getNeighbors();
}

} // namespace Algorithms
} // namespace Thea
