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
#include <fstream>
#include <sstream>

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

void
SampleGraph::clear()
{
  has_normals = false;
  samples.clear();
  dense_samples.clear();
  avg_separation = 0;
  initialized = false;
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

bool
SampleGraph::load(std::string const & graph_path, std::string const & samples_path)
{
  clear();

  // Load samples
  std::ifstream sin(samples_path.c_str());
  if (!sin)
  {
    THEA_ERROR << "SampleGraph: Could not open samples file '" << samples_path << "' for reading";
    return false;
  }

  has_normals = false;

  std::string line;
  Vector3 p, n;
  while (std::getline(sin, line))
  {
    std::istringstream line_in(line);
    if (!(line_in >> p[0] >> p[1] >> p[2]))
    {
      THEA_ERROR << "Could not read sample " << samples.size() << " from '" << samples_path << '\'';
      return false;
    }

    bool sample_has_normal = (bool)(line_in >> n[0] >> n[1] >> n[2]);
    if (samples.empty())
      has_normals = sample_has_normal;
    else if (has_normals != sample_has_normal)
    {
      THEA_ERROR << "SampleGraph: Some samples in '" << samples_path << "' have normals and some don't";
      return false;
    }

    samples.push_back(SurfaceSample((long)samples.size(), 0));
    samples.back().setPosition(p);
    if (has_normals) samples.back().setNormal(n);
  }

  // Load graph
  std::ifstream gin(graph_path.c_str());
  if (!gin)
  {
    THEA_ERROR << "SampleGraph: Could not open graph file '" << graph_path << "' for reading";
    return false;
  }

  if (std::getline(gin, line))
  {
    std::istringstream line_in(line);
    long max_nbrs;
    if (!(line_in >> max_nbrs) || max_nbrs < 0)
    {
      THEA_ERROR << "SampleGraph: Could not read valid maximum degree from '" << graph_path << '\'';
      return false;
    }

    options.setMaxDegree(max_nbrs);
  }
  else
  {
    THEA_ERROR << "SampleGraph: Could not read maximum degree from '" << graph_path << '\'';
    return false;
  }

  TheaArray<SurfaceSample::Neighbor> nbrs((array_size_t)options.max_degree);
  long num_nbrs, nbr_index;
  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    if (!std::getline(gin, line))
    {
      THEA_ERROR << "SampleGraph: Could not read neighbors of vertex " << i << " from '" << graph_path << '\'';
      return false;
    }

    std::istringstream line_in(line);
    if (!(line_in >> num_nbrs) || num_nbrs < 0)
    {
      THEA_ERROR << "SampleGraph: Could not read valid degree of vertex " << i << " from '" << graph_path << '\'';
      return false;
    }

    for (long j = 0; j < num_nbrs; ++j)
    {
      if (!(line_in >> nbr_index) || nbr_index < 0 || nbr_index >= (long)samples.size())
      {
        THEA_ERROR << "SampleGraph: Could not read valid neighbor " << j << " of vertex " << i << " from '" << graph_path
                   << '\'';
        return false;
      }

      SurfaceSample * nbr_sample = &samples[(array_size_t)nbr_index];
      nbrs[(array_size_t)j] = SurfaceSample::Neighbor(nbr_sample, -1);
    }

    // Check if precomputed separations exist, else compute separations as Euclidean distances
    Real nbr_sep;
    if (line_in >> nbr_sep)
    {
      for (array_size_t j = 0; j < nbrs.size(); ++j)
      {
        if (j > 0)
        {
          if (!(line_in >> nbr_sep))
          {
            THEA_ERROR << "SampleGraph: Could not read separation of neighbor " << j << " of vertex " << i << " from '"
                       << graph_path << '\'';
            return false;
          }
        }

        nbrs[j].setSeparation(nbr_sep);
      }
    }
    else
    {
      for (array_size_t j = 0; j < nbrs.size(); ++j)
      {
        nbr_sep = (samples[i].getPosition() - nbrs[j].getSample()->getPosition()).length();
        nbrs[j].setSeparation(nbr_sep);
      }
    }

    samples[i].getNeighbors().setCapacity(options.max_degree);
    for (array_size_t j = 0; j < nbrs.size(); ++j)
      samples[i].getNeighbors().insert(nbrs[j]);
  }

  updateAverageSeparation();
  initialized = true;

  return true;
}

bool
SampleGraph::save(std::string const & graph_path, std::string const & samples_path, bool write_distances) const
{
  // Save graph
  std::ofstream gout(graph_path.c_str(), std::ios::binary);
  if (!gout)
  {
    THEA_ERROR << "SampleGraph: Could not open graph file '" << graph_path << "' for writing";
    return false;
  }

  gout << options.max_degree << '\n';

  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    SurfaceSample::NeighborSet const & nbrs = samples[i].getNeighbors();
    gout << nbrs.size();

    for (int j = 0; j < nbrs.size(); ++j)
      gout << ' ' << nbrs[j].getSample()->getIndex();

    if (write_distances)
    {
      for (int j = 0; j < nbrs.size(); ++j)
        gout << ' ' << nbrs[j].getSeparation();
    }

    gout << '\n';
  }

  // Save samples
  if (!samples_path.empty())
  {
    std::ofstream sout(samples_path.c_str(), std::ios::binary);
    if (!sout)
    {
      THEA_ERROR << "SampleGraph: Could not open samples file '" << samples_path << "' for writing";
      return false;
    }

    for (array_size_t i = 0; i < samples.size(); ++i)
    {
      Vector3 const & p = samples[i].getPosition();
      sout << p[0] << ' ' << p[1] << ' ' << p[2];
      if (has_normals)
      {
        Vector3 const & n = samples[i].getNormal();
        sout << ' ' << n[0] << ' ' << n[1] << ' ' << n[2];
      }

      sout << '\n';
    }
  }

  return true;
}

} // namespace Algorithms
} // namespace Thea
