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

#include "PointSet3.hpp"
#include "ShortestPaths.hpp"
#include <fstream>
#include <sstream>

namespace Thea {
namespace Algorithms {

PointSet3::PointSet3(PointSet3 const & src)
: graph(samples)
{
  *this = src;
}

namespace PointSet3Internal {

void
updateNeighborPointers(PointSet3::SampleArray & samples, PointSet3::SampleArray const & src_samples)
{
  alwaysAssertM(samples.size() == src_samples.size(),
                "PointSet3: Can't update sample neighbor pointers from source array of different size");

  for (size_t i = 0; i < samples.size(); ++i)
  {
    auto & nbrs = samples[i].getNeighbors();
    for (int j = 0; j < nbrs.size(); ++j)
    {
      size_t index = nbrs[j].getSample() - &src_samples[0];  // take advantage of array storage (this is NOT getIndex())
      debugAssertM(index >= 0 && index < samples.size(), "PointSet3: Can't get array index of neighboring sample");

      // Again, because of array storage, this should not break the relative ordering of neighbors with equal separation, since
      // pointer less-than is preserved
      const_cast<PointSet3::Sample::Neighbor &>(nbrs[j]).setSample(&samples[index]);
    }
  }
}

} // namespace PointSet3Internal

PointSet3 &
PointSet3::operator=(PointSet3 const & src)
{
  options = src.options;
  has_normals = src.has_normals;
  samples = src.samples;
  dense_samples = src.dense_samples;
  valid_graph = src.valid_graph;
  avg_separation = src.avg_separation;
  valid_scale = src.valid_scale;
  scale = src.scale;
  valid_kdtree = false;  // this must be recomputed

  // Update neighbor pointers
  PointSet3Internal::updateNeighborPointers(samples, src.samples);
  PointSet3Internal::updateNeighborPointers(dense_samples, src.dense_samples);

  return *this;
}

void
PointSet3::clear()
{
  has_normals = false;
  samples.clear();
  dense_samples.clear();
  valid_graph = true;
  avg_separation = 0;
  valid_scale = true;
  scale = 0;
  valid_kdtree = false;
}

namespace PointSet3Internal {

// A graph on samples specified as pointers, using the adjacency information already in the samples.
class SamplePointerGraph
{
  public:
    typedef PointSet3::Sample Sample;
    typedef Sample * VertexHandle;
    typedef Sample const * VertexConstHandle;
    typedef Sample ** VertexIterator;
    typedef Sample const * const * VertexConstIterator;
    typedef Sample::Neighbor * NeighborIterator;
    typedef Sample::Neighbor const * NeighborConstIterator;

    SamplePointerGraph(intx num_nodes_, Sample ** nodes_) : num_nodes(num_nodes_), nodes(nodes_) {}
    intx numVertices() const { return num_nodes; }

    VertexIterator verticesBegin() { return nodes; }
    VertexConstIterator verticesBegin() const { return nodes; }
    VertexIterator verticesEnd() { return nodes + num_nodes; }
    VertexConstIterator verticesEnd() const { return nodes + num_nodes; }
    VertexHandle getVertex(VertexIterator vi) { return *vi; }
    VertexConstHandle getVertex(VertexConstIterator vi) const { return *vi; }

    intx numNeighbors(VertexConstHandle vertex) const { return vertex->getNeighbors().size(); }
    NeighborIterator neighborsBegin(VertexHandle vertex)
    {
      return vertex->getNeighbors().isEmpty() ? nullptr : const_cast<Sample::Neighbor *>(&vertex->getNeighbors()[0]);
    }
    NeighborConstIterator neighborsBegin(VertexConstHandle vertex) const
    {
      return vertex->getNeighbors().isEmpty() ? nullptr : &vertex->getNeighbors()[0];
    }
    NeighborIterator neighborsEnd(VertexHandle vertex) { return neighborsBegin(vertex) + numNeighbors(vertex); }
    NeighborConstIterator neighborsEnd(VertexConstHandle vertex) const { return neighborsBegin(vertex) + numNeighbors(vertex); }
    VertexHandle getVertex(NeighborIterator ni) { return ni->getSample(); }
    VertexConstHandle getVertex(NeighborConstIterator ni) const { return ni->getSample(); }
    double distance(VertexConstHandle v, NeighborConstIterator ni) const { return ni->getSeparation(); }

  private:
    intx num_nodes;
    Sample ** nodes;

}; // class SamplePointerGraph

// Callback for shortest paths algorithm.
struct DijkstraCallback
{
  DijkstraCallback(PointSet3::Sample * sample_, intx num_orig_samples_, int max_nbrs_)
  : sample(sample_), num_orig_samples(num_orig_samples_), max_nbrs(max_nbrs_)
  {
    sample->getNeighbors().clear();
  }

  bool operator()(SamplePointerGraph::VertexHandle vertex, double distance, bool has_pred,
                  SamplePointerGraph::VertexHandle pred)
  {
    // Assume original samples are at the head of the list
    if (vertex->getIndex() != sample->getIndex() && vertex->getIndex() < num_orig_samples)
      sample->getNeighbors().insert(PointSet3::Sample::Neighbor(vertex, (Real)distance));

    return sample->getNeighbors().size() >= max_nbrs;
  }

  PointSet3::Sample * sample;
  intx num_orig_samples;
  int max_nbrs;
};

} // namespace PointSet3Internal

void
PointSet3::extractOriginalAdjacencies(intx num_samples, Sample ** sample_ptrs) const
{
  using namespace PointSet3Internal;

  SamplePointerGraph graph(num_samples, sample_ptrs);
  ShortestPaths<SamplePointerGraph> shortest_paths;

  // The graph is considered a mutable property of the point set so the const_casts below are ok because only neighbors are
  // modified
  SampleArray samples_with_new_nbrs(samples.size());
  for (size_t i = 0; i < samples.size(); ++i)
  {
    samples_with_new_nbrs[i] = samples[i];
    DijkstraCallback callback(&samples_with_new_nbrs[i], (intx)samples.size(), options.max_degree);
    shortest_paths.dijkstraWithCallback(graph, const_cast<Sample *>(&samples[i]), callback);
  }

  for (size_t i = 0; i < samples.size(); ++i)
    const_cast<Sample &>(samples[i]).getNeighbors() = samples_with_new_nbrs[i].getNeighbors();
}

bool
PointSet3::load(std::string const & samples_path, std::string const & graph_path)
{
  clear();

  // Load samples
  std::ifstream sin(samples_path.c_str());
  if (!sin)
  {
    THEA_ERROR << "PointSet3: Could not open samples file '" << samples_path << "' for reading";
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
      THEA_ERROR << "PointSet3: Some samples in '" << samples_path << "' have normals and some don't";
      return false;
    }

    samples.push_back(Sample((intx)samples.size(), 0));
    samples.back().setPosition(p);
    if (has_normals) samples.back().setNormal(n);
  }

  invalidateAll();  // graph, scale etc must be recomputed

  // Load graph, if a non-empty path was specified. Else, it must be computed with updateGraph().
  if (graph_path.empty())
    return true;

  std::ifstream gin(graph_path.c_str());
  if (!gin)
  {
    THEA_ERROR << "PointSet3: Could not open graph file '" << graph_path << "' for reading";
    return false;
  }

  if (std::getline(gin, line))
  {
    std::istringstream line_in(line);
    intx max_nbrs;
    if (!(line_in >> max_nbrs) || max_nbrs < 0)
    {
      THEA_ERROR << "PointSet3: Could not read valid maximum degree from '" << graph_path << '\'';
      return false;
    }

    options.setMaxDegree(max_nbrs);
  }
  else
  {
    THEA_ERROR << "PointSet3: Could not read maximum degree from '" << graph_path << '\'';
    return false;
  }

  Array<Sample::Neighbor> nbrs((size_t)options.max_degree);
  intx num_nbrs, nbr_index;
  for (size_t i = 0; i < samples.size(); ++i)
  {
    if (!std::getline(gin, line))
    {
      THEA_ERROR << "PointSet3: Could not read neighbors of vertex " << i << " from '" << graph_path << '\'';
      return false;
    }

    std::istringstream line_in(line);
    if (!(line_in >> num_nbrs) || num_nbrs < 0)
    {
      THEA_ERROR << "PointSet3: Could not read valid degree of vertex " << i << " from '" << graph_path << '\'';
      return false;
    }

    for (intx j = 0; j < num_nbrs; ++j)
    {
      if (!(line_in >> nbr_index) || nbr_index < 0 || nbr_index >= (intx)samples.size())
      {
        THEA_ERROR << "PointSet3: Could not read valid neighbor " << j << " of vertex " << i << " from '" << graph_path
                   << '\'';
        return false;
      }

      Sample * nbr_sample = &samples[(size_t)nbr_index];
      nbrs[(size_t)j] = Sample::Neighbor(nbr_sample, -1);
    }

    // Check if precomputed separations exist, else compute separations as Euclidean distances
    Real nbr_sep;
    if (line_in >> nbr_sep)
    {
      for (size_t j = 0; j < nbrs.size(); ++j)
      {
        if (j > 0)
        {
          if (!(line_in >> nbr_sep))
          {
            THEA_ERROR << "PointSet3: Could not read separation of neighbor " << j << " of vertex " << i << " from '"
                       << graph_path << '\'';
            return false;
          }
        }

        nbrs[j].setSeparation(nbr_sep);
      }
    }
    else
    {
      for (size_t j = 0; j < nbrs.size(); ++j)
      {
        nbr_sep = (samples[i].getPosition() - nbrs[j].getSample()->getPosition()).norm();
        nbrs[j].setSeparation(nbr_sep);
      }
    }

    samples[i].getNeighbors().setCapacity(options.max_degree);
    for (size_t j = 0; j < nbrs.size(); ++j)
      samples[i].getNeighbors().insert(nbrs[j]);
  }

  updateAverageSeparation();
  valid_graph = true;

  return true;
}

bool
PointSet3::save(std::string const & graph_path, std::string const & samples_path, bool write_distances) const
{
  // Save graph
  std::ofstream gout(graph_path.c_str(), std::ios::binary);
  if (!gout)
  {
    THEA_ERROR << "PointSet3: Could not open graph file '" << graph_path << "' for writing";
    return false;
  }

  gout << options.max_degree << '\n';

  for (size_t i = 0; i < samples.size(); ++i)
  {
    auto const & nbrs = samples[i].getNeighbors();
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
      THEA_ERROR << "PointSet3: Could not open samples file '" << samples_path << "' for writing";
      return false;
    }

    for (size_t i = 0; i < samples.size(); ++i)
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
