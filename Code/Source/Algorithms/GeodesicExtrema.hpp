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
// First version: 2012
//
//============================================================================

#ifndef __Thea_Algorithms_GeodesicExtrema_hpp__
#define __Thea_Algorithms_GeodesicExtrema_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../GraphType.hpp"
#include "ShortestPaths.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Compute vertices that represent extrema of a graph, in the sense that they are, on average, further away from other vertices
 * in the graph than their neighbors. Formlly, extreme vertices are defined as local maxima of the function measuring average
 * geodesic distance of a vertex to all other vertices.
 *
 * GraphT must satisfy IsAdjacencyGraph and define a meaningful distance between adjacent vertices.
 */
template <typename GraphT>
class /* THEA_API */ GeodesicExtrema
{
  public:
    typedef GraphT Graph;  ///< The graph type.
    typedef typename GraphT::Vertex Vertex;  ///< A vertex of the graph.
    typedef typename GraphT::Vertex const * VertexConstHandle;  ///< Handle to a vertex in the graph.

    /**
     * Get vertices that represent extrema of a graph, in the sense that they are, on average, further away from other vertices
     * in the graph than their neighbors. The extreme vertices are required to have a minimum prominence over their
     * neighborhood. We define the "prominence" of a vertex analogously to topographic prominence
     * (http://en.wikipedia.org/wiki/Topographic_prominence) as follows. Let d(i) be the average geodesic distance of vertex i
     * to other vertices in the graph. Consider the highest path connecting d(i) to another vertex j such that d(i) <= d(j),
     * i.e. the path whose vertices have the greatest minimum d() value d'. The prominence of vertex i is then the difference
     * d(i) - d'.
     *
     * @param graph The graph whose extrema are to be found.
     * @param min_prominence Minimum prominence for a vertex to be considered an extremum.
     * @param max_results The maximum number of results to return, from the list of extrema sorted by decreasing average
     * geodesic distance.
     * @param extrema Used to return the extreme vertices. All prior data is cleared.
     *
     * @return The number of extreme vertices found.
     */
    intx getExtrema(Graph const & graph, double min_prominence, intx max_results, Array<VertexConstHandle> & extrema)
    {
      // TODO
    }

  private:
    ShortestPaths<GraphT> shortest_paths;

}; // GeodesicExtrema

} // namespace Algorithms
} // namespace Thea

#endif
