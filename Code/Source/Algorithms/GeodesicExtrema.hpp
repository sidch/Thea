//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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
    long getExtrema(Graph const & graph, double min_prominence, long max_results, TheaArray<VertexConstHandle> & extrema)
    {
    }

  private:
    ShortestPaths<GraphT> shortest_paths;

}; // GeodesicExtrema

} // namespace Algorithms
} // namespace Thea

#endif
