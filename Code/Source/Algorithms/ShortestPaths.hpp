//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_ShortestPaths_hpp__
#define __Thea_Algorithms_ShortestPaths_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../GraphType.hpp"
#include "../UnorderedMap.hpp"
#include "fibheap/fibheap.hpp"
#include <limits>

// #define THEA_SHORTEST_PATHS_DO_STATS
// #define THEA_SHORTEST_PATHS_TIMER

#ifdef THEA_SHORTEST_PATHS_TIMER
#  include "../Stopwatch.hpp"
#endif

// Forward declarations
struct fibheap_el;

namespace Thea {
namespace Algorithms {

/** Compute shortest paths on graphs. GraphT must satisfy IsAdjacencyGraph. */
template <typename GraphT>
class /* THEA_API */ ShortestPaths
{
  public:
    typedef GraphT Graph;  ///< The graph type.
    typedef typename GraphT::VertexHandle VertexHandle;  ///< Handle to a vertex in the graph.

    /** Stores shortest path information for a vertex. */
    class ShortestPathInfo
    {
      public:
        /** Default constructor. */
        ShortestPathInfo() {}

        /** Constructor. */
        ShortestPathInfo(double dist_, bool has_pred_, VertexHandle pred_)
        : dist(dist_), has_pred(has_pred_), pred(pred_) {}

        /** Get the distance of this vertex from the source point. */
        double getDistance() const { return dist; }

        /** Check if the vertex has a precedessor in the shortest-paths graph. */
        bool hasPredecessor() const { return pred; }

        /**
         * Get the predecessor of this vertex in the shortest path to the source. Return value is undefined is vertex does not
         * have a predecessor (see hasPredecessor()).
         */
        VertexHandle getPredecessor() const { return pred; }

      private:
        double dist;        ///< Distance to source.
        bool has_pred;      ///< Does this vertex have a predecessor in the shortest path from the source?
        VertexHandle pred;  ///< Predecessor in the shortest path from the source.

    }; // class ShortestPathInfo

  private:
    /** A callback that adds vertices to a map. */
    class MapCallback
    {
      public:
        MapCallback(TheaUnorderedMap<VertexHandle, ShortestPathInfo> & result_) : result(result_) {}

        bool operator()(VertexHandle vertex, double distance, bool has_pred, VertexHandle pred)
        {
          result[vertex] = ShortestPathInfo(distance, has_pred, pred);
          return false;
        }

      private:
        TheaUnorderedMap<VertexHandle, ShortestPathInfo> & result;

    }; // class MapCallback

  public:
    /** Destructor. */
    ~ShortestPaths()
    {
      THEA_CONCEPT_CHECK(IsAdjacencyGraph<GraphT>);
    }

    /**
     * Compute the shortest paths in a graph from a source vertex (or a set of source vertices) to every other vertex, using
     * Dijkstra's algorithm [1959]. A Fibonacci heap accelerates the computation to O(|E| + |V| log |V|)
     * [Fredman & Tarjan 1984]. <b>All edge lengths must be non-negative.</b>
     *
     * @param graph The graph to process. <b>Must have non-negative edge lengths.</b>
     * @param src The source vertex from which to measure distances. This vertex is initialized to distance zero and no
     *   predecessor. If \a src_region is non-empty, the \a src argument is ignored,
     * @param result The computed shortest-path information. Maps each vertex to the length of its shortest path from the
     *   source(s), as well its predecessor in the shortest path. Any prior data in the map is discarded.
     * @param limit If set to a non-negative value, only returns vertices that are closer than this value.
     * @param src_region If non-empty, specifies a set of source vertices and the initial distances
     *   <b>(must be non-negative)</b> to them. In this case the \a src argument is ignored. The predecessor of each such vertex
     *   is absent, unless a shorter path to the vertex is found.
     * @param include_unreachable If true, vertices unreachable from the source vertex are also returned, mapped to negative
     *   distances and without predecessors.
     */
    void dijkstra(Graph & graph, VertexHandle src, TheaUnorderedMap<VertexHandle, ShortestPathInfo> & result, double limit = -1,
                  TheaUnorderedMap<VertexHandle, double> const * src_region = NULL, bool include_unreachable = false)
    {
      result.clear();
      MapCallback callback(result);
      dijkstraWithCallback(graph, src, &callback, limit, src_region, include_unreachable);
    }

    /**
     * Compute the shortest paths in a graph from a source vertex (or a set of source vertices) to every other vertex, using
     * Dijkstra's algorithm [1959]. A Fibonacci heap accelerates the computation to O(|E| + |V| log |V|)
     * [Fredman & Tarjan 1984]. <b>All edge lengths must be non-negative.</b>
     *
     * This version of the function calls a callback operation on each processed vertex once its distance from the source has
     * been determined. The callback must be equivalent to the following function signature:
     *
     * \code
     *   //
     *   // vertex: The visited vertex in the graph.
     *   // distance: Length of shortest path to the vertex from the source.
     *   // has_pred: Does the vertex have a predecessor in the shortest path to the source?
     *   //           (False if the vertex is (part of) the source.)
     *   // pred: Predecessor of the vertex in the shortest path, if has_pred is true, else an undefined value.
     *   //
     *   bool operator()(VertexHandle vertex, double distance, bool has_pred, VertexHandle pred);
     * \endcode
     *
     * The callback should normally return false, unless it wants to terminate the search, in which case it should return true.
     *
     * @param graph The graph to process. <b>Must have non-negative edge lengths.</b>
     * @param src The source vertex from which to measure distances. This vertex is initialized to distance zero and no
     *   predecessor. If \a src_region is non-empty, the \a src argument is ignored,
     * @param callback Called for every visited vertex along with the distance to its shortest path from the source. If the
     *   function returns true, the search is terminated at this point.
     * @param limit If set to a non-negative value, only returns vertices that are closer than this value.
     * @param src_region If non-empty, specifies a set of source vertices and the initial distances
     *   <b>(must be non-negative)</b> to them. In this case the \a src argument is ignored. The predecessor of each such vertex
     *   is absent, unless a shorter path to the vertex is found.
     * @param include_unreachable If true, vertices unreachable from the source vertex are also returned, mapped to negative
     *   distances and without predecessors.
     */
    template <typename CallbackT>
    void dijkstraWithCallback(Graph & graph, VertexHandle src, CallbackT * callback, double limit = -1,
                              TheaUnorderedMap<VertexHandle, double> const * src_region = NULL,
                              bool include_unreachable = false);

  private:
    /** Status of vertex during Dijkstra traversal. */
    enum VisitStatus
    {
      BLACK,
      GREY,
      WHITE

    }; // enum VisitStatus

    /** Holds information about a vertex during Dijkstra traversal. */
    struct ScratchElement
    {
      VertexHandle vertex;
      double dist;
      bool has_pred;
      VertexHandle pred;
      VisitStatus flag;
      fibheap_el * heap_elem;

    }; // struct ScratchElement

    typedef TheaUnorderedMap<VertexHandle, ScratchElement> Scratch;  ///< Scratch space for all vertices.

    /** Compare the distances to two vertices. */
    static int compareDistances(void * vx_data1, void * vx_data2);

    Scratch scratch;  ///< Scratch space for all vertices.

};  // class ShortestPaths

template <typename GraphT>
template <typename CallbackT>
void
ShortestPaths<GraphT>::dijkstraWithCallback(Graph & graph, VertexHandle src, CallbackT * callback, double limit,
                                            TheaUnorderedMap<VertexHandle, double> const * src_region, bool include_unreachable)
{
  if (graph.numVertices() <= 0 || !callback)
    return;

#ifdef THEA_SHORTEST_PATHS_TIMER
  Stopwatch timer;
  timer.tick();
#endif

  // Sanity checks
  typedef TheaUnorderedMap<VertexHandle, double> DistanceMap;
  bool has_src_region = (src_region && !src_region->empty());
  if (has_src_region)
  {
    for (typename DistanceMap::const_iterator di = src_region->begin(); di != src_region->end(); ++di)
    {
      alwaysAssertM(di->second >= 0, "ShortestPaths: Dijkstra's algorithm requires non-negative distances")
    }
  }

  // Deallocate scratch space only if absolutely necessary. With this policy, the hash table will have at most 2.5 times the
  // number of vertices (1.5 earlier + 1 on this call) and repeated calls with the same graph will not cause any reallocations.
  long num_verts = graph.numVertices();
  if (scratch.size() > 1.5 * num_verts)
    scratch = Scratch();  // make sure the memory is actually freed, better than calling scratch.clear()

  // Create a scratch element representing negative infinity, for the Fibonacci heap
  static double const INF_DIST = (std::numeric_limits<double>::has_infinity ? std::numeric_limits<double>::infinity()
                                                                            : std::numeric_limits<double>::max());
  ScratchElement NEG_INF;
  NEG_INF.dist = -INF_DIST;

  // Initialize the heap
  fibheap * dijkstra_queue = fh_makeheap();
  fh_setcmp(dijkstra_queue, compareDistances);
  fh_setneginf(dijkstra_queue, &NEG_INF);

  // Initialize the vertex scratch data
  for (typename GraphT::VertexIterator vi = graph.verticesBegin(); vi != graph.verticesEnd(); ++vi)
  {
    VertexHandle vertex = graph.getVertex(vi);
    ScratchElement & data = scratch[vertex];  // no new allocs if the element already exists

    data.vertex = vertex;
    data.has_pred = false;

    if (has_src_region)
    {
      typename DistanceMap::const_iterator existing = src_region->find(vertex);
      if (existing != src_region->end())
      {
        data.dist = existing->second;
        data.flag = GREY;
        data.heap_elem = fh_insert(dijkstra_queue, (void *)&data);
        continue;
      }
    }
    else
    {
      if (vertex == src)
      {
        data.dist = 0;
        data.flag = GREY;
        data.heap_elem = fh_insert(dijkstra_queue, (void *)&data);
        continue;
      }
    }

    data.dist = INF_DIST;
    data.flag = WHITE;
  }

#ifdef THEA_SHORTEST_PATHS_TIMER
  timer.tock();
  THEA_CONSOLE << "ShortestPaths: Setting up scratch data and initial heap took " << 1000 * timer.elapsedTime() << "ms";
  timer.tick();
#endif

#ifdef THEA_SHORTEST_PATHS_DO_STATS
  long num_enqueued = (has_src_region ? (long)src_region->size() : 1);
  long num_iters = 0;
#endif

  ScratchElement tmp_data;
  while (true)
  {
    ScratchElement * data = (ScratchElement *)fh_extractmin(dijkstra_queue);
    if (!data)
      break;

    for (typename GraphT::NeighborIterator ni = graph.neighborsBegin(data->vertex), nbrs_end = graph.neighborsEnd(data->vertex);
         ni != nbrs_end; ++ni)
    {
      VertexHandle nbr = graph.getVertex(ni);
      double nbr_dist = graph.distance(data->vertex, ni);

      typename Scratch::iterator nbr_loc = scratch.find(nbr);
      debugAssertM(nbr_loc != scratch.end(), "ShortestPaths: No scratch entry found for neighboring vertex");

      ScratchElement * nbr_data = &nbr_loc->second;
      if (nbr_data->flag == BLACK)  // this should not be necessary, but let's do it just to be safe
        continue;

      double test_dist = nbr_dist + data->dist;
      if (test_dist < nbr_data->dist)
      {
        // Update the predecessor
        nbr_data->has_pred = true;
        nbr_data->pred = data->vertex;

        // Update the shortest-path distance and reorder the heap
        if (nbr_data->flag == WHITE)
        {
          nbr_data->dist = test_dist;
          nbr_data->flag = GREY;
          nbr_data->heap_elem = fh_insert(dijkstra_queue, (void *)nbr_data);

#ifdef THEA_SHORTEST_PATHS_DO_STATS
          num_enqueued++;
#endif
        }
        else if (nbr_data->flag == GREY)
        {
          // Hack the element to a point to a different sample, so old and new data have different pointer values (and hence
          // aren't seen as identical) but the same distance (so the heap ordering is still valid)
          tmp_data.dist = nbr_data->dist;
          fhe_setdata(nbr_data->heap_elem, (void *)&tmp_data);

          nbr_data->dist = test_dist;
          fh_replacedata(dijkstra_queue, nbr_data->heap_elem, (void *)nbr_data);
        }
      }
    }

#ifdef THEA_SHORTEST_PATHS_DO_STATS
    num_iters++;
#endif

    data->flag = BLACK;

    if (limit >= 0 && data->dist > limit)  // all remaining distances will be greater than this
      break;

    if ((*callback)(data->vertex, data->dist, data->has_pred, data->pred))
      break;
  }

  if (include_unreachable)
  {
    for (typename GraphT::VertexIterator vi = graph.verticesBegin(); vi != graph.verticesEnd(); ++vi)
    {
      VertexHandle vertex = graph.getVertex(vi);
      ScratchElement & data = scratch[vertex];  // no new allocs if the element already exists
      if (data.flag != BLACK)
      {
        if ((*callback)(data.vertex, -1, false, NULL))
          break;
      }
    }
  }

  fh_deleteheap(dijkstra_queue);

#ifdef THEA_SHORTEST_PATHS_TIMER
  timer.tock();
  THEA_CONSOLE << "ShortestPaths: Dijkstra iterations took " << 1000 * timer.elapsedTime() << "ms";
#endif

#ifdef THEA_SHORTEST_PATHS_DO_STATS
  THEA_CONSOLE << "ShortestPaths: Enqueued " << num_enqueued << " samples after " << num_iters << " Dijkstra iterations";
#endif
}

template <typename GraphT>
int
ShortestPaths<GraphT>::compareDistances(void * vx_data1, void * vx_data2)
{
  double d1 = ((ScratchElement const *)vx_data1)->dist;
  double d2 = ((ScratchElement const *)vx_data2)->dist;

  if (d1 < d2)
    return -1;
  else if (d1 > d2)
    return 1;

  return 0;
}

} // namespace Algorithms
} // namespace Thea

#endif
