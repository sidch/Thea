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
// First version: 2009
//
//============================================================================

#ifndef __Thea_Algorithms_Manifold_hpp__
#define __Thea_Algorithms_Manifold_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../MatVec.hpp"
#include "../Set.hpp"
#include "../Stack.hpp"
#include "../UnionFind.hpp"
#include <utility>

namespace Thea {
namespace Algorithms {

/** Adjust mesh-based surfaces to have manifold topology. */
class THEA_API Manifold
{
  public:
    /**
     * Apply the cutting-and-stitching algorithm from A. Gueziec and G. Taubin, "Cutting and Stitching: Converting Sets of
     * Polygons to Manifold Surfaces", IEEE Trans. Vis. Comp. Graph., 7(2), pp. 136--151, 2001. Currently only the "cutting"
     * part is implemented (via IBM's OpenDX source code).
     */
    static bool makeOrientedManifold(Array<Vector3> const & in_vertices, Array< Array<intx> > const & in_faces,
                                     Array<Vector3> & out_vertices, Array< Array<intx> > & out_faces,
                                     Array<intx> & vertex_map, Array<intx> & face_map);

    /**
     * Convert a non-manifold mesh to manifold form. Each vertex with a non-manifold neighborhood is split up into copies with
     * manifold neighborhoods. Note that this function is not guaranteed to work. In particular, it throws an error if any face
     * has a repeated vertex. Such faces are easy to detect and should be removed from the list of input faces before this
     * function is called.
     *
     * FaceT should be a sequence holding values of type FaceT::value_type.
     *
     * @param faces Each entry is the sequence of vertex indices of the corresponding face. These indices are updated by the
     *          function when copying vertices.
     * @param num_vertices Total number of vertices.
     * @param vertex_map Maps output vertices (including ones created as copies of original vertices) to the set of original
     *          vertices
     *
     * @return True if some vertices needed to be duplicated, false if no changes were required.
     */
    template <typename FaceT>
    static bool makeManifold(Array<FaceT> & faces, intx num_vertices, Array<intx> & vertex_map)
    {
      if (num_vertices <= 0 || faces.size() <= 0)
        return false;

      size_t nf = faces.size();
      size_t nv = (size_t)num_vertices;

      // Initialize every vertex to map to itself
      vertex_map.resize(nv);
      for (size_t i = 0; i < nv; ++i)
        vertex_map[i] = (intx)i;

      // Associate each vertex with its incident faces
      Array< Array<size_t> > v2f(nv);
      for (size_t i = 0; i < nf; ++i)
        for (typename FaceT::const_iterator vi = faces[i].begin(); vi != faces[i].end(); ++vi)
        {
          debugAssertM(*vi >= 0 && (intx)*vi < num_vertices,
                      format("Manifold: Vertex index %ld out of range [0, %ld)", (intx)*vi, num_vertices));
          v2f[(size_t)*vi].push_back(i);
        }

      // Queue the set of vertices
      Stack<size_t> vertex_stack;
      for (intx i = (intx)nv - 1; i >= 0; --i)
        vertex_stack.push((size_t)i);

      // Split the faces at each vertex into manifold groups
      while (!vertex_stack.empty())
      {
        size_t i = vertex_stack.top();
        vertex_stack.pop();

        // Set of edges (represented by their further vertices) incident at this vertex that have already been observed to be
        // shared by two faces
        Set<size_t> shared_edges;

        // Group the faces into maximal edge-connected components
        typedef UnionFind<size_t> UnionFind;
        UnionFind uf(v2f[i].begin(), v2f[i].end());

        for (size_t j = 0; j < v2f[i].size(); ++j)
          for (size_t k = j + 1; k < v2f[i].size(); ++k)
            if (shareEdgeAtVertex(faces[v2f[i][j]], faces[v2f[i][k]], i, shared_edges))
              uf.merge((intx)j, (intx)k);

        // Retain only the faces edge-connected to the first one, assigning the rest to a copy of the vertex
        bool created_new_vertex = false;
        for (size_t j = 1; j < v2f[i].size(); ++j)
        {
          if (!uf.sameSet(0, j))
          {
            size_t face = v2f[i][j];
            size_t copy_index = v2f.size() - 1;

            if (!created_new_vertex) // create a copy of the vertex
            {
              copy_index++;
              v2f.push_back(Array<size_t>());
              vertex_map.push_back(vertex_map[i]);
              vertex_stack.push(copy_index);

              created_new_vertex = true;
            }

            v2f[copy_index].push_back(face);

            // Update the vertex index of this face to point to the copy
            for (typename FaceT::iterator vi = faces[face].begin(); vi != faces[face].end(); ++vi)
              if (*vi == (typename FaceT::value_type)i)
              {
                *vi = (typename FaceT::value_type)copy_index;
                break;
              }
          }
        }

        // If we made a copy of this vertex, we can get rid of faces reassigned to the copy
        if (created_new_vertex)
        {
          Array<size_t> nbd;
          for (size_t j = 0; j < v2f[i].size(); ++j)
            if (uf.sameSet(0, j))
              nbd.push_back(v2f[i][j]);

          v2f[i] = nbd;
        }
      }

      if ((intx)vertex_map.size() > num_vertices)
      {
        THEA_DEBUG << "Manifold: Non-manifold vertices found and fixed";
        return true;
      }
      else
        return false;
    }

  private:
    /** Get the last index of a face. */
    template <typename FaceT> static size_t getLastIndex(FaceT const & face)
    { return (size_t)*face.rbegin(); }

    /**
     * Get the two neighboring vertices of a given vertex on the boundary of a face. Throws an error if the face has repeated
     * vertices.
     */
    template <typename FaceT>
    static bool getNeighboringVertices(FaceT const & face, size_t vertex, size_t & nbr1, size_t & nbr2)
    {
      // Find the neighbors of the vertex in the first face
      size_t last_vertex = getLastIndex(face);
      size_t next_vertex;
      bool found = false;
      for (typename FaceT::const_iterator vi = face.begin(); vi != face.end(); ++vi)
      {
        next_vertex = (size_t)(vi + 1 == face.end() ? *face.begin() : *(vi + 1));
        size_t this_vertex = (size_t)*vi;

        if (this_vertex == next_vertex || last_vertex == next_vertex)
          throw Error("Manifold: Face has repeated vertices");

        if (this_vertex == vertex)
        {
          found = true;
          break;
        }
        last_vertex = this_vertex;
      }

      if (!found)
        return false;

      nbr1 = (size_t)last_vertex;
      nbr2 = (size_t)next_vertex;
      return true;
    }

    /**
     * Check if two faces share an edge that is incident on a given vertex, and is not already in a set of shared edges. If the
     * condition holds, then all such shared edges are added to the set of shared edges.
     */
    template <typename FaceT>
    static bool shareEdgeAtVertex(FaceT const & face1, FaceT const & face2, size_t vertex,
                                  Set<size_t> & shared_edges)
    {
      size_t u1, u2, w1, w2;
      if (!getNeighboringVertices(face1, vertex, u1, u2))
        throw Error("Manifold: Vertex does not belong to first face");

      if (!getNeighboringVertices(face2, vertex, w1, w2))
        throw Error("Manifold: Vertex does not belong to second face");

      bool ret = false;
      if ((u1 == w1 || u1 == w2) && shared_edges.find(u1) == shared_edges.end())
      {
        ret = true;
        shared_edges.insert(u1);
      }

      if ((u2 == w1 || u2 == w2) && shared_edges.find(u2) == shared_edges.end())
      {
        ret = true;
        shared_edges.insert(u2);
      }

      return ret;
    }

}; // class Manifold

} // namespace Algorithms
} // namespace Thea

#endif
