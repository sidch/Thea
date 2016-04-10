//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_Manifold_hpp__
#define __Thea_Algorithms_Manifold_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Set.hpp"
#include "../Stack.hpp"
#include "../UnionFind.hpp"
#include "../Vector3.hpp"
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
    static bool makeOrientedManifold(TheaArray<Vector3> const & in_vertices,
                                     TheaArray< TheaArray<long> > const & in_faces,
                                     TheaArray<Vector3> & out_vertices, TheaArray< TheaArray<long> > & out_faces,
                                     TheaArray<long> & vertex_map, TheaArray<long> & face_map);

    /**
     * Convert a non-manifold mesh to manifold form. Each vertex with a non-manifold neighbourhood is split up into copies with
     * manifold neighbourhoods. Note that this function is not guaranteed to work. In particular, it throws an error if any face
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
    static bool makeManifold(TheaArray<FaceT> & faces, long num_vertices, TheaArray<long> & vertex_map)
    {
      if (num_vertices <= 0 || faces.size() <= 0)
        return false;

      array_size_t nf = (array_size_t)faces.size();
      array_size_t nv = (array_size_t)num_vertices;

      // Initialize every vertex to map to itself
      vertex_map.resize(nv);
      for (array_size_t i = 0; i < nv; ++i)
        vertex_map[i] = (long)i;

      // Associate each vertex with its incident faces
      TheaArray< TheaArray<array_size_t> > v2f(nv);
      for (array_size_t i = 0; i < nf; ++i)
        for (typename FaceT::const_iterator vi = faces[i].begin(); vi != faces[i].end(); ++vi)
        {
          debugAssertM(*vi >= 0 && (long)*vi < num_vertices,
                      format("Manifold: Vertex index %ld out of range [0, %ld)", (long)*vi, num_vertices));
          v2f[(array_size_t)*vi].push_back(i);
        }

      // Queue the set of vertices
      TheaStack<array_size_t> vertex_stack;
      for (long i = (long)nv - 1; i >= 0; --i)
        vertex_stack.push((array_size_t)i);

      // Split the faces at each vertex into manifold groups
      while (!vertex_stack.empty())
      {
        array_size_t i = vertex_stack.top();
        vertex_stack.pop();

        // Set of edges (represented by their further vertices) incident at this vertex that have already been observed to be
        // shared by two faces
        TheaSet<array_size_t> shared_edges;

        // Group the faces into maximal edge-connected components
        typedef UnionFind<array_size_t> UnionFind;
        UnionFind uf(v2f[i].begin(), v2f[i].end());

        for (array_size_t j = 0; j < v2f[i].size(); ++j)
          for (array_size_t k = j + 1; k < v2f[i].size(); ++k)
            if (shareEdgeAtVertex(faces[v2f[i][j]], faces[v2f[i][k]], i, shared_edges))
              uf.merge((long)j, (long)k);

        // Retain only the faces edge-connected to the first one, assigning the rest to a copy of the vertex
        bool created_new_vertex = false;
        for (array_size_t j = 1; j < v2f[i].size(); ++j)
        {
          if (!uf.sameSet(0, j))
          {
            array_size_t face = v2f[i][j];
            array_size_t copy_index = v2f.size() - 1;

            if (!created_new_vertex) // create a copy of the vertex
            {
              copy_index++;
              v2f.push_back(TheaArray<array_size_t>());
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
          TheaArray<array_size_t> nbd;
          for (array_size_t j = 0; j < v2f[i].size(); ++j)
            if (uf.sameSet(0, j))
              nbd.push_back(v2f[i][j]);

          v2f[i] = nbd;
        }
      }

      if ((long)vertex_map.size() > num_vertices)
      {
        THEA_DEBUG << "Manifold: Non-manifold vertices found and fixed";
        return true;
      }
      else
        return false;
    }

  private:
    /** Get the last index of a face. */
    template <typename FaceT> static array_size_t getLastIndex(FaceT const & face)
    { return (array_size_t)*face.rbegin(); }

    /**
     * Get the two neighbouring vertices of a given vertex on the boundary of a face. Throws an error if the face has repeated
     * vertices.
     */
    template <typename FaceT>
    static bool getNeighbouringVertices(FaceT const & face, array_size_t vertex, array_size_t & nbr1, array_size_t & nbr2)
    {
      // Find the neighbours of the vertex in the first face
      array_size_t last_vertex = getLastIndex(face);
      array_size_t next_vertex;
      bool found = false;
      for (typename FaceT::const_iterator vi = face.begin(); vi != face.end(); ++vi)
      {
        next_vertex = (array_size_t)(vi + 1 == face.end() ? *face.begin() : *(vi + 1));
        array_size_t this_vertex = (array_size_t)*vi;

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

      nbr1 = (array_size_t)last_vertex;
      nbr2 = (array_size_t)next_vertex;
      return true;
    }

    /**
     * Check if two faces share an edge that is incident on a given vertex, and is not already in a set of shared edges. If the
     * condition holds, then all such shared edges are added to the set of shared edges.
     */
    template <typename FaceT>
    static bool shareEdgeAtVertex(FaceT const & face1, FaceT const & face2, array_size_t vertex,
                                  TheaSet<array_size_t> & shared_edges)
    {
      array_size_t u1, u2, w1, w2;
      if (!getNeighbouringVertices(face1, vertex, u1, u2))
        throw Error("Manifold: Vertex does not belong to first face");

      if (!getNeighbouringVertices(face2, vertex, w1, w2))
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

};  // class Manifold

} // namespace Algorithms
} // namespace Thea

#endif
