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
//============================================================================

#ifndef __Thea_Polygon2_hpp__
#define __Thea_Polygon2_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "Vector2.hpp"

namespace Thea {

class Polygon3;

/** A polygon in 2-space. */
class THEA_API Polygon2
{
  public:
    THEA_DEF_POINTER_TYPES(Polygon2, shared_ptr, weak_ptr)

    /** %Options controlling interior triangulation. */
    struct THEA_API TriangulationOptions
    {
      Real area_bound;  ///< The maximum area of any output triangle (default -1, indicating no upper bound).
      long max_steiner_points;  /**< The maximum number of Steiner points that can be inserted (default -1, indicating no
                                     limit). */

      /** Constructor. */
      TriangulationOptions();

      /** Get the default set of triangulation options. */
      static TriangulationOptions const & defaults() { static TriangulationOptions const def; return def; }

    }; // struct TriangulationOptions

    /** A vertex plus an index. */
    struct THEA_API IndexedVertex
    {
      /** Default constructor. */
      IndexedVertex() {}

      /** Initializing constructor. */
      IndexedVertex(Vector2 const & position_, long index_) : position(position_), index(index_) {}

      Vector2 position;  ///< The position of the vertex.
      long index;  ///< The index of the vertex.
    };

    /** Construct an empty polygon. */
    Polygon2();

    /** Destructor. */
    ~Polygon2();

    /**
     * Add a vertex to the polygon. The vertex is inserted at the end of the current sequence of vertices, and by default is
     * assigned an index that is one more than the maximum index in the polygon so far (or zero if this is the first vertex).
     */
    void addVertex(Vector2 const & p);

    /** Add an indexed vertex to the polygon. The vertex is inserted at the end of the current sequence of vertices. */
    void addVertex(Vector2 const & p, long index);

    /** Get the number of vertices in the polygon. */
    long numVertices() const;

    /**
     * Get the vertex at position \a poly_index in the sequence of vertices around the polygon boundary.
     *
     * @note \a poly_index is determined by the sequence of addVertex() calls, <b>NOT</b> by the index supplied in
     *   addVertex(Vector2 const &, long)!
     */
    IndexedVertex getVertex(long poly_index) const;

    /** Delete all vertices from the polygon. */
    void clear();

    /**
     * Triangulate the polygon and return the set of triangle indices (in successive groups of 3). All prior data in the
     * supplied array are cleared.
     *
     * @return The number of triangles created.
     */
    long triangulate(TheaArray<long> & tri_indices) const;

    /**
     * Triangulate the polygon, inserting Steiner vertices as necessary in the interior of the polygon for a well-conditioned
     * result. All prior data in the supplied arrays are cleared.
     *
     * @param tri_verts Used to return the vertices of the output triangulation.
     * @param tri_indices Used to return the vertex indices of output triangles (w.r.t. \a tri_verts), in successive groups of
     *   3.
     * @param tri_vert_is_boundary If not-null, used to return a flag for each output vertex indicating whether it's on the
     *   boundary of the triangulated domain or not.
     * @param options %Options controlling the quality of the triangulation.
     *
     * @return The number of triangles created.
     */
    long triangulateInterior(TheaArray<Vector2> & tri_verts, TheaArray<long> & tri_indices,
                             TheaArray<bool> * tri_vert_is_boundary = NULL,
                             TriangulationOptions const & options = TriangulationOptions::defaults()) const;

    /** Compute the area of the polygon. */
    Real computeArea() const;

  private:
    Polygon3 * impl;

}; // class Polygon2

} // namespace Thea

#endif
