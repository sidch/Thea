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
//============================================================================

#ifndef __Thea_Polygon2_hpp__
#define __Thea_Polygon2_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "MatVec.hpp"

namespace Thea {

class Polygon3;

/** A polygon in 2-space. */
class THEA_API Polygon2
{
  public:
    THEA_DECL_SMART_POINTERS(Polygon2)

    /** %Options controlling interior triangulation. */
    struct THEA_API TriangulationOptions
    {
      Real area_bound;  ///< The maximum area of any output triangle (default -1, indicating no upper bound).
      intx max_steiner_points;  /**< The maximum number of Steiner points that can be inserted (default -1, indicating no
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
      IndexedVertex(Vector2 const & position_, intx index_) : position(position_), index(index_) {}

      Vector2 position;  ///< The position of the vertex.
      intx index;  ///< The index of the vertex.
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
    void addVertex(Vector2 const & p, intx index);

    /** Get the number of vertices in the polygon. */
    intx numVertices() const;

    /**
     * Get the vertex at position \a poly_index in the sequence of vertices around the polygon boundary.
     *
     * @note \a poly_index is determined by the sequence of addVertex() calls, <b>NOT</b> by the index supplied in
     *   addVertex(Vector2 const &, intx)!
     */
    IndexedVertex getVertex(intx poly_index) const;

    /** Delete all vertices from the polygon. */
    void clear();

    /**
     * Triangulate the polygon and return the set of triangle indices (in successive groups of 3). All prior data in the
     * supplied array are cleared.
     *
     * @return The number of triangles created.
     */
    intx triangulate(Array<intx> & tri_indices) const;

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
    intx triangulateInterior(Array<Vector2> & tri_verts, Array<intx> & tri_indices,
                             Array<bool> * tri_vert_is_boundary = nullptr,
                             TriangulationOptions const & options = TriangulationOptions::defaults()) const;

    /** Compute the area of the polygon. */
    Real computeArea() const;

  private:
    Polygon3 * impl;

}; // class Polygon2

} // namespace Thea

#endif
