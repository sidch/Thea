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
// First version: 2017
//
//============================================================================

#ifndef __Thea_Algorithms_GeodesicSphere3_hpp__
#define __Thea_Algorithms_GeodesicSphere3_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Generates a unit geodesic sphere, obtained by recursive subdivision of the faces of an icosahedron. This is useful for
 * getting a (nearly) uniformly separated distribution of points over the sphere.
 */
class THEA_API GeodesicSphere3
{
  public:
    /**
     * Recursively subdivides the faces of the unit icosahedron, and projects each vertex to the surface of the unit sphere.
     *
     * @param num_subdivs The number of recursive subdivisions. 0 returns the icosahedron.
     * @param vertices Used to return the vertices of the geodesic sphere.
     * @param triangles If not null, used to return the faces of the sphere, as triplets of vertex indices.
     *
     * @return True on success, false on error.
     */
    static bool compute(intx num_subdivs, Array<Vector3> & vertices, Array<intx> * triangles = nullptr);

    /**
     * Recursively subdivides a given set of triangles formed from vertices on the unit sphere, projecting each new vertex to
     * the surface of the sphere. This function is useful for generating only a part of the full geodesic sphere.
     *
     * @param num_subdivs The number of recursive subdivisions. 0 is a no-op.
     * @param vertices Used to pass the initial set of vertices, and to return the new vertices, pushed onto the end of the
     *   list.
     * @param old_triangles The existing set of triangles to subdivide, as triplets of vertex indices.
     * @param new_triangles If not null, used to return the new faces, as triplets of vertex indices.
     *
     * @return True on success, false on error.
     */
    static bool compute(intx num_subdivs, Array<Vector3> & vertices, Array<intx> const & old_triangles,
                        Array<intx> * new_triangles = nullptr);

}; // class GeodesicSphere3

} // namespace Algorithms
} // namespace Thea

#endif
