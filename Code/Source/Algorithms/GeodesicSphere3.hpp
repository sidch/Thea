//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2017, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_GeodesicSphere3_hpp__
#define __Thea_Algorithms_GeodesicSphere3_hpp__

#include "../Common.hpp"
#include "../Vector3.hpp"

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
    static bool compute(long num_subdivs, TheaArray<Vector3> & vertices, TheaArray<long> * triangles = NULL);

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
    static bool compute(long num_subdivs, TheaArray<Vector3> & vertices, TheaArray<long> const & old_triangles,
                        TheaArray<long> * new_triangles = NULL);

}; // class GeodesicSphere3

} // namespace Algorithms
} // namespace Thea

#endif
