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

#include "GeodesicSphere3.hpp"
#include "../UnorderedMap.hpp"
#include <utility>

namespace Thea {
namespace Algorithms {

namespace GeodesicSphere3Internal {

typedef std::pair<array_size_t, array_size_t> IndexPair;
typedef TheaUnorderedMap<IndexPair, array_size_t> MidpointMap;

long
midpoint(array_size_t i, array_size_t j, TheaArray<Vector3> & vertices, MidpointMap & midpoints)
{
  IndexPair edge = (i < j ? IndexPair(i, j) : IndexPair(j, i));
  MidpointMap::const_iterator existing = midpoints.find(edge);
  if (existing == midpoints.end())
  {
    Vector3 p = (vertices[i] + vertices[j]).unit();
    midpoints[edge] = vertices.size();
    vertices.push_back(p);

    return vertices.size() - 1;
  }
  else
    return existing->second;
}

void
geodesicSphereSubdivide(long depth, array_size_t i0, array_size_t i1, array_size_t i2,
                        TheaArray<Vector3> & vertices, TheaArray<long> * triangles, MidpointMap & midpoints)
{
  if (depth <= 0)
  {
    if (triangles)
    {
      triangles->push_back((long)i0);
      triangles->push_back((long)i1);
      triangles->push_back((long)i2);
    }

    return;
  }

  array_size_t m01 = midpoint(i0, i1, vertices, midpoints);
  array_size_t m12 = midpoint(i1, i2, vertices, midpoints);
  array_size_t m20 = midpoint(i2, i0, vertices, midpoints);

  geodesicSphereSubdivide(depth - 1,  i0, m01, m20, vertices, triangles, midpoints);
  geodesicSphereSubdivide(depth - 1,  i1, m12, m01, vertices, triangles, midpoints);
  geodesicSphereSubdivide(depth - 1,  i2, m20, m12, vertices, triangles, midpoints);
  geodesicSphereSubdivide(depth - 1, m01, m12, m20, vertices, triangles, midpoints);
}

} // namespace GeodesicSphere3Internal

bool
GeodesicSphere3::compute(long num_subdivs, TheaArray<Vector3> & vertices, TheaArray<long> * triangles)
{
  // http://www.opengl.org.ru/docs/pg/0208.html

  if (num_subdivs < 0)
  {
    THEA_ERROR << "GeodesicSphere3: Number of subdivisions must be non-negative";
    return false;
  }

  static Real const ICO_X = 0.525731112119133606;
  static Real const ICO_Z = 0.850650808352039932;

  static Vector3 const ICO_VERTS[12] = {
    Vector3(-ICO_X, 0.0,  ICO_Z), Vector3(ICO_X, 0.0,  ICO_Z),
    Vector3(-ICO_X, 0.0, -ICO_Z), Vector3(ICO_X, 0.0, -ICO_Z),

    Vector3(0.0,  ICO_Z, ICO_X), Vector3(0.0,  ICO_Z, -ICO_X),
    Vector3(0.0, -ICO_Z, ICO_X), Vector3(0.0, -ICO_Z, -ICO_X),

    Vector3(ICO_Z,  ICO_X, 0.0), Vector3(-ICO_Z,  ICO_X, 0.0),
    Vector3(ICO_Z, -ICO_X, 0.0), Vector3(-ICO_Z, -ICO_X, 0.0)
  };

  static array_size_t const ICO_TRIS[20][3] = {
    {0,  4,  1},  {0, 9,  4},  {9,  5, 4},  { 4, 5, 8},  {4, 8,  1},
    {8, 10,  1},  {8, 3, 10},  {5,  3, 8},  { 5, 2, 3},  {2, 7,  3},
    {7, 10,  3},  {7, 6, 10},  {7, 11, 6},  {11, 0, 6},  {0, 1,  6},
    {6,  1, 10},  {9, 0, 11},  {9, 11, 2},  { 9, 2, 5},  {7, 2, 11}
  };

  vertices.resize(12);
  for (array_size_t i = 0; i < 12; ++i)
    vertices[i] = ICO_VERTS[i];

  if (triangles && num_subdivs == 0)
  {
    triangles->resize(20 * 3);
    for (array_size_t i = 0; i < 20; ++i)
    {
      (*triangles)[3 * i    ] = (long)ICO_TRIS[i][0];
      (*triangles)[3 * i + 1] = (long)ICO_TRIS[i][1];
      (*triangles)[3 * i + 2] = (long)ICO_TRIS[i][2];
    }
  }

  if (num_subdivs > 0)
  {
    GeodesicSphere3Internal::MidpointMap midpoints;
    for (array_size_t i = 0; i < 20; ++i)
    {
      GeodesicSphere3Internal::geodesicSphereSubdivide(num_subdivs, ICO_TRIS[i][0], ICO_TRIS[i][1], ICO_TRIS[i][2],
                                                       vertices, triangles, midpoints);
    }
  }

  return true;
}

} // namespace Algorithms
} // namespace Thea
