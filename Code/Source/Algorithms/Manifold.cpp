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

#include "Manifold.hpp"
#include "../Polygon3.hpp"
#include "dx/simplesurf.h"
#include "dx/memory.h"
#include <cstring>

namespace Thea {
namespace Algorithms {

bool
Manifold::makeOrientedManifold(TheaArray<Vector3> const & in_vertices,
                               TheaArray< TheaArray<long> > const & in_faces, TheaArray<Vector3> & out_vertices,
                               TheaArray< TheaArray<long> > & out_faces, TheaArray<long> & vertex_map,
                               TheaArray<long> & face_map)
{
  // Triangulate all faces, since OpenDX can only take triangles as input
  TheaArray<int> in_triangles;
  Polygon3 poly;
  TheaArray<long> tri_indices;
  for (array_size_t i = 0; i < in_faces.size(); ++i)
  {
    poly.clear();
    for (array_size_t j = 0; j < in_faces[i].size(); ++j)
      poly.addVertex(in_vertices[in_faces[i][j]]);

    poly.triangulate(tri_indices);
    for (array_size_t j = 0; j < tri_indices.size(); ++j)
      in_triangles.push_back((int)in_faces[i][tri_indices[j]]);
  }

  // Convert the set of input vertices to a packed sequence of floats in groups of 3
  TheaArray<float> flat_in_vertices(3 * in_vertices.size());
  for (array_size_t i = 0, i3 = 0; i < in_vertices.size(); ++i, i3 += 3)
  {
    flat_in_vertices[i3    ] = in_vertices[i].x();
    flat_in_vertices[i3 + 1] = in_vertices[i].y();
    flat_in_vertices[i3 + 2] = in_vertices[i].z();
  }

  // Call the OpenDX function to convert the input mesh to a collection of oriented manifolds
  int     num_vertices       =  (int)in_vertices.size();                    // number of input vertices
  int     num_triangles      =  (int)(in_triangles.size() / 3);             // number of input vertices
  float * in_v               =  const_cast<float *>(&flat_in_vertices[0]);  // array of input vertices
  int   * in_t               =  const_cast<int *>(&in_triangles[0]);        // array of input triangles
  int     new_num_vertices   =  0;                                          // number of generated vertices
  int     new_num_triangles  =  0;                                          // number of generated triangles
  float * new_vertices       =  NULL;                                       // array of generated vertices
  int   * new_triangles      =  NULL;                                       // array of generated vertices
  int   * vertex_lut         =  NULL;                                       // maps new vertex indices to old vertex indices
  int   * triangle_lut       =  NULL;                                       // maps new triangle indices to old triangle indices
  int     new_num_edges      =  0;                                          // number of generated edges
  EdgeS * new_edges          =  NULL;                                       // set of generated edges

  Htable edge_table[1];                                                    // used to internally track edges
  std::memset((void *)edge_table, 0, sizeof(Htable));

  int success = _dxfToOrientedManifold(num_vertices,
                                       in_v,
                                       num_triangles,
                                       in_t,
                                       &new_num_vertices,
                                       &new_vertices,
                                       &new_num_triangles,
                                       &new_triangles,
                                       &vertex_lut,
                                       &triangle_lut,
                                       &new_num_edges,
                                       edge_table,
                                       &new_edges);

  if (success)
  {
#ifdef THEA_DEBUG_BUILD
    THEA_DEBUG << "Manifold: Conversion to oriented manifold succeeded";
    THEA_DEBUG << "    Input: " << num_vertices << " vertices, " << in_faces.size() << " faces (" << num_triangles
               << " triangles)";
    THEA_DEBUG << "    Output: " << new_num_vertices << " vertices, " << new_num_triangles << " triangles";
#endif

    out_vertices.resize((array_size_t)new_num_vertices);
    for (int i = 0, i3 = 0; i < new_num_vertices; ++i, i3 += 3)
      out_vertices[i] = Vector3(new_vertices[i3], new_vertices[i3 + 1], new_vertices[i3 + 2]);

    out_faces.resize((array_size_t)new_num_triangles);
    for (int i = 0, i3 = 0; i < new_num_triangles; ++i, i3 += 3)
    {
      TheaArray<long> & face = out_faces[(array_size_t)i];
      face.resize(3);
      face[0] = new_triangles[i3    ];
      face[1] = new_triangles[i3 + 1];
      face[2] = new_triangles[i3 + 2];
    }

    vertex_map.resize((array_size_t)new_num_vertices);
    for (int i = 0; i < new_num_vertices; ++i)
      if (vertex_lut[i] >= 0)
        vertex_map[i] = (long)vertex_lut[i];

    face_map.resize((array_size_t)new_num_triangles);
    for (int i = 0; i < new_num_triangles; ++i)
      face_map[(array_size_t)i] = (long)triangle_lut[i];
  }

  // Cleanup
  if (new_vertices  != in_v) DXFree((Pointer)new_vertices);
  if (new_triangles != in_t) DXFree((Pointer)new_triangles);
  DXFree((Pointer)vertex_lut);
  DXFree((Pointer)triangle_lut);
  _dxfFreeHashTable(edge_table);
  DXFree((Pointer)new_edges);

  return success != 0;
}

} // namespace Algorithms
} // namespace Thea
