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

#include "Manifold.hpp"
#include "../Polygon3.hpp"
#include "dx/simplesurf.h"
#include "dx/memory.h"
#include <cstring>

namespace Thea {
namespace Algorithms {

bool
Manifold::makeOrientedManifold(Array<Vector3> const & in_vertices, Array< Array<intx> > const & in_faces,
                               Array<Vector3> & out_vertices, Array< Array<intx> > & out_faces,
                               Array<intx> & vertex_map, Array<intx> & face_map)
{
  // Triangulate all faces, since OpenDX can only take triangles as input
  Array<int> in_triangles;
  Polygon3 poly;
  Array<intx> tri_indices;
  for (size_t i = 0; i < in_faces.size(); ++i)
  {
    poly.clear();
    for (size_t j = 0; j < in_faces[i].size(); ++j)
      poly.addVertex(in_vertices[in_faces[i][j]]);

    poly.triangulate(tri_indices);
    for (size_t j = 0; j < tri_indices.size(); ++j)
      in_triangles.push_back((int)in_faces[i][tri_indices[j]]);
  }

  // Convert the set of input vertices to a packed sequence of floats in groups of 3
  Array<float> flat_in_vertices(3 * in_vertices.size());
  for (size_t i = 0, i3 = 0; i < in_vertices.size(); ++i, i3 += 3)
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
  float * new_vertices       =  nullptr;                                    // array of generated vertices
  int   * new_triangles      =  nullptr;                                    // array of generated vertices
  int   * vertex_lut         =  nullptr;                                    // maps new vertex indices to old vertex indices
  int   * triangle_lut       =  nullptr;                                    // maps new triangle indices to old triangle indices
  int     new_num_edges      =  0;                                          // number of generated edges
  EdgeS * new_edges          =  nullptr;                                    // set of generated edges

  Htable edge_table[1];                                                     // used to internally track edges
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

    out_vertices.resize((size_t)new_num_vertices);
    for (int i = 0, i3 = 0; i < new_num_vertices; ++i, i3 += 3)
      out_vertices[i] = Vector3(new_vertices[i3], new_vertices[i3 + 1], new_vertices[i3 + 2]);

    out_faces.resize((size_t)new_num_triangles);
    for (int i = 0, i3 = 0; i < new_num_triangles; ++i, i3 += 3)
    {
      Array<intx> & face = out_faces[(size_t)i];
      face.resize(3);
      face[0] = new_triangles[i3    ];
      face[1] = new_triangles[i3 + 1];
      face[2] = new_triangles[i3 + 2];
    }

    vertex_map.resize((size_t)new_num_vertices);
    for (int i = 0; i < new_num_vertices; ++i)
      if (vertex_lut[i] >= 0)
        vertex_map[i] = (intx)vertex_lut[i];

    face_map.resize((size_t)new_num_triangles);
    for (int i = 0; i < new_num_triangles; ++i)
      face_map[(size_t)i] = (intx)triangle_lut[i];
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
