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

#include "Polygon2.hpp"
#include "Polygon3.hpp"

#ifdef SINGLE
#  define REAL float
#else
#  define REAL double
#endif

#include "Triangle/triangle.h"

#include <cmath>
#include <cstring>

namespace Thea {

Polygon2::TriangulationOptions::TriangulationOptions()
: area_bound(-1), max_steiner_points(-1)
{}

Polygon2::Polygon2()
: impl(new Polygon3)
{}

Polygon2::~Polygon2()
{
  delete impl;
}

void
Polygon2::addVertex(Vector2 const & p)
{
  impl->addVertex(Vector3(p, 0));
}

void
Polygon2::addVertex(Vector2 const & p, long index)
{
  impl->addVertex(Vector3(p, 0), index);
}

long
Polygon2::numVertices() const
{
  return impl->numVertices();
}

Polygon2::IndexedVertex
Polygon2::getVertex(long poly_index) const
{
  Polygon3::IndexedVertex const & vx3 = impl->getVertex(poly_index);
  return IndexedVertex(vx3.position.xy(), vx3.index);
}

void
Polygon2::clear()
{
  impl->clear();
}

long
Polygon2::triangulate(TheaArray<long> & tri_indices) const
{
  return impl->triangulate(tri_indices);
}

long
Polygon2::triangulateInterior(TheaArray<Vector2> & tri_verts, TheaArray<long> & tri_indices,
                              TheaArray<bool> * tri_vert_is_boundary,TriangulationOptions const & options) const
{
  struct triangulateio in, out;

  tri_verts.clear();
  tri_indices.clear();

  if (impl->vertices.size() < 3)
    return 0;

  in.numberofpointattributes = 0;
  in.pointlist = new REAL[2 * (int)impl->vertices.size()];
  in.pointmarkerlist = NULL;

  in.numberofpoints = 0;
  for (array_size_t i = 0, j = 0; i < impl->vertices.size(); ++i)
  {
    Vector3 const & p = impl->vertices[i].position;

    // Explicitly remove successive duplicate vertices from input, Triangle seems to not like them very much
    if (i <= 0 || (p - impl->vertices[i - 1].position).squaredLength() > 1.0e-10f)
    {
      in.pointlist[j++] = p.x();
      in.pointlist[j++] = p.y();
      in.numberofpoints++;
    }
  }

  if (in.numberofpoints < 3)
    return 0;

  // THEA_CONSOLE << "Triangulating polygon with " << in.numberofpoints << " vertices";
  // for (int i = 0; i < in.numberofpoints; ++i)
  //   THEA_CONSOLE << "  (" << in.pointlist[2 * i] << ", " << in.pointlist[2 * i + 1] << ")";

  in.numberofsegments = in.numberofpoints;
  in.segmentlist = new int[2 * in.numberofsegments];
  in.segmentmarkerlist = NULL;
  in.numberofholes = 0;
  in.numberofregions = 0;

  for (int i = 0, j = 1, k = 0; i < in.numberofpoints; ++i, ++j, k += 2)
  {
    if (j == in.numberofpoints)
      j = 0;

    in.segmentlist[k    ] = i;
    in.segmentlist[k + 1] = j;
  }

  out.pointlist = NULL;
  out.trianglelist = NULL;
  out.pointmarkerlist = NULL;

  std::string opt_str = format("p"         // triangulate planar straight-line graph (PSLG)
                               "q"         // quality mesh generation by Delaunay refinement, adding Steiner points
                               "a%0.32lf"  // area bound on output triangles
                               "j"         // remove unused vertices from output (e.g. duplicate vertices in input)
                               "P"         // don't output segments
                               "z"         // index everything from zero
                               "Y"         // no new vertices on boundary
                               "Q"         // quiet mode
                               , options.area_bound
                               );

  if (options.max_steiner_points >= 0)
    opt_str += format("S%ld", options.max_steiner_points);

  if (!tri_vert_is_boundary)
    opt_str += "B";  // don't output boundary markers

  char * opt_c_str = new char[opt_str.size() + 1];
  std::strcpy(opt_c_str, opt_str.c_str());

  ::triangulate(opt_c_str, &in, &out, NULL);

  // Free input arrays
  delete [] in.pointlist;
  delete [] in.segmentlist;
  delete [] opt_c_str;

  THEA_DEBUG << "Polygon2: " << in.numberofpoints << " vertices triangulated into " << out.numberofpoints << " vertices and "
             << out.numberoftriangles << " triangles";

  if (out.numberofpoints < 0 || out.numberoftriangles < 0)  // should never happen
    return 0;

  tri_verts.resize((array_size_t)out.numberofpoints);
  for (array_size_t i = 0, j = 0; i < tri_verts.size(); ++i, j += 2)
    tri_verts[i] = Vector2((Real)out.pointlist[j], (Real)out.pointlist[j + 1]);

  tri_indices.resize(3 * (array_size_t)out.numberoftriangles);
  for (array_size_t i = 0; i < tri_indices.size(); ++i)
    tri_indices[i] = (long)out.trianglelist[i];

  if (tri_vert_is_boundary)
  {
    tri_vert_is_boundary->resize((array_size_t)out.numberofpoints);
    for (array_size_t i = 0; i < tri_vert_is_boundary->size(); ++i)
      (*tri_vert_is_boundary)[i] = (out.pointmarkerlist[i] != 0);
  }

  // Free arrays allocated by Triangle
  trifree(out.pointlist);
  trifree(out.pointmarkerlist);
  trifree(out.trianglelist);

  return (long)tri_indices.size() / 3;
}

Real
Polygon2::computeArea() const
{
  return impl->computeArea();
}

} // namespace Thea
