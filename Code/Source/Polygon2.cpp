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

#include "Polygon2.hpp"
#include "Polygon3.hpp"

// Be careful not to use NumericType::REAL after this!!!
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
  impl->addVertex(Vector3(p[0], p[1], 0));
}

void
Polygon2::addVertex(Vector2 const & p, intx index)
{
  impl->addVertex(Vector3(p[0], p[1], 0), index);
}

intx
Polygon2::numVertices() const
{
  return impl->numVertices();
}

Polygon2::IndexedVertex
Polygon2::getVertex(intx poly_index) const
{
  Polygon3::IndexedVertex const & vx3 = impl->getVertex(poly_index);
  return IndexedVertex(vx3.position.head<2>(), vx3.index);
}

void
Polygon2::clear()
{
  impl->clear();
}

intx
Polygon2::triangulate(Array<intx> & tri_indices) const
{
  return impl->triangulate(tri_indices);
}

intx
Polygon2::triangulateInterior(Array<Vector2> & tri_verts, Array<intx> & tri_indices,
                              Array<bool> * tri_vert_is_boundary,TriangulationOptions const & options) const
{
  struct triangulateio in, out;

  tri_verts.clear();
  tri_indices.clear();

  if (impl->vertices.size() < 3)
    return 0;

  in.numberofpointattributes = 0;
  in.pointlist = new REAL[2 * (int)impl->vertices.size()];
  in.pointmarkerlist = nullptr;

  in.numberofpoints = 0;
  for (size_t i = 0, j = 0; i < impl->vertices.size(); ++i)
  {
    Vector3 const & p = impl->vertices[i].position;

    // Explicitly remove successive duplicate vertices from input, Triangle seems to not like them very much
    if (i <= 0 || (p - impl->vertices[i - 1].position).squaredNorm() > 1.0e-10f)
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
  in.segmentmarkerlist = nullptr;
  in.numberofholes = 0;
  in.numberofregions = 0;

  for (int i = 0, j = 1, k = 0; i < in.numberofpoints; ++i, ++j, k += 2)
  {
    if (j == in.numberofpoints)
      j = 0;

    in.segmentlist[k    ] = i;
    in.segmentlist[k + 1] = j;
  }

  out.pointlist = nullptr;
  out.trianglelist = nullptr;
  out.pointmarkerlist = nullptr;

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

  ::triangulate(opt_c_str, &in, &out, nullptr);

  // Free input arrays
  delete [] in.pointlist;
  delete [] in.segmentlist;
  delete [] opt_c_str;

  THEA_DEBUG << "Polygon2: " << in.numberofpoints << " vertices triangulated into " << out.numberofpoints << " vertices and "
             << out.numberoftriangles << " triangles";

  if (out.numberofpoints < 0 || out.numberoftriangles < 0)  // should never happen
    return 0;

  tri_verts.resize((size_t)out.numberofpoints);
  for (size_t i = 0, j = 0; i < tri_verts.size(); ++i, j += 2)
    tri_verts[i] = Vector2((Real)out.pointlist[j], (Real)out.pointlist[j + 1]);

  tri_indices.resize(3 * (size_t)out.numberoftriangles);
  for (size_t i = 0; i < tri_indices.size(); ++i)
    tri_indices[i] = (intx)out.trianglelist[i];

  if (tri_vert_is_boundary)
  {
    tri_vert_is_boundary->resize((size_t)out.numberofpoints);
    for (size_t i = 0; i < tri_vert_is_boundary->size(); ++i)
      (*tri_vert_is_boundary)[i] = (out.pointmarkerlist[i] != 0);
  }

  // Free arrays allocated by Triangle
  trifree(out.pointlist);
  trifree(out.pointmarkerlist);
  trifree(out.trianglelist);

  return (intx)tri_indices.size() / 3;
}

Real
Polygon2::computeArea() const
{
  return impl->computeArea();
}

} // namespace Thea
