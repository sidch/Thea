//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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
//============================== Original License ============================
//
// Copyright (c) 2009 by John W. Ratcliff mailto:jratcliffscarab@gmail.com
//
// Portions of this source has been released with the PhysXViewer application, as well as
// Rocket, CreateDynamics, ODF, and as a number of sample code snippets.
//
// If you find this code useful or you are feeling particularily generous I would
// ask that you please go to http://www.amillionpixels.us and make a donation
// to Troy DeMolay.
//
// DeMolay is a youth group for young men between the ages of 12 and 21.
// It teaches strong moral principles, as well as leadership skills and
// public speaking.  The donations page uses the 'pay for pixels' paradigm
// where, in this case, a pixel is only a single penny.  Donations can be
// made for as small as $4 or as high as a $100 block.  Each person who donates
// will get a link to their own site as well as acknowledgement on the
// donations blog located here http://www.amillionpixels.blogspot.com/
//
// If you wish to contact me you can use the following methods:
//
// Skype ID: jratcliff63367
// Yahoo: jratcliff63367
// AOL: jratcliff1961
// email: jratcliffscarab@gmail.com
//
//
// The MIT license:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//============================================================================

#include "Polygon3.hpp"
#include <cmath>
#include <limits>

namespace Thea {

Polygon3::Polygon3()
: max_index(-1)
{}

void
Polygon3::clear()
{
  vertices.clear();
  max_index = -1;
  bounds = AxisAlignedBox3();
}

void
Polygon3::addVertex(Vector3 const & p)
{
  addVertex(p, max_index + 1);
}

void
Polygon3::addVertex(Vector3 const & p, long index)
{
  // update bounding box...
  bounds.merge(p);

  vertices.push_back(IndexedVertex(p, index));
  if (index > max_index) max_index = index;
}

long
Polygon3::numVertices() const
{
  return (long)vertices.size();
}

Polygon3::IndexedVertex const &
Polygon3::getVertex(long poly_index) const
{
  debugAssertM(poly_index > 0 && poly_index < numVertices(), "Polygon3: Vertex index out of bounds");

  return vertices[(array_size_t)poly_index];
}

static bool
Polygon3_insideTriangle2(Vector2 const & A, Vector2 const & B, Vector2 const & C, Vector2 const & P)
{
  Real ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
  Real cCROSSap, bCROSScp, aCROSSbp;

  ax  = C.x() - B.x();  ay  = C.y() - B.y();
  bx  = A.x() - C.x();  by  = A.y() - C.y();
  cx  = B.x() - A.x();  cy  = B.y() - A.y();
  apx = P.x() - A.x();  apy = P.y() - A.y();
  bpx = P.x() - B.x();  bpy = P.y() - B.y();
  cpx = P.x() - C.x();  cpy = P.y() - C.y();

  aCROSSbp = ax * bpy - ay * bpx;
  cCROSSap = cx * apy - cy * apx;
  bCROSScp = bx * cpy - by * cpx;

  return ((aCROSSbp >= 0) && (bCROSScp >= 0) && (cCROSSap >= 0));
}

bool
Polygon3::snip(array_size_t u, array_size_t v, array_size_t w, array_size_t n, TheaArray<array_size_t> const & indices) const
{
  static Real const epsilon = 1e-10f;

  Vector2 const & A = proj_vertices[indices[u]];
  Vector2 const & B = proj_vertices[indices[v]];
  Vector2 const & C = proj_vertices[indices[w]];

  if (epsilon > (((B.x() - A.x()) * (C.y() - A.y())) - ((B.y() - A.y()) * (C.x() - A.x()))))
    return false;

  for (array_size_t p = 0; p < n; ++p)
  {
    if ((p == u) || (p == v) || (p == w))
      continue;

    Vector2 const & P = proj_vertices[indices[p]];
    if (Polygon3_insideTriangle2(A, B, C, P))
      return false;
  }

  return true;
}

// Triangulation happens in 2d.  We could inverse transform the polygon around the normal direction, or we just use the two most
// signficant axes. Here we find the two longest axes and use them to triangulate.  Inverse transforming them would introduce
// more doubling point error and isn't worth it.
long
Polygon3::triangulate(TheaArray<long> & tri_indices) const
{
  if (vertices.size() < 3)
  {
    tri_indices.clear();
  }
  else if (vertices.size() == 3)
  {
    tri_indices.resize(3);
    tri_indices[0] = vertices[0].index;
    tri_indices[1] = vertices[1].index;
    tri_indices[2] = vertices[2].index;
  }
  else if (vertices.size() > 3)
  {
    tri_indices.clear();

    // Compute a permutation of vertices that put the most significant axes of the triangle in X and Y
    Vector3 extent = bounds.getExtent();
    int i1 = 0, i2 = 0, i3 = 0;
    if (extent[0] > extent[1] && extent[0] > extent[2])
    {
      i1 = 0;
      if (extent[1] > extent[2])
      {
        i2 = 1;
        i3 = 2;
      }
      else
      {
        i2 = 2;
        i3 = 1;
      }
    }
    else if (extent[1] > extent[0] && extent[1] > extent[2])
    {
      i1 = 1;
      if (extent[0] > extent[2])
      {
        i2 = 0;
        i3 = 2;
      }
      else
      {
        i2 = 2;
        i3 = 0;
      }
    }
    else
    {
      i1 = 2;
      if (extent[0] > extent[1])
      {
        i2 = 0;
        i3 = 1;
      }
      else
      {
        i2 = 1;
        i3 = 0;
      }
    }

    array_size_t n = vertices.size();
    proj_vertices.resize(n);
    for (array_size_t i = 0; i < n; ++i)
    {
      Vector3 const & v = vertices[i].position;
      proj_vertices[i] = Vector2(v[i1], v[i2]);
    }

    TheaArray<array_size_t> indices(n);
    bool flipped = false;
    if (projArea() > 0)
    {
      for (array_size_t v = 0; v < n; ++v)
        indices[v] = v;
    }
    else
    {
      for (array_size_t v = 0; v < n; ++v)
        indices[v] = (n - 1) - v;

      flipped = true;
    }

    array_size_t nv = n;
    array_size_t count = 2 * nv;
    for (array_size_t v = nv - 1; nv > 2; )
    {
      if ((count--) <= 0)
        break;

      array_size_t u = v;
      if (nv <= u) u = 0;

      v = u + 1;
      if (nv <= v) v = 0;

      array_size_t w = v + 1;
      if (nv <= w) w = 0;

      if (snip(u, v, w, nv, indices))
      {
        array_size_t a = indices[u];
        array_size_t b = indices[v];
        array_size_t c = indices[w];
        if (flipped)
        {
          tri_indices.push_back(vertices[c].index);
          tri_indices.push_back(vertices[b].index);
          tri_indices.push_back(vertices[a].index);
        }
        else
        {
          tri_indices.push_back(vertices[a].index);
          tri_indices.push_back(vertices[b].index);
          tri_indices.push_back(vertices[c].index);
        }

        array_size_t s = v, t = v + 1;
        for ( ; t < nv; ++s, ++t)
          indices[s] = indices[t];

        nv--;
        count = 2 * nv;
      }
    }
  }

  return (long)tri_indices.size() / 3;
}

Real
Polygon3::projArea() const
{
  array_size_t n = proj_vertices.size();
  Real a = 0;
  for (array_size_t p = n - 1, q = 0; q < n; p = q++)
  {
    Vector2 const & pval = proj_vertices[p];
    Vector2 const & qval = proj_vertices[q];
    a += pval.x() * qval.y() - qval.x() * pval.y();
  }

  return 0.5f * a;
}

Real
Polygon3::area() const
{
  array_size_t n = vertices.size();
  if (n < 3)
    return 0;

  Vector3 normal = (vertices[2].position - vertices[1].position).cross(vertices[0].position - vertices[1].position);
  if (vertices.size() == 3)
    return 0.5f * normal.length();

  normal.unitize();

  Real a = 0;
  Vector3 const & last = vertices[n - 1].position;
  Vector3 prev = last;
  for (array_size_t p = 0; p < n; ++p)
  {
    Vector3 const & v = vertices[p].position;
    a += (prev.cross(v)).dot(normal);
    prev = v;
  }

  return std::fabs(0.5f * a);
}

} // namespace Thea
