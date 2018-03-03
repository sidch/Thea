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

#include "ConvexHull3.hpp"
#include "../Polygon3.hpp"
#include "StanHull/hull.h"

namespace Thea {
namespace Algorithms {

ConvexHull3::ConvexHull3(Options const & options_)
: options(options_), approx_updated(true), exact_updated(true)
{}

void
ConvexHull3::addPoint(Vector3 const & point)
{
  points.push_back(point);
  approx_updated = exact_updated = false;
}

void
ConvexHull3::clear()
{
  points.clear();
  approx_updated = exact_updated = false;
}

void
ConvexHull3::releaseMemoryWithoutUpdate()
{
  points.clear();
}

void
ConvexHull3::updateApprox() const
{
  if (approx_updated)
    return;

  THEA_DEBUG << "ConvexHull3: Computing convex hull of " << points.size() << " points";

  // Convert the input points to a packed array of floats
  TheaArray<float> coords(3 * points.size());
  array_size_t base = 0;
  for (array_size_t i = 0; i < points.size(); ++i, base += 3)
  {
    Vector3 const & p = points[i];
    coords[base    ] = (float)p.x();
    coords[base + 1] = (float)p.y();
    coords[base + 2] = (float)p.z();
  }

  approx_vertices.clear();
  approx_indices.clear();

  if (points.size() == 3)
  {
    approx_vertices.resize(3);
    approx_vertices[0] = points[0];
    approx_vertices[1] = points[1];
    approx_vertices[2] = points[2];

    approx_indices.resize(3);
    approx_indices[0] = 0;
    approx_indices[1] = 1;
    approx_indices[2] = 2;
  }
  else if (points.size() > 3)
  {
    HullDesc desc(QF_TRIANGLES,
                  (unsigned int)(points.size()),
                  &coords[0],
                  3 * sizeof(float));

    desc.SetHullFlag(QF_SKIN_WIDTH);
    desc.mNormalEpsilon  =  0.001f;
    desc.mMaxVertices    =  (unsigned int)(options.approx.max_vertices_hint > 0 ? options.approx.max_vertices_hint : 100);
    desc.mSkinWidth      =  (float)options.approx.skin_width;

    HullLibrary hull_lib;
    HullResult result;
    HullError err = hull_lib.CreateConvexHull(desc, result);
    THEA_DEBUG << "ConvexHull3: Computed convex hull with " << result.mNumOutputVertices << " vertices, " << result.mNumFaces
               << " faces, " << result.mNumIndices << " indices and error " << err;

    approx_vertices.resize((array_size_t)result.mNumOutputVertices);
    array_size_t j = 0;
    for (array_size_t i = 0; i < approx_vertices.size(); ++i, j += 3)
      approx_vertices[i] = Vector3(result.mOutputVertices[j], result.mOutputVertices[j + 1], result.mOutputVertices[j + 2]);

    if (result.mPolygons)  // shouldn't happen since we requested triangles, but still...
    {
      THEA_DEBUG << "ConvexHull3: Requested triangles, but higher-degree polygons obtained";

      Polygon3 poly;
      TheaArray<long> tri_indices;
      array_size_t num_face_vertices = 0, index;

      for (array_size_t i = 0; i < (array_size_t)result.mNumIndices; i += (1 + num_face_vertices))
      {
        num_face_vertices = (array_size_t)result.mIndices[i];

        if (num_face_vertices == 3)  // again, shouldn't happen
        {
          for (int j = 1; j <= 3; ++j)
          {
            index = (array_size_t)result.mIndices[i + j];
            debugAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
            approx_indices.push_back(index);
          }
        }
        else if (num_face_vertices > 3)  // triangulate
        {
          poly.clear();
          for (array_size_t j = 1; j <= num_face_vertices; ++j)
          {
            index = (array_size_t)result.mIndices[i + j];
            debugAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
            poly.addVertex(approx_vertices[index], (long)index);
          }

          tri_indices.clear();
          poly.triangulate(tri_indices);
          for (array_size_t j = 0; j < tri_indices.size(); ++j)
            approx_indices.push_back((array_size_t)tri_indices[j]);
        }
      }
    }
    else
    {
      approx_indices.resize((array_size_t)result.mNumIndices);
      array_size_t index;

      for (array_size_t i = 0; i < approx_indices.size(); ++i)
      {
        index = (array_size_t)result.mIndices[i];
        debugAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
        approx_indices[i] = index;
      }
    }

    hull_lib.ReleaseResult(result);
  }

  approx_updated = true;
}

} // namespace Algorithms
} // namespace Thea

// Quick-and-dirty test to see if the export-to-mesh code compiles. Comment this out for normal use.
// #include "../Graphics/GeneralMesh.hpp"
//
// void test()
// {
//   Thea::Algorithms::ConvexHull3 hull;
//   Thea::Graphics::GeneralMesh<> mesh;
//   hull.computeApprox(mesh);
// }
