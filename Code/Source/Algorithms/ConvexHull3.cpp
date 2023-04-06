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

#include "ConvexHull3.hpp"
#include "../Polygon3.hpp"
#include "../ThirdParty/StanHull/hull.h"

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
  Array<float> coords(3 * points.size());
  size_t base = 0;
  for (size_t i = 0; i < points.size(); ++i, base += 3)
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

    approx_vertices.resize((size_t)result.mNumOutputVertices);
    size_t j = 0;
    for (size_t i = 0; i < approx_vertices.size(); ++i, j += 3)
      approx_vertices[i] = Vector3(result.mOutputVertices[j], result.mOutputVertices[j + 1], result.mOutputVertices[j + 2]);

    if (result.mPolygons)  // shouldn't happen since we requested triangles, but still...
    {
      THEA_DEBUG << "ConvexHull3: Requested triangles, but higher-degree polygons obtained";

      Polygon3 poly;
      Array<intx> tri_indices;
      size_t num_face_vertices = 0, index;

      for (size_t i = 0; i < (size_t)result.mNumIndices; i += (1 + num_face_vertices))
      {
        num_face_vertices = (size_t)result.mIndices[i];

        if (num_face_vertices == 3)  // again, shouldn't happen
        {
          for (int j = 1; j <= 3; ++j)
          {
            index = (size_t)result.mIndices[i + j];
            theaAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
            approx_indices.push_back(index);
          }
        }
        else if (num_face_vertices > 3)  // triangulate
        {
          poly.clear();
          for (size_t j = 1; j <= num_face_vertices; ++j)
          {
            index = (size_t)result.mIndices[i + j];
            theaAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
            poly.addVertex(approx_vertices[index], (intx)index);
          }

          tri_indices.clear();
          poly.triangulate(tri_indices);
          for (size_t j = 0; j < tri_indices.size(); ++j)
            approx_indices.push_back((size_t)tri_indices[j]);
        }
      }
    }
    else
    {
      approx_indices.resize((size_t)result.mNumIndices);
      size_t index;

      for (size_t i = 0; i < approx_indices.size(); ++i)
      {
        index = (size_t)result.mIndices[i];
        theaAssertM(index >= 0 && index < approx_vertices.size(), "ConvexHull3: Vertex index out of bounds");
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
