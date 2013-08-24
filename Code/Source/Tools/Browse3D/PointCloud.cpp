//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#include "PointCloud.hpp"
#include "Util.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../Math.hpp"
#include <fstream>

namespace Browse3D {

PointCloud::PointCloud(std::string const & path)
: has_normals(false), normals_are_normalized(false)
{
  if (!path.empty())
    load(path);
}

PointCloud::~PointCloud()
{
}

void
PointCloud::clear()
{
  points.clear();
  has_normals = false;
  bounds.setNull();
}

bool
PointCloud::load(std::string const & path)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Couldn't load points from file '" << path << '\'';
    return false;
  }

  clear();
  has_normals = false;

  std::string line;
  Vector3 p, n;
  while (getline(in, line))
  {
    std::istringstream line_in(line);
    if (!(line_in >> p[0] >> p[1] >> p[2]))
      continue;

    if (line_in >> n[0] >> n[1] >> n[2])
      has_normals = true;
    else
      n = Vector3::zero();

    points.push_back(Point(p, n));
  }

  if (has_normals)
  {
    normals_are_normalized = true;
    for (array_size_t i = 0; i < points.size(); ++i)
      if (!Math::fuzzyEq(points[i].n.squaredLength(), (Real)1, (Real)0.001))
      {
        normals_are_normalized = false;
        break;
      }
  }

  updateBounds();

  setName(FilePath::baseName(path));

  THEA_CONSOLE << getName() << ": Loaded " << points.size() << " points with bounding box " << bounds.toString() << " from '"
               << path << '\'';

  return true;
}

AxisAlignedBox3 const &
PointCloud::getBounds() const
{
  return bounds;
}

void
PointCloud::invalidateBounds()
{
  bounds.setNull();
}

void
PointCloud::updateBounds()
{
  if (bounds.isNull())
  {
    for (array_size_t i = 0; i < points.size(); ++i)
      bounds.merge(points[i].p);
  }
}

void
PointCloud::uploadToGraphicsSystem(Graphics::RenderSystem & render_system)
{
}

void
PointCloud::draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
{
  if (isEmpty())
    return;

  const_cast<PointCloud *>(this)->uploadToGraphicsSystem(render_system);

  Real point_radius = 0.005f * getBounds().getExtent().length();
  for (array_size_t i = 0; i < points.size(); ++i)
  {
    if (has_normals)
    {
      Vector3 n = points[i].n;
      if (!normals_are_normalized)
        n.unitize();

      render_system.setColor(ColorRGB(0.5f * (n[0] + 1), 0.5f * (n[1] + 1), 0.5f * (n[2] + 1)));
    }

    drawSphere(render_system, points[i].p, point_radius);
  }

  if (has_normals)
  {
    render_system.pushShader();
    render_system.pushTextures();
    render_system.pushColorFlags();

      render_system.setShader(NULL);
      render_system.setTexture(0, NULL);
      render_system.setColor(ColorRGB(0, 0, 1));

      Real normal_scale = (normals_are_normalized ? 0.025f * getBounds().getExtent().length() : 1);

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

        for (array_size_t i = 0; i < points.size(); ++i)
        {
          render_system.sendVertex(points[i].p);
          render_system.sendVertex(points[i].p + normal_scale * points[i].n);
        }

      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popTextures();
    render_system.popShader();
  }
}

} // namespace Browse3D
