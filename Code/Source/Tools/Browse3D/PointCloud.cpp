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
#include "App.hpp"
#include "Util.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Math.hpp"
#include <algorithm>
#include <fstream>

namespace Browse3D {

PointCloud::PointCloud(std::string const & path)
: has_normals(false), normals_are_normalized(false), has_graph(false)
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
  graph.clear();

  has_normals = false;
  has_graph = false;

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

  std::string feat_path = getFeaturesFilename(path);
  if (!feat_path.empty() && loadFeatures(feat_path))
    THEA_CONSOLE << getName() << ": Loaded " << features.size() << " features from '" << feat_path << '\'';

  std::string graph_path = FilePath::concat(FilePath::parent(path), FilePath::completeBaseName(path) + ".graph");
  if (FileSystem::exists(graph_path))
  {
    has_graph = true;

    std::ifstream gin(graph_path.c_str());
    if (!gin)
    {
      THEA_WARNING << getName() << ": Error opening file '" << graph_path << "' with adjacency graph of points";
      has_graph = false;
    }
    else
    {
      graph.resize(points.size());

      for (array_size_t i = 0; i < points.size(); ++i)
      {
        long num_nbrs;
        if (!(gin >> num_nbrs) || num_nbrs < 0 || num_nbrs >= (long)points.size())
        {
          THEA_WARNING << getName() << ": Error reading number of neighbors of point " << i << " from '" << graph_path << '\'';
          has_graph = false;
          break;
        }

        graph[i].resize((array_size_t)num_nbrs);

        for (array_size_t j = 0; j < graph[i].size(); ++j)
        {
          if (!(gin >> graph[i][j]) || graph[i][j] < 0 || graph[i][j] >= (long)points.size())
          {
            THEA_WARNING << getName() << ": Error reading neighbor" << j << " of point " << i << " from '" << graph_path
                         << '\'';
            has_graph = false;
            break;
          }
        }

        if (!has_graph)
          break;
      }

      if (!has_graph)
        graph.clear();

      THEA_CONSOLE << getName() << ": Loaded sample adjacency graph from '" << graph_path << '\'';
    }
  }

  return true;
}

bool
PointCloud::loadFeatures(std::string const & filename_)
{
  bool status = true;
  try
  {
    std::ifstream in(filename_.c_str());
    if (!in)
      throw Error("Could not open file");

    features.resize(1);
    features[0].resize(points.size());

    std::string line;
    Vector3 p;
    double f;

    for (array_size_t i = 0; i < points.size(); ++i)
    {
      if (!std::getline(in, line))
        throw Error(format("Could not read feature for point %ld", (long)i));

      std::istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2] >> f))
        throw Error(format("Could not read first feature for point %ld", (long)i));

      features[0][i] = (Real)f;

      if (i == 0)
      {
        while (line_in >> f)
        {
          features.push_back(TheaArray<Real>(points.size()));
          features.back()[0] = (Real)f;
        }
      }
      else
      {
        for (array_size_t j = 1; j < features.size(); ++j)
        {
          if (!(line_in >> f))
            throw Error(format("Could not read feature %ld for point %ld", (long)j, (long)i));

          features[j][i] = (Real)f;
        }
      }
    }

    if (features[0].empty())
    {
      features.clear();
      return true;
    }

    if (app().options().accentuate_features)
    {
      for (array_size_t i = 0; i < features.size(); ++i)
      {
        TheaArray<Real> sorted = features[i];
        std::sort(features.begin(), features.end());

        array_size_t tenth = (int)(0.1 * sorted.size());
        Real lo = *(sorted.begin() + tenth);

        array_size_t ninetieth = (int)(0.9 * sorted.size());
        Real hi = *(sorted.begin() + ninetieth);

        Real range = hi - lo;

        if (range < 1e-20)
        {
          lo = sorted.front();
          hi = sorted.back();
          range = hi - lo;

          if (range < 1e-20)
            continue;
        }

        if (sorted[0] >= 0)  // make a guess if this is a [0, 1] feature (e.g. SDF) or a [-1, 1] feature (e.g. curvature)
        {
          for (array_size_t j = 0; j < features[i].size(); ++j)
            features[i][j] = Math::clamp((features[i][j] - lo) / range, (Real)0, (Real)1);
        }
        else
        {
          Real abs_max = std::max(std::fabs(lo), std::fabs(hi));
          for (array_size_t j = 0; j < features[i].size(); ++j)
            features[i][j] = Math::clamp((features[i][j] + abs_max) / (2 * abs_max), (Real)0, (Real)1);
        }
      }
    }
  }
  THEA_STANDARD_CATCH_BLOCKS({ features.clear(); status = false; }, WARNING, "Couldn't load point features from '%s'",
                             filename_.c_str())

  return status;
}

std::string
PointCloud::getFeaturesFilename(std::string const & filename) const
{
  std::string ffn = filename + ".features";
  if (FileSystem::exists(ffn))
    return ffn;
  else
  {
    ffn = FilePath::concat(FilePath::parent(filename), FilePath::completeBaseName(filename) + ".features");
    if (FileSystem::exists(ffn))
      return ffn;
    else
    {
      ffn = FilePath::concat(FilePath::parent(filename), FilePath::baseName(filename) + ".features");
      if (FileSystem::exists(ffn))
        return ffn;
    }
  }

  return "";
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

  Real scale = getBounds().getExtent().length();
  Real point_radius = Math::clamp(10.0f / points.size(), 0.002f, 0.005f) * scale;
  for (array_size_t i = 0; i < points.size(); ++i)
  {
    if (!features.empty())
    {
      switch (features.size())
      {
        case 1:   render_system.setColor(ColorRGB::jetColorMap(features[0][i])); break;
        case 2:   render_system.setColor(ColorRGB(features[0][i], features[1][i], 1.0f)); break;
        default:  render_system.setColor(ColorRGB(features[0][i], features[1][i], features[2][i])); break;
      }
    }
    else if (has_normals)
    {
      Vector3 n = points[i].n;
      if (!normals_are_normalized)
        n.unitize();

      render_system.setColor(ColorRGB(0.5f * (n[0] + 1), 0.5f * (n[1] + 1), 0.5f * (n[2] + 1)));
    }

    drawSphere(render_system, points[i].p, point_radius, 8);
  }

  if (has_graph && app().options().show_graph)
  {
    render_system.pushShader();
    render_system.pushTextures();
    render_system.pushColorFlags();

      render_system.setShader(NULL);
      render_system.setTexture(0, NULL);
      render_system.setColor(ColorRGB(1, 1, 0));

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

        for (array_size_t i = 0; i < graph.size(); ++i)
        {
          for (array_size_t j = 0; j < graph[i].size(); ++j)
          {
            render_system.sendVertex(points[i].p);
            render_system.sendVertex(points[(array_size_t)graph[i][j]].p);
          }
        }

      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popTextures();
    render_system.popShader();
  }
  else if (has_normals)
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
