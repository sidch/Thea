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
#include "../../Graphics/VARArea.hpp"
#include "../../Graphics/VAR.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Math.hpp"
#include <algorithm>
#include <fstream>

namespace Browse3D {

PointCloud::PointCloud(std::string const & path, std::string const & features_path)
: has_normals(false), normals_are_normalized(false), has_graph(false), changed_buffers(0), var_area(NULL), vertices_var(NULL),
  colors_var(NULL)
{
  if (!path.empty())
    load(path, features_path);
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

  colors.clear();
  features.clear();

  has_graph = false;
  graph.clear();

  changed_buffers = 0xFFFF;
}

bool
PointCloud::load(std::string const & path, std::string const & features_path)
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
    for (size_t i = 0; i < points.size(); ++i)
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

  std::string features_file_path = getDefaultFeaturesFilename(path, features_path);
  if (!features_file_path.empty() && loadFeatures(features_file_path))
    THEA_CONSOLE << getName() << ": Loaded " << features.size() << " feature(s) from '" << features_file_path << '\'';

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

      std::string line;
      if (!std::getline(gin, line))
      {
        THEA_WARNING << getName() << ": Error reading maximum number of neighbors in graph '" << graph_path << '\'';
        has_graph = false;
      }

      long max_nbrs;
      if (has_graph)
      {
        std::istringstream line_in(line);
        if (!(line_in >> max_nbrs) || max_nbrs < 0)
        {
          THEA_WARNING << getName() << ": Error reading maximum number of neighbors in graph '" << graph_path << '\'';
          has_graph = false;
        }
      }

      if (has_graph)
      {
        for (size_t i = 0; i < points.size(); ++i)
        {
          if (!std::getline(gin, line))
          {
            THEA_WARNING << getName() << ": Error reading neighbors of point " << i << " in graph '" << graph_path << '\'';
            has_graph = false;
            break;
          }

          std::istringstream line_in(line);
          long num_nbrs;
          if (!(line_in >> num_nbrs) || num_nbrs < 0 || num_nbrs > max_nbrs || num_nbrs >= (long)points.size())
          {
            THEA_WARNING << getName() << ": Error reading valid number of neighbors of point " << i << " from '" << graph_path
                         << '\'';
            has_graph = false;
            break;
          }

          graph[i].resize((size_t)num_nbrs);

          for (size_t j = 0; j < graph[i].size(); ++j)
          {
            if (!(line_in >> graph[i][j]) || graph[i][j] < 0 || graph[i][j] >= (long)points.size())
            {
              THEA_WARNING << getName() << ": Error reading valid neighbor" << j << " of point " << i << " from '" << graph_path
                           << '\'';
              has_graph = false;
              break;
            }
          }

          if (!has_graph)
            break;
        }
      }

      if (!has_graph)
        graph.clear();

      THEA_CONSOLE << getName() << ": Loaded sample adjacency graph from '" << graph_path << '\'';
    }
  }

  changed_buffers = 0xFFFF;
  return true;
}

namespace PointCloudInternal {

bool
readFeaturesTXT(std::string const & path, long num_points, TheaArray< TheaArray<Real> > & features, bool has_point_prefix)
{
  std::ifstream in(path.c_str());
  if (!in)
    throw Error("TXT: Could not open file '" + path + '\'');

  features.resize(1);
  features[0].resize((size_t)num_points);

  std::string line;
  Vector3 p;
  double f;

  for (long i = 0; i < num_points; ++i)
  {
    if (!std::getline(in, line))
      throw Error(format("TXT: Could not read feature for point %ld", i) + " from '" + path + '\'');

    std::istringstream line_in(line);
    if (has_point_prefix) line_in >> p[0] >> p[1] >> p[2];
    if (!(line_in >> f))
      throw Error(format("TXT: Could not read first feature for point %ld", i) + " from '" + path + '\'');

    features[0][(size_t)i] = (Real)f;

    if (i == 0)
    {
      while (line_in >> f)
      {
        features.push_back(TheaArray<Real>((long)num_points));
        features.back()[0] = (Real)f;
      }
    }
    else
    {
      for (size_t j = 1; j < features.size(); ++j)
      {
        if (!(line_in >> f))
          throw Error(format("TXT: Could not read feature %ld for point %ld", (long)j, i) + " from '" + path + '\'');

        features[j][(size_t)i] = (Real)f;
      }
    }
  }

  return true;
}

bool
readFeaturesARFF(std::string const & path, long num_points, TheaArray< TheaArray<Real> > & features)
{
  std::ifstream in(path.c_str());
  if (!in)
    throw Error("ARFF: Could not open file '" + path + '\'');

  features.clear();

  std::string line;
  long num_features = 0;
  while (std::getline(in, line))
  {
    line = trimWhitespace(line);
    if (line.empty())
      continue;

    line = toLower(line);
    if (beginsWith(line, "$data"))
      break;
    else if (beginsWith(line, "@attribute"))
    {
      std::string field;
      std::istringstream line_in(line);
      line_in >> field; line_in >> field;
      if (field == "class")
        break;
      else if (!endsWith(line, "numeric"))
      {
        THEA_ERROR << "ARFF: Non-numeric attribute not supported (" << path << ')';
        return false;
      }
      else
        num_features++;
    }
  }

  if (num_features <= 0)
  {
    THEA_WARNING << "ARFF: No features found in '" << path << '\'';
    return true;
  }

  features.resize((size_t)num_features);
  for (size_t i = 0; i < features.size(); ++i)
    features[i].resize((size_t)num_points);

  TheaArray<std::string> fields;
  double f;
  for (long i = 0; i < num_points; ++i)
  {
    do
    {
      if (!std::getline(in, line))
        throw Error(format("ARFF: Could not read features for point %ld", i) + " from '" + path + '\'');

      line = trimWhitespace(line);
    } while (line.empty());

    stringSplit(line, ',', fields);
    if ((long)fields.size() < num_features)
    {
      THEA_ERROR << "ARFF: Point " << i << " does not have enough features (" << path << ')';
      return false;
    }

    for (size_t j = 0; j < features.size(); ++j)
    {
      std::istringstream field_in(fields[j]);
      if (!(field_in >> f))
        throw Error(format("ARFF: Could not read feature %ld for point %ld", (long)j, i) + " from '" + path + '\'');

      features[j][(size_t)i] = (Real)f;
    }
  }

  return true;
}

} // namespace PointCloudInternal

bool
PointCloud::loadFeatures(std::string const & filename_)
{
  bool status = true;
  try
  {
    std::string filename_lc = toLower(filename_);
    bool status = true;
    if (endsWith(filename_lc, ".arff"))
      status = PointCloudInternal::readFeaturesARFF(filename_, (long)points.size(), features);
    else if (endsWith(filename_lc, ".features") || endsWith(filename_lc, ".feat"))
      status = PointCloudInternal::readFeaturesTXT(filename_, (long)points.size(), features, true);
    else
      status = PointCloudInternal::readFeaturesTXT(filename_, (long)points.size(), features, false);

    if (!status)
      return false;

    if (features[0].empty())
    {
      features.clear();
      return true;
    }

    if (app().options().accentuate_features)
    {
      if (app().options().color_cube_features && features.size() == 3)
      {
        Real abs_max = -1;
        for (size_t i = 0; i < features.size(); ++i)
        {
          for (size_t j = 0; j < features[i].size(); ++j)
          {
            Real abs_feat_val = std::fabs(features[i][j]);
            if (abs_feat_val > abs_max)
              abs_max = abs_feat_val;
          }
        }

        if (abs_max > 0)
        {
          for (size_t i = 0; i < features.size(); ++i)
            for (size_t j = 0; j < features[i].size(); ++j)
              features[i][j] = Math::clamp(0.5 * (features[i][j]  / abs_max + 1), (Real)0, (Real)1);
        }
      }
      else
      {
        for (size_t i = 0; i < features.size(); ++i)
        {
          TheaArray<Real> sorted = features[i];
          std::sort(sorted.begin(), sorted.end());

          size_t tenth = (int)(0.1 * sorted.size());
          Real lo = *(sorted.begin() + tenth);

          size_t ninetieth = (int)(0.9 * sorted.size());
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
            for (size_t j = 0; j < features[i].size(); ++j)
              features[i][j] = Math::clamp((features[i][j] - lo) / range, (Real)0, (Real)1);
          }
          else
          {
            Real abs_max = std::max(std::fabs(lo), std::fabs(hi));
            for (size_t j = 0; j < features[i].size(); ++j)
              features[i][j] = Math::clamp((features[i][j] + abs_max) / (2 * abs_max), (Real)0, (Real)1);
          }
        }
      }
    }
  }
  THEA_STANDARD_CATCH_BLOCKS({ features.clear(); status = false; }, WARNING, "Couldn't load point features from '%s'",
                             filename_.c_str())

  return status;
}

std::string
PointCloud::getDefaultFeaturesFilename(std::string const & filename, std::string const & features_path) const
{
  if (FileSystem::fileExists(features_path))
    return features_path;

  std::string app_features_path = app().options().features;
  if (FileSystem::fileExists(app_features_path))
    return app_features_path;

  static std::string const EXTS[] = { ".arff", ".features", ".feat" };  // in order of decreasing priority
  static size_t NUM_EXTS = sizeof(EXTS) / sizeof(std::string);
  static int const NUM_DIRS = 3;

  for (int i = 0; i < NUM_DIRS; ++i)
  {
    std::string dir;
    switch (i)
    {
      case 0: dir = features_path; break;
      case 1: dir = app_features_path; break;
      default: dir = FilePath::parent(filename);
    }

    if (i < NUM_DIRS - 1 && !FileSystem::directoryExists(dir))
      continue;

    for (size_t j = 0; j < NUM_EXTS; ++j)
    {
      std::string ffn = FilePath::concat(dir, filename + EXTS[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }

    for (size_t j = 0; j < NUM_EXTS; ++j)
    {
      std::string ffn = FilePath::concat(dir, FilePath::completeBaseName(filename) + EXTS[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }

    for (size_t j = 0; j < NUM_EXTS; ++j)
    {
      std::string ffn = FilePath::concat(dir, FilePath::baseName(filename) + EXTS[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }
  }

  return "";
}

bool
PointCloud::setPointColors(TheaArray<ColorRGBA> const & colors_)
{
  if (colors_.size() < points.size())
  {
    THEA_ERROR << getName() << ": Number of colors < number of points";
    return false;
  }

  colors = colors_;

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
    for (size_t i = 0; i < points.size(); ++i)
      bounds.merge(points[i].p);
  }
}

ColorRGBA
PointCloud::getColor(size_t point_index) const
{
  alwaysAssertM(point_index >= 0 && point_index < points.size(),
                format("%s: Index %ld out of bounds", getName(), (long)point_index));

  if (!colors.empty())
  {
    return colors[point_index];
  }
  else if (!features.empty())
  {
    switch (features.size())
    {
      case 1:   return ColorRGB::jetColorMap(0.2 + 0.6 * features[0][point_index]);
      case 2:   return ColorRGB(features[0][point_index], features[1][point_index], 1.0f);
      default:  return ColorRGB(features[0][point_index], features[1][point_index], features[2][point_index]);
    }
  }
  else if (has_normals)
  {
    Vector3 n = points[point_index].n;
    if (!normals_are_normalized)
      n.unitize();

    return ColorRGB(0.5f * (n[0] + 1), 0.5f * (n[1] + 1), 0.5f * (n[2] + 1));
  }
  else if (app().options().fancy_colors)
  {
    static Real const MIN_COLOR = 0.1;
    static Real const MAX_COLOR = 1.0;

    Vector3 ext = getBounds().getExtent().max(Vector3(1e-20f, 1e-20f, 1e-20f));
    Vector3 v = (points[point_index].p - getBounds().getLow()) / ext;
    ColorRGB c((MAX_COLOR - MIN_COLOR) * v + Vector3(MIN_COLOR, MIN_COLOR, MIN_COLOR));
    Vector3 hsv(c.toHSV().xy(), 1);
    return ColorRGB::fromHSV(hsv);
  }

  return ColorRGB::zero();
}

void
PointCloud::uploadToGraphicsSystem(Graphics::RenderSystem & render_system)
{
  if (app().options().fancy_points) return;
  if (changed_buffers == 0) return;

  vertices_var = colors_var = NULL;

  if (points.empty())
  {
    if (var_area)
    {
      render_system.destroyVARArea(var_area);
      var_area = NULL;
    }

    changed_buffers = 0;
    return;
  }

  bool has_colors = (!features.empty() || has_normals || app().options().fancy_colors);

  static int const PADDING = 32;
  long vertex_bytes  =  3 * 4 * (long)points.size() + PADDING;  // 3 * float
  long color_bytes   =  has_colors ?  4 * 4 * (long)points.size() + PADDING : 0;  // 4 * float

  long num_bytes = vertex_bytes + color_bytes + PADDING;

  if (var_area)
  {
    if (var_area->getCapacity() <= num_bytes || var_area->getCapacity() > (long)(1.5 * num_bytes))
    {
      render_system.destroyVARArea(var_area);

      std::string vararea_name = getNameStr() + " VAR area";
      var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, Graphics::VARArea::Usage::WRITE_OCCASIONALLY,
                                             true);
      if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
    }
    else
      var_area->reset();
  }
  else
  {
    std::string vararea_name = getNameStr() + " VAR area";
    var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, Graphics::VARArea::Usage::WRITE_OCCASIONALLY, true);
    if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
  }

  if (!points.empty())
  {
    vertices_var = var_area->createArray(vertex_bytes);
    if (!vertices_var) throw Error(getNameStr() + ": Couldn't create vertices VAR");

    TheaArray<Vector3> vbuf(points.size());
    for (size_t i = 0; i < points.size(); ++i)
      vbuf[i] = points[i].p;

    vertices_var->updateVectors(0, (long)vbuf.size(), &vbuf[0]);
  }

  if (has_colors)
  {
    colors_var = var_area->createArray(color_bytes);
    if (!colors_var) throw Error(getNameStr() + ": Couldn't create colors VAR");

    TheaArray<ColorRGBA> cbuf(points.size());
    for (size_t i = 0; i < points.size(); ++i)
      cbuf[i] = getColor(i);

    colors_var->updateColors(0, (long)cbuf.size(), &cbuf[0]);
  }

  changed_buffers = 0;
}

void
PointCloud::draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
{
  if (isEmpty())
    return;

  const_cast<PointCloud *>(this)->uploadToGraphicsSystem(render_system);

  render_system.pushShader();
  render_system.pushTextures();
  render_system.pushColorFlags();

    render_system.setTexture(0, NULL);

    if (app().options().fancy_points)
    {
      setPhongShader(render_system);

      bool has_colors = (!features.empty()
                       || has_normals
                       || app().options().fancy_colors);
      Real scale = getBounds().getExtent().length();
      Real point_radius = app().options().point_scale * Math::clamp(10.0f / points.size(), 0.002f, 0.005f) * scale;
      for (size_t i = 0; i < points.size(); ++i)
      {
        if (has_colors)
          render_system.setColor(getColor(i));

        drawSphere(render_system, points[i].p, point_radius, 8);
      }
    }
    else
    {
      render_system.setShader(NULL);

      render_system.beginIndexedPrimitives();

        render_system.setVertexArray(vertices_var);
        if (colors_var) render_system.setColorArray(colors_var);

        render_system.pushShapeFlags();
        render_system.setPointSize(2 * app().options().point_scale);

          render_system.sendSequentialIndices(Graphics::RenderSystem::Primitive::POINTS, 0, (long)points.size());

        render_system.popShapeFlags();

      render_system.endIndexedPrimitives();
    }

    if (has_graph && app().options().show_graph)
    {
      render_system.setShader(NULL);
      render_system.setColor(ColorRGB(1, 1, 0));

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

        for (size_t i = 0; i < graph.size(); ++i)
        {
          for (size_t j = 0; j < graph[i].size(); ++j)
          {
            render_system.sendVertex(points[i].p);
            render_system.sendVertex(points[(size_t)graph[i][j]].p);
          }
        }

      render_system.endPrimitive();
    }
    else if (has_normals && app().options().show_normals)
    {
      render_system.setShader(NULL);
      render_system.setColor(ColorRGB(0, 0, 1));

      Real normal_scale = (normals_are_normalized ? 0.025f * getBounds().getExtent().length() : 1);

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

        for (size_t i = 0; i < points.size(); ++i)
        {
          render_system.sendVertex(points[i].p);
          render_system.sendVertex(points[i].p + normal_scale * points[i].n);
        }

      render_system.endPrimitive();
    }

  render_system.popColorFlags();
  render_system.popTextures();
  render_system.popShader();
}

} // namespace Browse3D
