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
// First version: 2012
//
//============================================================================

#include "PointCloud.hpp"
#include "App.hpp"
#include "Util.hpp"
#include "../../Graphics/IBufferPool.hpp"
#include "../../Graphics/IBuffer.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Math.hpp"
#include <algorithm>
#include <fstream>

namespace Browse3D {

PointCloud::PointCloud(std::string const & path, std::string const & features_path)
: has_normals(false), normals_are_normalized(false), has_graph(false), changed_buffers(0), buf_pool(nullptr),
  vertices_buf(nullptr), colors_buf(nullptr)
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
      n = Vector3::Zero();

    points.push_back(Point(p, n));
  }

  if (has_normals)
  {
    normals_are_normalized = true;
    for (size_t i = 0; i < points.size(); ++i)
      if (!Math::fuzzyEq(points[i].n.squaredNorm(), (Real)1, (Real)0.001))
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

      intx max_nbrs;
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
          intx num_nbrs;
          if (!(line_in >> num_nbrs) || num_nbrs < 0 || num_nbrs > max_nbrs || num_nbrs >= (intx)points.size())
          {
            THEA_WARNING << getName() << ": Error reading valid number of neighbors of point " << i << " from '" << graph_path
                         << '\'';
            has_graph = false;
            break;
          }

          graph[i].resize((size_t)num_nbrs);

          for (size_t j = 0; j < graph[i].size(); ++j)
          {
            if (!(line_in >> graph[i][j]) || graph[i][j] < 0 || graph[i][j] >= (intx)points.size())
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
readFeaturesTXT(std::string const & path, intx num_points, Array< Array<Real> > & features, bool has_point_prefix)
{
  std::ifstream in(path.c_str());
  if (!in)
    throw Error("TXT: Could not open file '" + path + '\'');

  features.resize(1);
  features[0].resize((size_t)num_points);

  std::string line;
  Vector3 p;
  double f;

  for (intx i = 0; i < num_points; ++i)
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
        features.push_back(Array<Real>((intx)num_points));
        features.back()[0] = (Real)f;
      }
    }
    else
    {
      for (size_t j = 1; j < features.size(); ++j)
      {
        if (!(line_in >> f))
          throw Error(format("TXT: Could not read feature %ld for point %ld", (intx)j, i) + " from '" + path + '\'');

        features[j][(size_t)i] = (Real)f;
      }
    }
  }

  return true;
}

bool
readFeaturesARFF(std::string const & path, intx num_points, Array< Array<Real> > & features)
{
  std::ifstream in(path.c_str());
  if (!in)
    throw Error("ARFF: Could not open file '" + path + '\'');

  features.clear();

  std::string line;
  intx num_features = 0;
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

  Array<std::string> fields;
  double f;
  for (intx i = 0; i < num_points; ++i)
  {
    do
    {
      if (!std::getline(in, line))
        throw Error(format("ARFF: Could not read features for point %ld", i) + " from '" + path + '\'');

      line = trimWhitespace(line);
    } while (line.empty());

    stringSplit(line, ',', fields);
    if ((intx)fields.size() < num_features)
    {
      THEA_ERROR << "ARFF: Point " << i << " does not have enough features (" << path << ')';
      return false;
    }

    for (size_t j = 0; j < features.size(); ++j)
    {
      std::istringstream field_in(fields[j]);
      if (!(field_in >> f))
        throw Error(format("ARFF: Could not read feature %ld for point %ld", (intx)j, i) + " from '" + path + '\'');

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
      status = PointCloudInternal::readFeaturesARFF(filename_, (intx)points.size(), features);
    else if (endsWith(filename_lc, ".features") || endsWith(filename_lc, ".feat"))
      status = PointCloudInternal::readFeaturesTXT(filename_, (intx)points.size(), features, true);
    else
      status = PointCloudInternal::readFeaturesTXT(filename_, (intx)points.size(), features, false);

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
          Array<Real> sorted = features[i];
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
PointCloud::setPointColors(Array<ColorRgba> const & colors_)
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

ColorRgba
PointCloud::getColor(size_t point_index) const
{
  alwaysAssertM(point_index >= 0 && point_index < points.size(),
                format("%s: Index %ld out of bounds", getName(), (intx)point_index));

  if (!colors.empty())
  {
    return colors[point_index];
  }
  else if (!features.empty())
  {
    switch (features.size())
    {
      case 1:   return ColorRgb::jetColorMap(0.2 + 0.6 * features[0][point_index]);
      case 2:   return ColorRgb(features[0][point_index], features[1][point_index], 1.0f);
      default:  return ColorRgb(features[0][point_index], features[1][point_index], features[2][point_index]);
    }
  }
  else if (has_normals)
  {
    Vector3 n = points[point_index].n;
    if (!normals_are_normalized)
      n.normalize();

    return ColorRgb(0.5f * (n[0] + 1), 0.5f * (n[1] + 1), 0.5f * (n[2] + 1));
  }
  else if (app().options().fancy_colors)
  {
    static Real const MIN_COLOR = 0.1;
    static Real const MAX_COLOR = 1.0;

    Vector3 ext = getBounds().getExtent().cwiseMax(Vector3(1e-20f, 1e-20f, 1e-20f));
    Vector3 v = (points[point_index].p - getBounds().getLow()).cwiseQuotient(ext);
    ColorRgb c((MAX_COLOR - MIN_COLOR) * v + Vector3(MIN_COLOR, MIN_COLOR, MIN_COLOR));
    Vector3 hsv; hsv << c.toHSV().head<2>(), 1;
    return ColorRgb::fromHSV(hsv);
  }

  return ColorRgb::zero();
}

bool
PointCloud::uploadToGraphicsSystem(Graphics::IRenderSystem & render_system)
{
  if (app().options().fancy_points) return true;
  if (changed_buffers == 0) return true;

  vertices_buf = colors_buf = nullptr;

  if (points.empty())
  {
    if (buf_pool)
    {
      render_system.destroyBufferPool(buf_pool);
      buf_pool = nullptr;
    }

    changed_buffers = 0;
    return true;
  }

  bool has_colors = (!features.empty() || has_normals || app().options().fancy_colors);

  static int const PADDING = 32;
  intx vertex_bytes  =  3 * 4 * (intx)points.size() + PADDING;  // 3 * float
  intx color_bytes   =  has_colors ?  4 * 4 * (intx)points.size() + PADDING : 0;  // 4 * float

  intx num_bytes = vertex_bytes + color_bytes + PADDING;

  if (buf_pool)
  {
    if (buf_pool->getCapacity() <= num_bytes || buf_pool->getCapacity() > (intx)(1.5 * num_bytes))
    {
      render_system.destroyBufferPool(buf_pool);

      std::string bufpool_name = getNameStr() + " buffer pool";
      buf_pool = render_system.createBufferPool(bufpool_name.c_str(), num_bytes,
                                                Graphics::IBufferPool::Usage::WRITE_OCCASIONALLY, true);
      if (!buf_pool) { THEA_ERROR << getName() << ": Couldn't create buffer pool"; return false; }
    }
    else
      buf_pool->reset();
  }
  else
  {
    std::string bufpool_name = getNameStr() + " buffer pool";
    buf_pool = render_system.createBufferPool(bufpool_name.c_str(), num_bytes, Graphics::IBufferPool::Usage::WRITE_OCCASIONALLY,
                                              true);
    if (!buf_pool) { THEA_ERROR << getName() << ": Couldn't create buffer pool"; return false; }
  }

  if (!points.empty())
  {
    vertices_buf = buf_pool->createBuffer(vertex_bytes);
    if (!vertices_buf) { THEA_ERROR << getName() << ": Couldn't create vertex buffer"; return false; }

    Array<Vector3> vbuf(points.size());
    for (size_t i = 0; i < points.size(); ++i)
      vbuf[i] = points[i].p;

    if (!vertices_buf->updateAttributes(0, (intx)vbuf.size(), 3, NumericType::REAL, &vbuf[0])) return false;
  }

  if (has_colors)
  {
    colors_buf = buf_pool->createBuffer(color_bytes);
    if (!colors_buf) { THEA_ERROR << getName() << ": Couldn't create color buffer"; return false; }

    Array<ColorRgba> cbuf(points.size());
    for (size_t i = 0; i < points.size(); ++i)
      cbuf[i] = getColor(i);

    if (!colors_buf->updateAttributes(0, (intx)cbuf.size(), 4, NumericType::REAL, &cbuf[0])) return false;
  }

  changed_buffers = 0;

  return true;
}

int8
PointCloud::draw(Graphics::IRenderSystem * render_system, Graphics::IRenderOptions const * options) const
{
  if (isEmpty())
    return true;

  if (!options) options = Graphics::RenderOptions::defaults();

  if (!const_cast<PointCloud *>(this)->uploadToGraphicsSystem(*render_system)) return false;

  render_system->pushShader();
  render_system->pushTextures();
  render_system->pushColorFlags();

    render_system->setTexture(0, nullptr);

    if (app().options().fancy_points)
    {
      setPhongShader(*render_system);

      bool has_colors = (!features.empty()
                       || has_normals
                       || app().options().fancy_colors);
      Real scale = getBounds().getExtent().norm();
      Real point_radius = app().options().point_scale * Math::clamp(10.0f / points.size(), 0.002f, 0.005f) * scale;
      for (size_t i = 0; i < points.size(); ++i)
      {
        if (has_colors)
          render_system->setColor(getColor(i).data());

        drawSphere(*render_system, points[i].p, point_radius, 8);
      }
    }
    else
    {
      render_system->setShader(nullptr);

      render_system->beginIndexedPrimitives();

        render_system->setVertexBuffer(vertices_buf);
        if (colors_buf) render_system->setColorBuffer(colors_buf);

        render_system->pushShapeFlags();
          render_system->setPointSmooth(true);
          render_system->setPointSize(2 * app().options().point_scale);

          render_system->sendSequentialIndices(Graphics::IRenderSystem::Primitive::POINTS, 0, (intx)points.size());

        render_system->popShapeFlags();

      render_system->endIndexedPrimitives();
    }

    if (has_graph && app().options().show_graph)
    {
      render_system->setShader(nullptr);
      render_system->setColor(ColorRgb(1, 1, 0).data());

      render_system->beginPrimitive(Graphics::IRenderSystem::Primitive::LINES);

        for (size_t i = 0; i < graph.size(); ++i)
        {
          for (size_t j = 0; j < graph[i].size(); ++j)
          {
            render_system->sendVertex(3, points[i].p.data());
            render_system->sendVertex(3, points[(size_t)graph[i][j]].p.data());
          }
        }

      render_system->endPrimitive();
    }
    else if (has_normals && app().options().show_normals)
    {
      render_system->setShader(nullptr);
      render_system->setColor(ColorRgb(0, 0, 1).data());

      Real normal_scale = (normals_are_normalized ? 0.025f * getBounds().getExtent().norm() : 1);

      render_system->beginPrimitive(Graphics::IRenderSystem::Primitive::LINES);

        for (size_t i = 0; i < points.size(); ++i)
        {
          render_system->sendVertex(3, points[i].p.data());
          render_system->sendVertex(3, Vector3(points[i].p + normal_scale * points[i].n).data());
        }

      render_system->endPrimitive();
    }

  render_system->popColorFlags();
  render_system->popTextures();
  render_system->popShader();

  char const * err = nullptr;
  if ((err = render_system->getAndClearError()))
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

} // namespace Browse3D
