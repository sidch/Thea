//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#include "Model.hpp"
#include "App.hpp"
#include "MainWindow.hpp"
#include "Math.hpp"
#include "Mesh.hpp"
#include "Util.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/RayIntersectionTester.hpp"
#include "../../Graphics/MeshCodec.hpp"
#include <QFileDialog>
#include <QFileInfo>
#include <QMouseEvent>
#include <fstream>

namespace Browse3D {

namespace ModelInternal {

bool first_file_dialog = true;

QString
getWorkingDir()
{
  if (first_file_dialog)
    if (!app().options().working_dir.isEmpty() && QFileInfo(app().options().working_dir).isDir())
      return app().options().working_dir;

  return QString();
}

bool averageNormals(Mesh & mesh)
{
  if (!mesh.hasNormals())
    mesh.computeAveragedVertexNormals();

  return false;
}

bool enableWireframe(Mesh & mesh)
{
  mesh.setWireframeEnabled(true);
  return false;
}

} // namespace ModelInternal

Model::Model(QString const & initial_mesh)
: valid_pick(false),
  selected_sample(-1),
  valid_kdtree(true),
  kdtree(new KDTree),
  valid_vertex_kdtree(true),
  vertex_kdtree(new VertexKDTree)
{
  load(initial_mesh);

  picked_sample.type = "Picked";
}

Model::~Model()
{
  delete vertex_kdtree;
  delete kdtree;
}

QString
Model::getName() const
{
  return toQString(mesh_group->getName());
}

bool
Model::isEmpty() const
{
  return (!mesh_group || mesh_group->isEmpty()) && points.empty();
}

void
Model::clear()
{
  clearMesh();
  clearPoints();
  invalidateAll();
}

void
Model::clearMesh()
{
  if (mesh_group) mesh_group->clear();
  samples.clear();
}

void
Model::clearPoints()
{
  points.clear();
  point_bounds.setNull();
}

bool
Model::load(QString const & filename_)
{
  if (filename_.isEmpty())
    return false;

  QFileInfo info(filename_);
  if (!info.exists() || info.canonicalFilePath() == QFileInfo(filename).canonicalFilePath())
    return false;

  if (filename_.toLower().endsWith(".pts"))
  {
    std::ifstream in(toStdString(filename_).c_str());
    if (!in)
    {
      THEA_ERROR << "Couldn't load points from file '" << filename_ << '\'';
      return false;
    }

    clear();
    filename = filename_;

    std::string line;
    Vector3 p;
    while (getline(in, line))
    {
      std::istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2]))
        continue;

      points.push_back(p);
      point_bounds.merge(p);
    }

    bounds = point_bounds;

    THEA_CONSOLE << "Loaded " << points.size() << " points with bounding box " << point_bounds.toString() << " from '"
                 << filename_ << '\'';
  }
  else
  {
    MeshGroupPtr new_mesh_group(new MeshGroup("Mesh Group"));
    Mesh::nextFaceIndex(Mesh::Face(), true);  // reset counting

    static CodecOBJ<Mesh> const obj_codec(NULL, CodecOBJ<Mesh>::ReadOptions().setIgnoreTexCoords(true));
    try
    {
      if (filename_.endsWith("obj", Qt::CaseInsensitive))
        new_mesh_group->load(toStdString(filename_), obj_codec);
      else
        new_mesh_group->load(toStdString(filename_));
    }
    THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "Couldn't load model '%s'", toStdString(filename_).c_str())

    invalidateAll();

    new_mesh_group->forEachMeshUntil(ModelInternal::averageNormals);
    new_mesh_group->forEachMeshUntil(ModelInternal::enableWireframe);

    mesh_group = new_mesh_group;
    clearPoints();
    filename = filename_;

    bounds = new_mesh_group->getBounds();

    THEA_CONSOLE << "Loaded model '" << filename_ << " with bounding box " << mesh_group->getBounds().toString();

    loadSamples(getSamplesFilename());
  }

  emit filenameChanged(filename);
  emit geometryChanged(this);

  return true;
}

bool
Model::selectAndLoad()
{
  QString filename_ = QFileDialog::getOpenFileName(app().getMainWindow(), tr("Load model"), ModelInternal::getWorkingDir(),
                                                   tr("Model files (*.3ds *.obj *.off *.off.bin *.pts)"));

  bool success = load(filename_);
  if (success)
    ModelInternal::first_file_dialog = false;

  return success;
}

void
Model::invalidateAll()
{
  invalidateVertexKDTree();
  invalidateKDTree();
}

void
Model::invalidateKDTree()
{
  valid_kdtree = false;
}

void
Model::updateKDTree() const
{
  if (valid_kdtree) return;

  kdtree->clear(false);

  if (mesh_group)
  {
    kdtree->add(*mesh_group);
    kdtree->init();

    THEA_CONSOLE << getName() << ": Updated kd-tree";
  }

  valid_kdtree = true;
}

void
Model::invalidateVertexKDTree()
{
  valid_vertex_kdtree = false;
}

namespace ModelInternal {

struct CollectVerticesFunctor
{
  CollectVerticesFunctor(TheaArray<IndexedMeshVertex> * verts_) : verts(verts_) {}

  bool operator()(Mesh & mesh)
  {
    Mesh::VertexArray const & mv = mesh.getVertices();
    for (array_size_t i = 0; i < mv.size(); ++i)
      verts->push_back(IndexedMeshVertex(&mesh, (long)i, mv[i]));

    return false;
  }

  TheaArray<IndexedMeshVertex> * verts;

}; // struct CollectVerticesFunctor

} // namespace ModelInternal

void
Model::updateVertexKDTree() const
{
  if (valid_vertex_kdtree) return;

  vertex_kdtree->clear(false);

  TheaArray<IndexedMeshVertex> verts;
  ModelInternal::CollectVerticesFunctor func(&verts);
  mesh_group->forEachMeshUntil(&func);
  vertex_kdtree->init(verts.begin(), verts.end());

  valid_vertex_kdtree = true;
}

Model::KDTree const &
Model::getKDTree(bool recompute_if_invalid) const
{
  if (recompute_if_invalid)
    updateKDTree();

  return *kdtree;
}

Model::VertexKDTree const &
Model::getVertexKDTree(bool recompute_if_invalid) const
{
  if (recompute_if_invalid)
    updateVertexKDTree();

  return *vertex_kdtree;
}

bool
Model::rayIntersects(Ray3 const & ray, Real max_time) const
{
  return getKDTree().rayIntersects<Algorithms::RayIntersectionTester>(ray, max_time);
}

Real
Model::rayIntersectionTime(Ray3 const & ray, Real max_time) const
{
  return getKDTree().rayIntersectionTime<Algorithms::RayIntersectionTester>(ray, max_time);
}

Model::RayStructureIntersection3
Model::rayIntersection(Ray3 const & ray, Real max_time) const
{
  return getKDTree().rayStructureIntersection<Algorithms::RayIntersectionTester>(ray, max_time);
}

long
Model::closestPoint(Vector3 const & query, Real distance_bound, Real * min_dist, Vector3 * closest_pt,
                    Vector3 * closest_pt_normal, bool accelerate_with_vertices) const
{
  using namespace Algorithms;

  if (!isEmpty())
  {
    if (accelerate_with_vertices)
    {
      // Tighten the bound as much as we can with a fast initial query on the set of vertices
      updateVertexKDTree();
      double fast_distance_bound = 0;
      long vertex_index = vertex_kdtree->closestElement<MetricL2>(query, distance_bound, &fast_distance_bound);
      if (vertex_index >= 0)
        distance_bound = (Real)fast_distance_bound;
    }

    updateKDTree();
    double d = 0;
    long index = kdtree->closestElement<MetricL2>(query, distance_bound, &d, closest_pt);
    if (index >= 0)
    {
      if (min_dist) *min_dist = (Real)d;
      if (closest_pt_normal) *closest_pt_normal = kdtree->getElements()[index].getNormal();
      return index;
    }
  }

  return -1;
}

Real
Model::pick(Ray3 const & ray)
{
  RayStructureIntersection3 isec = rayIntersection(ray);
  if (isec.isValid())
  {
    array_size_t index = (array_size_t)isec.getElementIndex();
    KDTree::VertexTriple const & triple = kdtree->getElements()[index].getVertices();
    picked_sample.mesh = triple.getMesh();
    picked_sample.face_index = triple.meshFaceIsTriangle() ? picked_sample.mesh->getTriFaceIndex(triple.getMeshFaceIndex())
                                                           : picked_sample.mesh->getQuadFaceIndex(triple.getMeshFaceIndex());
    picked_sample.position = ray.getPoint(isec.getTime());

    valid_pick = true;

    emit needsRedraw(this);
  }

  return isec.getTime();
}

void
Model::invalidatePick()
{
  valid_pick = false;
  emit needsRedraw(this);
}

void
Model::mousePressEvent(QMouseEvent * event)
{
  event->accept();
}

void
Model::mouseMoveEvent(QMouseEvent * event)
{
  // Currently no-op
}

void
Model::mouseReleaseEvent(QMouseEvent * event)
{
  // Currently no-op
}

void
Model::addSample(Sample const & sample)
{
  samples.push_back(sample);
  emit needsRedraw(this);
}

bool
Model::addPickedSample(QString const & label, bool snap_to_vertex)
{
  if (valid_pick)
  {
    Sample sample = picked_sample;
    sample.label = label;

    if (snap_to_vertex)
    {
      Mesh::Face const & face = Mesh::mapIndexToFace(sample.face_index);
      if (!face)
      {
        THEA_ERROR << getName() << ": Mesh face with index " << sample.face_index << " not found";
        return false;
      }

      Mesh * mesh = (Mesh *)face.getMesh();
      Mesh::VertexArray const & verts = mesh->getVertices();
      Vector3 face_verts[3];
      if (face.hasTriangles())
      {
        Mesh::IndexTriple triple = mesh->getTriangle(face.getFirstTriangle());
        face_verts[0] = verts[(array_size_t)triple[0]];
        face_verts[1] = verts[(array_size_t)triple[1]];
        face_verts[2] = verts[(array_size_t)triple[2]];
      }
      else if (face.hasQuads())  // use only the first vertices, so all barycentric coords will be non-negative
      {
        Mesh::IndexQuad quad = mesh->getQuad(face.getFirstQuad());
        face_verts[0] = verts[(array_size_t)quad[0]];
        face_verts[1] = verts[(array_size_t)quad[1]];
        face_verts[2] = verts[(array_size_t)quad[2]];
      }
      else
      {
        THEA_ERROR << getName() << ": Face " << sample.face_index << " has neither triangles nor quads";
        return false;
      }

      int nn = 0;
      Real min_sqdist = (face_verts[0] - sample.position).squaredLength();
      Real sqdist = (face_verts[1] - sample.position).squaredLength(); if (sqdist < min_sqdist) { min_sqdist = sqdist; nn = 1; }
      sqdist      = (face_verts[2] - sample.position).squaredLength(); if (sqdist < min_sqdist) { min_sqdist = sqdist; nn = 2; }

      sample.position = face_verts[nn];
    }

    samples.push_back(sample);

    saveSamples(getSamplesFilename());
  }

  return valid_pick;
}

void
Model::removeSample(long index)
{
  if (index >= 0 && index < (long)samples.size())
  {
    samples.erase(samples.begin() + index);
    saveSamples(getSamplesFilename());
    emit needsRedraw(this);
  }
}

void
Model::selectSample(long index)
{
  selected_sample = index;
  emit needsRedraw(this);
}

bool
Model::loadSamples(QString const & filename_)
{
  samples.clear();
  bool status = true;
  try
  {
    std::ifstream in(toStdString(filename_).c_str());
    if (!in)
      throw Error("Could not open file");

    std::string line;
    if (!std::getline(in, line))
      throw Error("Could not read first line");

    std::istringstream header_iss(line);
    long n;
    header_iss >> n;
    if (n < 0 || !header_iss)
      throw Error("Could not read valid number of samples");

    std::string type, label;
    long face_index;
    double bary[3], coords[3];
    Vector3 face_verts[3];
    Vector3 pos;

    samples.resize((array_size_t)n);

    for (long i = 0; i < n; ++i)
    {
      if (!std::getline(in, line))
        throw Error("Could not read line");

      std::istringstream iss(line);

      iss >> type >> face_index >> bary[0] >> bary[1] >> bary[2];

      if (!(iss >> label))
        label = "";

      iss >> coords[0] >> coords[1] >> coords[2];
      bool read_pos = (bool)iss;

      Mesh::Face const & face = Mesh::mapIndexToFace(face_index);
      if (!face)
        throw Error(format("Mesh face with index %ld not found", face_index));

      Mesh * mesh = (Mesh *)face.getMesh();

      if (read_pos)
        pos = Vector3((Real)coords[0], (Real)coords[1], (Real)coords[2]);
      else
      {
        Mesh::VertexArray const & verts = mesh->getVertices();
        if (face.hasTriangles())
        {
          Mesh::IndexTriple triple = mesh->getTriangle(face.getFirstTriangle());
          face_verts[0] = verts[(array_size_t)triple[0]];
          face_verts[1] = verts[(array_size_t)triple[1]];
          face_verts[2] = verts[(array_size_t)triple[2]];
        }
        else if (face.hasQuads())
        {
          Mesh::IndexQuad quad = mesh->getQuad(face.getFirstQuad());
          face_verts[0] = verts[(array_size_t)quad[0]];
          face_verts[1] = verts[(array_size_t)quad[1]];
          face_verts[2] = verts[(array_size_t)quad[2]];
        }
        else
          throw Error(format("Face %ld has neither triangles nor quads", face_index));

        pos = bary[0] * face_verts[0] + bary[1] * face_verts[1] + bary[2] * face_verts[2];
      }

      samples[(array_size_t)i] = Sample(mesh, face_index, pos, toQString(type), toQString(label));
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(status = false;, WARNING, "Couldn't load model samples from '%s'",
                             toStdString(filename_).c_str())

  emit needsSyncSamples(this);
  return status;
}

QString
Model::getSamplesFilename() const
{
  return getFilename() + ".featpts";
}

bool
Model::saveSamples(QString const & filename_) const
{
  std::ofstream out(toStdString(filename_).c_str(), std::ios::binary);
  if (!out)
    return false;

  out << samples.size() << std::endl;

  Real bary[3];
  array_size_t vx_index, v[3];

  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    Sample const & sample = samples[i];
    Mesh::Face const & face = Mesh::mapIndexToFace(sample.face_index);
    if (!face)
    {
      THEA_ERROR << getName() << ": Mesh face with index " << sample.face_index << " not found, aborting saving samples";
      return false;
    }

    if (face.getMesh() != sample.mesh)
    {
      THEA_ERROR << getName() << ": Face " << sample.face_index << " belongs to wrong mesh";
      return false;
    }

    Mesh const * mesh = (Mesh const *)face.getMesh();
    if (face.hasTriangles())
    {
      vx_index = 3 * (array_size_t)face.getFirstTriangle();
      v[0] = (array_size_t)mesh->getTriangleIndices()[vx_index    ];
      v[1] = (array_size_t)mesh->getTriangleIndices()[vx_index + 1];
      v[2] = (array_size_t)mesh->getTriangleIndices()[vx_index + 2];
    }
    else if (face.hasQuads())
    {
      vx_index = 4 * (array_size_t)face.getFirstQuad();
      v[0] = (array_size_t)mesh->getTriangleIndices()[vx_index    ];
      v[1] = (array_size_t)mesh->getTriangleIndices()[vx_index + 1];
      v[2] = (array_size_t)mesh->getTriangleIndices()[vx_index + 2];
    }
    else
    {
      THEA_ERROR << getName() << ": Face " << sample.face_index << " has neither triangles nor quads";
      return false;
    }

    Mesh::VertexArray const & verts = mesh->getVertices();
    getBarycentricCoordinates3(sample.position, verts[v[0]], verts[v[1]], verts[v[2]], bary[0], bary[1], bary[2]);

    // Sometimes 0 gets written as -0
    if (bary[0] <= 0 && bary[0] >= -1.0e-06) bary[0] = 0;
    if (bary[1] <= 0 && bary[1] >= -1.0e-06) bary[1] = 0;
    if (bary[2] <= 0 && bary[2] >= -1.0e-06) bary[2] = 0;

    out << toStdString(sample.type) << ' ' << sample.face_index << ' '
        << (double)bary[0] << ' ' << (double)bary[1] << ' ' << (double)bary[2] << ' ' << toStdString(sample.label) << ' '
        << (double)sample.position[0] << ' ' << (double)sample.position[1] << ' ' << (double)sample.position[2] << std::endl;
  }

  return true;
}

AxisAlignedBox3 const &
Model::getBounds() const
{
  return bounds;
}

void
Model::updateBounds()
{
  bounds.setNull();

  if (mesh_group)
  {
    mesh_group->updateBounds();
    bounds.merge(mesh_group->getBounds());
  }

  if (!points.empty())
  {
    point_bounds.setNull();
    for (array_size_t i = 0; i < points.size(); ++i)
      point_bounds.merge(points[i]);

    bounds.merge(point_bounds);
  }
}

void
Model::uploadToGraphicsSystem(Graphics::RenderSystem & render_system)
{
  if (mesh_group)
    mesh_group->uploadToGraphicsSystem(render_system);
}

void
Model::draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
{
  if (isEmpty())
    return;

  const_cast<Model *>(this)->uploadToGraphicsSystem(render_system);

  static Color4 const DEFAULT_COLOR(1.0f, 0.9f, 0.8f, 1.0f);
  static Color4 const POINT_COLOR(1.0f, 1.0f, 0.5f, 1.0f);
  GraphicsWidget::setLight(Vector3(-1, -1, -2), Color3(1, 1, 1), Color3(1, 0.8f, 0.7f));

  render_system.pushShader();
    render_system.pushColorFlags();

      setPhongShader(render_system);

      if (app().getMainWindow()->pickPoints())
      {
        Real sample_radius = 0.005f * getBounds().getExtent().length();
        if (valid_pick)
        {
          render_system.setColor(Color3::red());
          drawSphere(render_system, picked_sample.position, sample_radius);
        }

        for (array_size_t i = 0; i < samples.size(); ++i)
        {
          render_system.setColor(getLabelColor(samples[i].label));

          if ((long)i == selected_sample)
            drawSphere(render_system, samples[i].position, 3 * sample_radius);
          else
            drawSphere(render_system, samples[i].position, sample_radius);
        }
      }

      render_system.setColor(DEFAULT_COLOR);
      if (mesh_group) mesh_group->draw(render_system, options);

      render_system.setColor(POINT_COLOR);
      Real point_radius = 0.002f * point_bounds.getExtent().length();
      for (array_size_t i = 0; i < points.size(); ++i)
        drawSphere(render_system, points[i], point_radius);

    render_system.popColorFlags();
  render_system.popShader();
}

} // namespace Browse3D
