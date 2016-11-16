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
#include "ModelDisplay.hpp"
#include "PointCloud.hpp"
#include "Util.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/RayIntersectionTester.hpp"
#include "../../Graphics/MeshCodec.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include <wx/filedlg.h>
#include <algorithm>
#include <fstream>

wxDEFINE_EVENT(EVT_MODEL_PATH_CHANGED,         wxCommandEvent);
wxDEFINE_EVENT(EVT_MODEL_GEOMETRY_CHANGED,     wxCommandEvent);
wxDEFINE_EVENT(EVT_MODEL_NEEDS_REDRAW,         wxCommandEvent);
wxDEFINE_EVENT(EVT_MODEL_NEEDS_SYNC_SAMPLES,   wxCommandEvent);
wxDEFINE_EVENT(EVT_MODEL_NEEDS_SYNC_SEGMENTS,  wxCommandEvent);

namespace Browse3D {

namespace ModelInternal {

bool first_file_dialog = true;

std::string
getWorkingDir()
{
  if (first_file_dialog)
    if (!app().options().working_dir.empty() && FileSystem::directoryExists(app().options().working_dir))
      return app().options().working_dir;

  return std::string();
}

bool
enableGPURendering(Mesh & mesh)
{
  mesh.setGPUBufferedRendering(true);
  mesh.setGPUBufferedWireframe(true);
  return false;
}

bool
disableGPURendering(Mesh & mesh)
{
  mesh.setGPUBufferedRendering(false);
  mesh.setGPUBufferedWireframe(false);
  return false;
}

void
linkMeshesToParent(MeshGroupPtr mesh_group)
{
  for (MeshGroup::MeshConstIterator mi = mesh_group->meshesBegin(); mi != mesh_group->meshesEnd(); ++mi)
    (*mi)->setParent(mesh_group.get());

  for (MeshGroup::GroupConstIterator ci = mesh_group->childrenBegin(); ci != mesh_group->childrenEnd(); ++ci)
    linkMeshesToParent(*ci);
}

static ColorRGBA const DEFAULT_COLOR(1.0f, 0.9f, 0.8f, 1.0f);
static ColorRGBA const PICKED_SEGMENT_COLOR(0.4f, 0.69f, 0.21f, 1.0f);

} // namespace ModelInternal

Model::Model(std::string const & initial_mesh)
: has_features(false),
  has_face_labels(false),
  color(ModelInternal::DEFAULT_COLOR),
  valid_pick(false),
  selected_sample(-1),
  segment_depth_promotion(0),
  selected_segment(-1),
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

std::string
Model::getName() const
{
  if (mesh_group)
    return mesh_group->getName();
  else if (point_cloud)
    return point_cloud->getName();
  else
    return "Untitled";
}

bool
Model::isEmpty() const
{
  return (!mesh_group || mesh_group->isEmpty()) && (!point_cloud || point_cloud->isEmpty());
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
  has_features = false;
  has_face_labels = false;
  samples.clear();
  segments.clear();
}

void
Model::clearPoints()
{
  if (point_cloud) point_cloud->clear();
}

bool
Model::load(std::string path_)
{
  if (path_.empty())
    return false;

  path_ = FileSystem::resolve(path_);
  if (!FileSystem::fileExists(path_) || path_ == FileSystem::resolve(path))
    return false;

  if (endsWith(toLower(path_), ".pts"))
  {
    clear();

    point_cloud = PointCloudPtr(new PointCloud);
    if (!point_cloud->load(path_))
      return false;

    bounds = point_cloud->getBounds();
    path = path_;
  }
  else
  {
    MeshGroupPtr new_mesh_group(new MeshGroup("Mesh Group"));

    Mesh::resetVertexIndices();  // reset counting
    Mesh::resetFaceIndices();

    static CodecOBJ<Mesh> const obj_codec(CodecOBJ<Mesh>::ReadOptions().setIgnoreTexCoords(true));
    try
    {
      if (endsWith(toLower(path_), ".obj"))
        new_mesh_group->load(path_, obj_codec);
      else
        new_mesh_group->load(path_);
    }
    THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "Couldn't load model '%s'", path_.c_str())

    invalidateAll();

    ModelInternal::linkMeshesToParent(new_mesh_group);

    mesh_group = new_mesh_group;
    clearPoints();

    bounds = new_mesh_group->getBounds();

    THEA_CONSOLE << "Loaded model '" << path_ << "' with bounding box " << mesh_group->getBounds().toString();

    path = path_;

    loadSamples(getSamplesPath());
    loadSegments(getSegmentsPath());
    loadFeatures(getDefaultFeaturesPath());
    loadFaceLabels(getDefaultFaceLabelsPath());
  }

  wxPostEvent(this, wxCommandEvent(EVT_MODEL_PATH_CHANGED));
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_GEOMETRY_CHANGED));

  return true;
}

bool
Model::selectAndLoad()
{
  wxFileDialog file_dialog(app().getMainWindow(), "Load model", "", "",
                           "Model files (*.3ds *.obj *.off *.off.bin *.ply *.pts)|*.3ds;*.obj;*.off;*.off.bin;*.ply;*.pts",
                           wxFD_OPEN | wxFD_FILE_MUST_EXIST);
  if (file_dialog.ShowModal() == wxID_CANCEL)
      return false;

  bool success = load(file_dialog.GetPath().ToStdString());
  if (success)
    ModelInternal::first_file_dialog = false;

  return success;
}

void
Model::setTransform(AffineTransform3 const & trans_)
{
  TransformableBaseT::setTransform(trans_);

  if (valid_kdtree)
    kdtree->setTransform(trans_);

  if (valid_vertex_kdtree)
    vertex_kdtree->setTransform(trans_);
}

void
Model::clearTransform()
{
  TransformableBaseT::clearTransform();

  if (valid_kdtree)
    kdtree->clearTransform();

  if (valid_vertex_kdtree)
    vertex_kdtree->clearTransform();
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

    if (hasTransform())
      kdtree->setTransform(getTransform());

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
  CollectVerticesFunctor(TheaArray<MeshVertex *> * verts_) : verts(verts_) {}

  bool operator()(Mesh & mesh)
  {
    for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      verts->push_back(&(*vi));

    return false;
  }

  TheaArray<MeshVertex *> * verts;

}; // struct CollectVerticesFunctor

} // namespace ModelInternal

void
Model::updateVertexKDTree() const
{
  if (valid_vertex_kdtree) return;

  vertex_kdtree->clear(false);

  TheaArray<MeshVertex *> verts;
  ModelInternal::CollectVerticesFunctor func(&verts);
  mesh_group->forEachMeshUntil(&func);
  vertex_kdtree->init(verts.begin(), verts.end());

  if (hasTransform())
    vertex_kdtree->setTransform(getTransform());

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
    picked_sample.mesh = const_cast<Mesh *>(triple.getMesh());
    picked_sample.face_index = triple.getMeshFace()->attr().getIndex();
    picked_sample.position = ray.getPoint(isec.getTime());

    valid_pick = true;

    wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
  }

  return isec.getTime();
}

void
Model::invalidatePick()
{
  valid_pick = false;
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

void
Model::mousePressEvent(wxMouseEvent & event)
{
  event.StopPropagation();
}

void
Model::mouseMoveEvent(wxMouseEvent & event)
{
  // Currently no-op
}

void
Model::mouseReleaseEvent(wxMouseEvent & event)
{
  // Currently no-op
}

namespace ModelInternal {

std::istream &
getNextNonBlankLine(std::istream & in, std::string & line)
{
  while (std::getline(in, line))
  {
    if (!trimWhitespace(line).empty())
      break;
  }

  return in;
}

} // namespace ModelInternal

void
Model::addSample(Sample const & sample)
{
  samples.push_back(sample);
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

bool
Model::addPickedSample(std::string const & label, bool snap_to_vertex)
{
  if (valid_pick)
  {
//     if (label.isEmpty())
//     {
//       THEA_WARNING << getName() << ": Empty label, cannot add sample";
//       return false;
//     }

    Sample sample = picked_sample;
    sample.label = label;

    if (snap_to_vertex)
    {
      Mesh::Face const * face = Mesh::mapIndexToFace(sample.face_index);
      if (!face)
      {
        THEA_ERROR << getName() << ": Mesh face with index " << sample.face_index << " not found";
        return false;
      }

      MeshVertex const * nnv = NULL;
      Real min_sqdist = -1;
      for (MeshFace::VertexConstIterator fvi = face->verticesBegin(); fvi != face->verticesEnd(); ++fvi)
      {
        Real sqdist = ((*fvi)->getPosition() - sample.position).squaredLength();
        if (!nnv || sqdist < min_sqdist)
        {
          min_sqdist = sqdist;
          nnv = *fvi;
        }
      }

      sample.position = nnv->getPosition();
    }

    samples.push_back(sample);

    saveSamples(getSamplesPath());
  }

  return valid_pick;
}

void
Model::removeSample(long index)
{
  if (index >= 0 && index < (long)samples.size())
  {
    samples.erase(samples.begin() + index);
    saveSamples(getSamplesPath());
    selected_sample = -1;  // since the removal in general messes up the indices
    wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
  }
}

void
Model::selectSample(long index)
{
  selected_sample = index;
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

bool
Model::loadSamples(std::string const & path_)
{
  using namespace ModelInternal;

  samples.clear();
  bool status = true;
  try
  {
    std::ifstream in(path_.c_str());
    if (!in)
      throw Error("Couldn't open samples file '" + path_ + "'");

    std::string line;
    if (!getNextNonBlankLine(in, line))
      throw Error("Couldn't read first line");

    std::istringstream header_iss(line);
    long n;
    header_iss >> n;
    if (n < 0 || !header_iss)
      throw Error("Couldn't read valid number of samples");

    std::string type, label;
    long face_index;
    double bary[3], coords[3];
    Vector3 pos;

    samples.resize((array_size_t)n);

    for (long i = 0; i < n; ++i)
    {
      if (!getNextNonBlankLine(in, line))
        throw Error("Couldn't read line");

      std::istringstream iss(line);

      iss >> type >> face_index >> bary[0] >> bary[1] >> bary[2];

      bool read_pos = false;
      if (!(iss >> label))
        label = "";
      else
        read_pos = (bool)(iss >> coords[0] >> coords[1] >> coords[2]);

      Mesh::Face const * face = Mesh::mapIndexToFace(face_index);
      if (!face)
        throw Error(format("Mesh face with index %ld not found", face_index));

      Mesh * mesh = face->attr().getParent();

      if (read_pos)
        pos = Vector3((Real)coords[0], (Real)coords[1], (Real)coords[2]);
      else
      {
        if (face->numVertices() < 3)
          throw Error(format("Face %ld has %d vertices", face->attr().getIndex(), face->numVertices()));

        MeshFace::VertexConstIterator v2 = face->verticesBegin();
        MeshFace::VertexConstIterator v0 = v2; ++v2;
        MeshFace::VertexConstIterator v1 = v2; ++v2;

        pos = bary[0] * ((*v0)->getPosition())
            + bary[1] * ((*v1)->getPosition())
            + bary[2] * ((*v2)->getPosition());
      }

      samples[(array_size_t)i] = Sample(mesh, face_index, pos, type, label);
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(status = false;, WARNING, "Couldn't load model samples from '%s'", path_.c_str())

  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_SYNC_SAMPLES));

  return status;
}

bool
Model::saveSamples(std::string const & path_) const
{
  std::ofstream out(path_.c_str(), std::ios::binary);
  if (!out)
    return false;

  out << samples.size() << std::endl;

  Real bary[3];

  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    Sample const & sample = samples[i];
    Mesh::Face const * face = Mesh::mapIndexToFace(sample.face_index);
    if (!face)
    {
      THEA_ERROR << getName() << ": Mesh face with index " << sample.face_index << " not found, aborting saving samples";
      return false;
    }

    if (face->attr().getParent() != sample.mesh)
    {
      THEA_ERROR << getName() << ": Face " << sample.face_index << " belongs to wrong mesh";
      return false;
    }

    if (face->numVertices() < 3)
      throw Error(format("Face %ld has %d vertices", face->attr().getIndex(), face->numVertices()));

    MeshFace::VertexConstIterator v2 = face->verticesBegin();
    MeshFace::VertexConstIterator v0 = v2; ++v2;
    MeshFace::VertexConstIterator v1 = v2; ++v2;

    getBarycentricCoordinates3(sample.position, (*v0)->getPosition(), (*v1)->getPosition(), (*v2)->getPosition(),
                               bary[0], bary[1], bary[2]);

    // Sometimes 0 gets written as -0
    if (bary[0] <= 0 && bary[0] >= -1.0e-06) bary[0] = 0;
    if (bary[1] <= 0 && bary[1] >= -1.0e-06) bary[1] = 0;
    if (bary[2] <= 0 && bary[2] >= -1.0e-06) bary[2] = 0;

    out << sample.type << ' ' << sample.face_index << ' '
        << (double)bary[0] << ' ' << (double)bary[1] << ' ' << (double)bary[2] << ' ' << sample.label << ' '
        << (double)sample.position[0] << ' ' << (double)sample.position[1] << ' ' << (double)sample.position[2] << std::endl;
  }

  return true;
}

std::string
Model::getSamplesPath() const
{
  std::string sfn = getPath() + ".picked";
  if (FileSystem::fileExists(sfn))
    return sfn;
  else
    return FilePath::concat(FilePath::parent(getPath()), FilePath::baseName(getPath()) + ".picked");
}

namespace ModelInternal {

struct SimilarComponentCollector
{
  SimilarComponentCollector() : query_mesh(NULL), query_group(NULL) {}

  void setQuery(Mesh const * mesh) { query_mesh = mesh; }
  void setQuery(MeshGroup const * group) { query_group = group; }

  void collectSimilarIn(MeshGroup const & root)
  {
    if (!query_mesh && !query_group)
      return;

    bool is_similar = (query_group ? isSimilarTo(root, *query_group) : isSimilarTo(root, *query_mesh));
    if (is_similar)
    {
      similar_groups.push_back(&root);
      return;
    }

    for (MeshGroup::GroupConstIterator ci = root.childrenBegin(); ci != root.childrenEnd(); ++ci)
    {
      is_similar = (query_group ? isSimilarTo(**ci, *query_group) : isSimilarTo(**ci, *query_mesh));
      if (is_similar)
        similar_groups.push_back(ci->get());
      else
        collectSimilarIn(**ci);
    }

    for (MeshGroup::MeshConstIterator mi = root.meshesBegin(); mi != root.meshesEnd(); ++mi)
    {
      is_similar = (query_group ? isSimilarTo(**mi, *query_group) : isSimilarTo(**mi, *query_mesh));
      if (is_similar)
        similar_meshes.push_back(mi->get());
    }
  }

  Mesh const * query_mesh;
  MeshGroup const * query_group;

  TheaArray<Mesh const *> similar_meshes;
  TheaArray<MeshGroup const *> similar_groups;
};

} // namespace ModelInternal

Real
Model::togglePickMesh(Ray3 const & ray, bool extend_to_similar)
{
  RayStructureIntersection3 isec = rayIntersection(ray);
  if (isec.isValid())
  {
    array_size_t index = (array_size_t)isec.getElementIndex();
    KDTree::VertexTriple const & triple = kdtree->getElements()[index].getVertices();
    Mesh * mesh = const_cast<Mesh *>(triple.getMesh());

    Segment const * existing = getSegment(mesh);
    if (existing)
    {
      THEA_WARNING << "Cannot pick mesh, it is already in another segment with label '" << existing->getLabel() << '\'';
      return -1;
    }

    bool add = true;
    if (picked_segment.hasMesh(mesh, segment_depth_promotion))
    {
      picked_segment.removeMesh(mesh, segment_depth_promotion);
      add = false;
      THEA_CONSOLE << "Removed mesh '" << mesh->getName() << "' from picked segment";
    }
    else
    {
      picked_segment.addMesh(mesh);
      THEA_CONSOLE << "Added mesh '" << mesh->getName() << "' to picked segment";
    }

    if (extend_to_similar)
    {
      ModelInternal::SimilarComponentCollector scc;
      if (segment_depth_promotion <= 0)
        scc.setQuery(mesh);
      else
      {
        MeshGroup const * ancestor = mesh->getAncestor(segment_depth_promotion);
        scc.setQuery(ancestor);
      }

      scc.collectSimilarIn(*mesh_group);

      for (array_size_t i = 0; i < scc.similar_meshes.size(); ++i)
      {
        if (add)
          picked_segment.addMesh(const_cast<Mesh *>(scc.similar_meshes[i]));
        else
          picked_segment.removeMesh(scc.similar_meshes[i]);
      }

      for (array_size_t i = 0; i < scc.similar_groups.size(); ++i)
      {
        if (add)
          picked_segment.addMeshGroup(const_cast<MeshGroup *>(scc.similar_groups[i]));
        else
          picked_segment.removeMeshGroup(scc.similar_groups[i]);
      }
    }

    wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
  }

  return isec.getTime();
}

void
Model::promotePickedSegment(long offset)
{
  segment_depth_promotion += offset;

  if (segment_depth_promotion < 0)
    segment_depth_promotion = 0;

  long min_depth = picked_segment.minDepth();
  if (min_depth >= 0 && segment_depth_promotion >= min_depth)
    segment_depth_promotion = std::max(min_depth - 1, 0L);

  THEA_CONSOLE << getName() << ": Segment depth promotion set to " << segment_depth_promotion;

  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

void
Model::addSegment(Segment const & segment)
{
  segments.push_back(segment);
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

bool
Model::addPickedSegment(std::string const & label)
{
  if (label.empty())
  {
    THEA_WARNING << getName() << ": Empty label, cannot add segment";
    return false;
  }

  if (picked_segment.numMeshes() <= 0)
  {
    THEA_WARNING << getName() << ": Empty selection, cannot add segment";
    return false;
  }

  picked_segment.setLabel(label);
  segments.push_back(picked_segment);
  saveSegments(getSegmentsPath());

  return true;
}

void
Model::removeSegment(long index)
{
  if (index >= 0 && index < (long)segments.size())
  {
    segments.erase(segments.begin() + index);
    saveSegments(getSegmentsPath());
    selected_segment = -1;  // since the removal in general messes up the indices
    wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
  }
}

Segment *
Model::getSegment(Mesh const * mesh, int * index)
{
  for (array_size_t i = 0; i < segments.size(); ++i)
    if (segments[i].hasMesh(mesh, segment_depth_promotion))
    {
      if (index) *index = (int)i;
      return &segments[i];
    }

  if (index) *index = -1;
  return NULL;
}

void
Model::selectSegment(long index)
{
  selected_segment = index;
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
}

bool
Model::loadSegments(std::string const & path_)
{
  using namespace ModelInternal;

  segments.clear();

  bool status = true;
  try
  {
    std::ifstream in(path_.c_str());
    if (!in)
      throw Error("Couldn't open file");

    segments.clear();

    std::string line;
    while (getNextNonBlankLine(in, line))
    {
      std::string label = trimWhitespace(line);
      Segment seg(label);

      if (!getNextNonBlankLine(in, line))
        throw Error("Couldn't read list of representative faces");

      std::istringstream iss(line);
      long face_index = -1;
      while (iss >> face_index)
      {
        Mesh::Face const * face = Mesh::mapIndexToFace(face_index);
        if (!face)
          throw Error(format("Mesh face with index %ld not found", face_index));

        Mesh * mesh = face->attr().getParent();
        seg.addMesh(mesh);
      }

      if (seg.numMeshes() > 0)
        segments.push_back(seg);
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(status = false;, WARNING, "Couldn't load model segments from '%s'", path_.c_str())

  if (!status)
    segments.clear();

  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_SYNC_SEGMENTS));

  return status;
}

bool
Model::saveSegments(std::string const & path_) const
{
  std::ofstream out(path_.c_str(), std::ios::binary);
  if (!out)
    return false;

  for (array_size_t i = 0; i < segments.size(); ++i)
  {
    Segment const & seg = segments[i];
    out << seg.getLabel() << '\n';

    Segment::MeshSet const & meshes = seg.getMeshes();
    for (Segment::MeshSet::const_iterator mj = meshes.begin(); mj != meshes.end(); ++mj)
    {
      Mesh const * mesh = *mj;
      if (!mesh || mesh->numFaces() <= 0)
        continue;

      long face_index = mesh->facesBegin()->attr().getIndex();  // the first face

      if (mj != meshes.begin())
        out << ' ';

      out << face_index;
    }

    if (!meshes.empty())
      out << '\n';

    if (i + 1 < segments.size())
      out << '\n';
  }

  return true;
}

std::string
Model::getSegmentsPath() const
{
  std::string sfn = getPath() + ".labels";
  if (FileSystem::fileExists(sfn))
    return sfn;
  else
    return FilePath::concat(FilePath::parent(getPath()), FilePath::baseName(getPath()) + ".labels");
}

namespace ModelInternal {

typedef Algorithms::KDTreeN<Vector3, 3> PointKDTree;

struct VertexFeatureVisitor
{
  VertexFeatureVisitor(PointKDTree * fkdtree_, Real const * feat_vals0_, Real const * feat_vals1_, Real const * feat_vals2_)
  : fkdtree(fkdtree_), feat_vals0(feat_vals0_), feat_vals1(feat_vals1_), feat_vals2(feat_vals2_) {}

  bool operator()(Mesh & mesh)
  {
    for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    {
      long nn_index = fkdtree->closestElement<Algorithms::MetricL2>(vi->getPosition());
      if (nn_index >= 0)
      {
        if (!feat_vals2)
        {
          if (!feat_vals1)
            vi->attr().setColor(ColorRGB::jetColorMap(0.2 + 0.6 * feat_vals0[nn_index]));
          else
            vi->attr().setColor(ColorRGB(feat_vals0[nn_index], feat_vals1[nn_index], 1.0f));
        }
        else
          vi->attr().setColor(ColorRGB(feat_vals0[nn_index], feat_vals1[nn_index], feat_vals2[nn_index]));
      }
      else
        THEA_CONSOLE << "No nearest neighbor found!";
    }

    mesh.invalidateGPUBuffers(Mesh::BufferID::VERTEX_COLOR);

    return false;
  }

  PointKDTree * fkdtree;
  Real const * feat_vals0;
  Real const * feat_vals1;
  Real const * feat_vals2;
};

} // namespace ModelInternal

bool
Model::loadFeatures(std::string const & path_)
{
  using namespace ModelInternal;

  features_path = path_;

  if (point_cloud)
  {
    has_features = point_cloud->loadFeatures(path_);
    return has_features;
  }

  if (!mesh_group)
  {
    has_features = false;
    return has_features;
  }

  TheaArray<Vector3> feat_pts;
  TheaArray< TheaArray<Real> > feat_vals(1);
  has_features = true;
  try
  {
    std::ifstream in(path_.c_str());
    if (!in)
      throw Error("Couldn't open file");

    std::string line;
    Vector3 p;
    Real f;
    while (getNextNonBlankLine(in, line))
    {
      std::istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2] >> f))
        throw Error("Couldn't read feature");

      feat_pts.push_back(p);
      feat_vals[0].push_back(f);

      if (feat_pts.size() == 1)
      {
        while (line_in >> f)
        {
          feat_vals.push_back(TheaArray<Real>());
          feat_vals.back().push_back(f);
        }
      }
      else
      {
        for (array_size_t i = 1; i < feat_vals.size(); ++i)
        {
          if (!(line_in >> f))
            throw Error("Couldn't read feature");

          feat_vals[i].push_back(f);
        }
      }
    }

    if (feat_pts.empty())
    {
      wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));
      return true;
    }

    if (app().options().accentuate_features)
    {
      if (app().options().color_cube_features && feat_vals.size() == 3)
      {
        Real abs_max = -1;
        for (array_size_t i = 0; i < feat_vals.size(); ++i)
        {
          for (array_size_t j = 0; j < feat_vals[i].size(); ++j)
          {
            Real abs_feat_val = std::fabs(feat_vals[i][j]);
            if (abs_feat_val > abs_max)
              abs_max = abs_feat_val;
          }
        }

        if (abs_max > 0)
        {
          for (array_size_t i = 0; i < feat_vals.size(); ++i)
            for (array_size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp(0.5 * (feat_vals[i][j]  / abs_max + 1), (Real)0, (Real)1);
        }
      }
      else
      {
        for (array_size_t i = 0; i < feat_vals.size(); ++i)
        {
          TheaArray<Real> sorted = feat_vals[i];
          std::sort(sorted.begin(), sorted.end());

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
            for (array_size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp((feat_vals[i][j] - lo) / range, (Real)0, (Real)1);
          }
          else
          {
            Real abs_max = std::max(std::fabs(lo), std::fabs(hi));
            for (array_size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp((feat_vals[i][j] + abs_max) / (2 * abs_max), (Real)0, (Real)1);
          }
        }
      }
    }

    PointKDTree fkdtree(feat_pts.begin(), feat_pts.end());
    VertexFeatureVisitor visitor(&fkdtree,
                                 &feat_vals[0][0],
                                 feat_vals.size() > 1 ? &feat_vals[1][0] : NULL,
                                 feat_vals.size() > 2 ? &feat_vals[2][0] : NULL);
    mesh_group->forEachMeshUntil(&visitor);
  }
  THEA_STANDARD_CATCH_BLOCKS(has_features = false;, WARNING, "Couldn't load model features from '%s'", path_.c_str())

  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));

  return has_features;
}

namespace ModelInternal {

class FaceLabeler
{
  public:
    FaceLabeler(TheaArray<ColorRGBA> const & face_colors_) : face_colors(face_colors_) {}

    bool operator()(Mesh & mesh) const
    {
      for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
      {
        long index = fi->attr().getIndex();
        if (index < 0 || index >= (long)face_colors.size())
          throw Error("Face index out of range of face labels array");

        fi->attr().setColor(face_colors[(array_size_t)index]);
      }

      return false;
    }

  private:
    TheaArray<ColorRGBA> const & face_colors;
};

} // namespace ModelInternal

bool
Model::loadFaceLabels(std::string const & path_)
{
  has_face_labels = false;

  if (!mesh_group)
    return has_face_labels;

  std::ifstream in(path_.c_str());
  if (!in)
  {
    THEA_WARNING << "Couldn't open face labels file '" << path_ << '\'';
    return has_face_labels;
  }

  TheaArray<ColorRGBA> face_colors;
  std::string line;
  while (in)
  {
    std::getline(in, line);
    line = trimWhitespace(line);

    // Check if this is an integer, if so, use as is.
    char * next;
    long n = strtol(line.c_str(), &next, 10);
    if (*next != 0)  // could not parse to end of string, not an integer
      face_colors.push_back(getLabelColor(line));
    else
      face_colors.push_back(getPaletteColor(n));
  }

  try
  {
    ModelInternal::FaceLabeler flab(face_colors);
    mesh_group->forEachMeshUntil(&flab);
  }
  THEA_STANDARD_CATCH_BLOCKS(return has_face_labels;, WARNING, "Couldn't load model face labels from '%s'", path_.c_str())

  has_face_labels = true;
  wxPostEvent(this, wxCommandEvent(EVT_MODEL_NEEDS_REDRAW));

  return has_face_labels;
}

namespace ModelInternal {

std::string
getDefaultPath(std::string model_path, std::string const & query_path, TheaArray<std::string> const & query_exts)
{
  if (FileSystem::fileExists(query_path))
    return query_path;

  int iter_begin = FileSystem::directoryExists(query_path) ? 0 : 1;

  for (int i = iter_begin; i < 2; ++i)
  {
    std::string dir = (i == 0 ? query_path : FilePath::parent(model_path));

    for (array_size_t j = 0; j < query_exts.size(); ++j)
    {
      std::string ffn = FilePath::concat(dir, model_path + query_exts[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }

    for (array_size_t j = 0; j < query_exts.size(); ++j)
    {
      std::string ffn = FilePath::concat(dir, FilePath::completeBaseName(model_path) + query_exts[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }

    for (array_size_t j = 0; j < query_exts.size(); ++j)
    {
      std::string ffn = FilePath::concat(dir, FilePath::baseName(model_path) + query_exts[j]);
      if (FileSystem::exists(ffn))
        return ffn;
    }
  }

  return "";
}

} // namespace ModelInternal

std::string
Model::getDefaultFeaturesPath() const
{
  TheaArray<std::string> exts;
  exts.push_back(".arff");
  exts.push_back(".features");

  return ModelInternal::getDefaultPath(path, app().options().features, exts);
}

std::string
Model::getDefaultFaceLabelsPath() const
{
  TheaArray<std::string> exts;
  exts.push_back(".seg");

  return ModelInternal::getDefaultPath(path, app().options().face_labels, exts);
}

void
Model::registerDisplay(ModelDisplay * display)
{
  if (!display) return;

  Bind(EVT_MODEL_GEOMETRY_CHANGED, &ModelDisplay::modelGeometryChanged, display);
  Bind(EVT_MODEL_NEEDS_REDRAW, &ModelDisplay::modelNeedsRedraw, display);
}

void
Model::deregisterDisplay(ModelDisplay * display)
{
  if (!display) return;

  Unbind(EVT_MODEL_GEOMETRY_CHANGED, &ModelDisplay::modelGeometryChanged, display);
  Unbind(EVT_MODEL_NEEDS_REDRAW, &ModelDisplay::modelNeedsRedraw, display);
}

AxisAlignedBox3 const &
Model::getBounds() const
{
  return bounds;
}

AxisAlignedBox3
Model::getTransformedBounds() const
{
  return hasTransform() ? bounds.transformAndBound(getTransform()) : bounds;
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

  if (point_cloud)
  {
    point_cloud->updateBounds();
    bounds.merge(point_cloud->getBounds());
  }
}

void
Model::uploadToGraphicsSystem(Graphics::RenderSystem & render_system)
{
  if (mesh_group)
    mesh_group->uploadToGraphicsSystem(render_system);

  if (point_cloud)
    point_cloud->uploadToGraphicsSystem(render_system);
}

void
Model::drawSegmentedMeshGroup(MeshGroupPtr mesh_group, int depth, int & node_index, Graphics::RenderSystem & render_system,
                              Graphics::RenderOptions const & options) const
{
  Graphics::RenderOptions ro = options;  // make a copy and tweak it

  ro.overrideEdgeColor() = true;

  for (MeshGroup::MeshConstIterator mi = mesh_group->meshesBegin(); mi != mesh_group->meshesEnd(); ++mi, ++node_index)
  {
    Mesh const * mesh = mi->get();
    if (!mesh) continue;

    int seg_index = -1;
    Segment const * seg = getSegment(mesh, &seg_index);
    if (seg)
    {
      if (seg_index >= 0 && seg_index == selected_segment)
      {
        ro.edgeColor() = ColorRGBA(1, 0, 0, 1);
        ro.drawEdges() = true;
      }
      else
        ro.drawEdges() = false;

      render_system.setColor(getLabelColor(seg->getLabel()));
    }
    else if (picked_segment.hasMesh(mesh, segment_depth_promotion))
    {
      ro.drawEdges() = false;
      render_system.setColor(ModelInternal::PICKED_SEGMENT_COLOR);
    }
    else
    {
      ro.drawEdges() = true;
      ro.edgeColor() = getPaletteColor(node_index);
      render_system.setColor(color);
    }

    mesh->draw(render_system, ro);
  }

  for (MeshGroup::GroupConstIterator ci = mesh_group->childrenBegin(); ci != mesh_group->childrenEnd(); ++ci)
    drawSegmentedMeshGroup(*ci, depth + 1, node_index, render_system, ro);
}

namespace ModelInternal {

struct DrawFaceNormals
{
  DrawFaceNormals(Graphics::RenderSystem * rs, Real normal_scale_) : render_system(rs), normal_scale(normal_scale_) {}

  bool operator()(Mesh const & mesh)
  {
    render_system->beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

      for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
      {
        if (fi->numVertices() <= 0)
          continue;

        Vector3 c(0, 0, 0);
        for (Mesh::Face::VertexConstIterator vi = fi->verticesBegin(); vi != fi->verticesEnd(); ++vi)
          c += (*vi)->getPosition();

        c /= fi->numVertices();

        render_system->sendVertex(c);
        render_system->sendVertex(c + normal_scale * fi->getNormal());
      }

    render_system->endPrimitive();

    return false;
  }

  Graphics::RenderSystem * render_system;
  Real normal_scale;
};

struct DrawVertexNormals
{
  DrawVertexNormals(Graphics::RenderSystem * rs, Real normal_scale_) : render_system(rs), normal_scale(normal_scale_) {}

  bool operator()(Mesh const & mesh)
  {
    render_system->beginPrimitive(Graphics::RenderSystem::Primitive::LINES);

      for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      {
        render_system->sendVertex(vi->getPosition());
        render_system->sendVertex(vi->getPosition() + normal_scale * vi->getNormal());
      }

    render_system->endPrimitive();

    return false;
  }

  Graphics::RenderSystem * render_system;
  Real normal_scale;
};

} // namespace ModelInternal

void
Model::draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
{
  if (isEmpty())
    return;

  const_cast<Model *>(this)->uploadToGraphicsSystem(render_system);

  GraphicsWidget::setLight(Vector3(-1, -1, -2), ColorRGB(1, 1, 1), ColorRGB(1, 0.8f, 0.7f));

  if (hasTransform())
  {
    render_system.setMatrixMode(Graphics::RenderSystem::MatrixMode::MODELVIEW); render_system.pushMatrix();
    render_system.multMatrix(getTransform().toHomMatrix());
  }

  render_system.pushShader();
  render_system.pushTextures();
  render_system.pushColorFlags();

    setPhongShader(render_system);
    render_system.setTexture(0, NULL);

    if (app().getMainWindow()->pickPoints())
    {
      Real sample_radius = 0.005f * getBounds().getExtent().length();
      if (valid_pick)
      {
        render_system.setColor(ColorRGB::red());
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

    render_system.setColor(color);

    if (mesh_group)
    {
      if (app().getMainWindow()->pickSegments())
      {
        int node_index = 0;
        drawSegmentedMeshGroup(mesh_group, 0, node_index, render_system, options);
      }
      else
      {
        Graphics::RenderOptions ro = options;  // make a copy

        if (has_features)
        {
          ro.sendColors() = true;
          ro.useVertexData() = true;
        }
        else if (has_face_labels)
        {
          ro.sendColors() = true;
          ro.useVertexData() = false;
        }

        bool smooth_shading = (ro.useVertexNormals() && ro.useVertexData());

        if (smooth_shading)
          mesh_group->forEachMeshUntil(&ModelInternal::enableGPURendering);
        else
          mesh_group->forEachMeshUntil(&ModelInternal::disableGPURendering);

        mesh_group->draw(render_system, ro);

        if (app().options().show_normals)
        {
          render_system.setShader(NULL);
          render_system.setColor(ColorRGB(0, 1, 0));

          Real normal_scale = 0.025f * getBounds().getExtent().length();

          if (smooth_shading)
          {
            ModelInternal::DrawVertexNormals drawer(&render_system, normal_scale);
            mesh_group->forEachMeshUntil(&drawer);
          }
          else
          {
            ModelInternal::DrawFaceNormals drawer(&render_system, normal_scale);
            mesh_group->forEachMeshUntil(&drawer);
          }
        }
      }
    }

  render_system.popColorFlags();
  render_system.popTextures();
  render_system.popShader();

  if (point_cloud)
    point_cloud->draw(render_system, options);

  if (hasTransform())
  {
    render_system.setMatrixMode(Graphics::RenderSystem::MatrixMode::MODELVIEW); render_system.popMatrix();
  }
}

} // namespace Browse3D
