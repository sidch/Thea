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

#include "DisplayMesh.hpp"
#include "../Polygon3.hpp"
#include "../UnorderedSet.hpp"

namespace Thea {
namespace Graphics {

DisplayMesh::DisplayMesh(std::string const & name)
: NamedObject(name),
  valid_bounds(true),
  wireframe_enabled(false),
  changed_buffers(BufferID::ALL),
  var_area(NULL),
  vertices_var(NULL),
  tris_var(NULL),
  quads_var(NULL),
  normals_var(NULL),
  colors_var(NULL),
  texcoords_var(NULL)
{}

DisplayMesh::DisplayMesh(DisplayMesh const & src)
: NamedObject(src),
  vertices(src.vertices),
  normals(src.normals),
  colors(src.colors),
  texcoords(src.texcoords),
  tris(src.tris),
  quads(src.quads),
  valid_bounds(src.valid_bounds),
  bounds(src.bounds),
  wireframe_enabled(src.wireframe_enabled),
  changed_buffers(BufferID::ALL),
  var_area(NULL),
  vertices_var(NULL),
  tris_var(NULL),
  quads_var(NULL),
  normals_var(NULL),
  colors_var(NULL),
  texcoords_var(NULL)
{}

void
DisplayMesh::clear()
{
  vertices.clear();
  normals.clear();
  colors.clear();
  texcoords.clear();
  tris.clear();
  quads.clear();
  edges.clear();

  vertex_source_indices.clear();
  tri_source_face_indices.clear();
  quad_source_face_indices.clear();

  face_vertex_indices.clear();
  triangulated_indices.clear();

  valid_bounds = true;
  bounds = AxisAlignedBox3();

  invalidateGPUBuffers();
}

DisplayMesh::Vertex
DisplayMesh::getVertex(long i)
{
  debugAssertM(i >= 0 && i < (long)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  size_t si = (size_t)i;
  return Vertex(this, vertices[si],
                hasNormals()    ?  &normals[si]    :  NULL,
                hasColors()     ?  &colors[si]     :  NULL,
                hasTexCoords()  ?  &texcoords[si]  :  NULL);
}

DisplayMesh::IndexTriple
DisplayMesh::getTriangle(long tri_index) const
{
  debugAssertM(tri_index >= 0 && 3 * tri_index < (long)tris.size(), getNameStr() + ": Triangle index out of bounds");

  size_t base_index = (size_t)(3 * tri_index);
  IndexTriple tri;
  tri[0] = (long)tris[base_index];
  tri[1] = (long)tris[base_index + 1];
  tri[2] = (long)tris[base_index + 2];

  return tri;
}

DisplayMesh::IndexQuad
DisplayMesh::getQuad(long quad_index) const
{
  debugAssertM(quad_index >= 0 && 4 * quad_index < (long)quads.size(), getNameStr() + ": Quad index out of bounds");

  size_t base_index = (size_t)(4 * quad_index);
  IndexQuad quad;
  quad[0] = (long)quads[base_index];
  quad[1] = (long)quads[base_index + 1];
  quad[2] = (long)quads[base_index + 2];
  quad[3] = (long)quads[base_index + 3];

  return quad;
}

void
DisplayMesh::addColors()
{
  if (colors.size() < vertices.size())
  {
    colors.resize(vertices.size(), ColorRGBA(0, 0, 0, 0));
    invalidateGPUBuffers();
  }
}

void
DisplayMesh::addNormals()
{
  if (normals.size() < vertices.size())
  {
    normals.resize(vertices.size(), Vector3::zero());
    invalidateGPUBuffers();
  }
}

void
DisplayMesh::addTexCoords()
{
  if (texcoords.size() < vertices.size())
  {
    texcoords.resize(vertices.size(), Vector2::zero());
    invalidateGPUBuffers();
  }
}

long
DisplayMesh::addVertex(Vector3 const & point, long source_index, Vector3 const * normal, ColorRGBA const * color,
                       Vector2 const * texcoord)
{
  alwaysAssertM((source_index >= 0 && vertex_source_indices.size() == vertices.size())
             || (source_index < 0 && vertex_source_indices.empty()),
                getNameStr() + ": Mesh must have all or no vertex source indices");
  alwaysAssertM((normal && normals.size() == vertices.size()) || (!normal && normals.empty()),
                getNameStr() + ": Mesh must have all or no normals");
  alwaysAssertM((color && colors.size() == vertices.size()) || (!color && colors.empty()),
                getNameStr() + ": Mesh must have all or no vertex colors");
  alwaysAssertM((texcoord && texcoords.size() == vertices.size()) || (!texcoord && texcoords.empty()),
                getNameStr() + ": Mesh must have all or no texture coordinates");

  long index = (long)vertices.size();

  if (valid_bounds)
    bounds.merge(point);

  vertices.push_back(point);
  if (source_index >= 0)  vertex_source_indices.push_back(source_index);
  if (normal)             normals.push_back(*normal);
  if (color)              colors.push_back(*color);
  if (texcoord)           texcoords.push_back(*texcoord);

  invalidateGPUBuffers();

  return index;
}

long
DisplayMesh::addTriangle(long vi0, long vi1, long vi2, long source_face_index)
{
  debugAssertM(vi0 >= 0 && vi1 >= 0 && vi2 >= 0
            && vi0 < (long)vertices.size()
            && vi1 < (long)vertices.size()
            && vi2 < (long)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  alwaysAssertM((source_face_index >= 0 && tri_source_face_indices.size() == tris.size())
             || (source_face_index < 0 && tri_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

  long index = (long)(tris.size() / 3);

  tris.push_back((uint32)vi0);
  tris.push_back((uint32)vi1);
  tris.push_back((uint32)vi2);

  if (source_face_index >= 0)
    tri_source_face_indices.push_back(source_face_index);

  invalidateGPUBuffers();

  return index;
}

long
DisplayMesh::addQuad(long vi0, long vi1, long vi2, long vi3, long source_face_index)
{
  debugAssertM(vi0 >= 0 && vi1 >= 0 && vi2 >= 0 && vi3 >= 0
            && vi0 < (long)vertices.size()
            && vi1 < (long)vertices.size()
            && vi2 < (long)vertices.size()
            && vi3 < (long)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  alwaysAssertM((source_face_index >= 0 && quad_source_face_indices.size() == quads.size())
             || (source_face_index < 0 && quad_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no quad face source indices");

  long index = (long)(quads.size() / 4);

  quads.push_back((uint32)vi0);
  quads.push_back((uint32)vi1);
  quads.push_back((uint32)vi2);
  quads.push_back((uint32)vi3);

  if (source_face_index >= 0)
    quad_source_face_indices.push_back(source_face_index);

  invalidateGPUBuffers();

  return index;
}

DisplayMesh::Face
DisplayMesh::addFace(int num_vertices, long const * face_vertex_indices_, long source_face_index)
{
  if (num_vertices < 3)
  {
    THEA_DEBUG << getName() << ": Skipping face -- too few vertices (" << num_vertices << ')';
    return Face();
  }

  if (num_vertices == 3)
    return Face(this, 3, true,
                addTriangle(face_vertex_indices_[0], face_vertex_indices_[1], face_vertex_indices_[2], source_face_index), 1);

  if (num_vertices == 4)
    return Face(this, 4, false,
                addQuad(face_vertex_indices_[0], face_vertex_indices_[1], face_vertex_indices_[2], face_vertex_indices_[3],
                        source_face_index),
                1);

  Polygon3 poly;
  for (int i = 0; i < num_vertices; ++i)
  {
    long vi = face_vertex_indices_[i];
    debugAssertM(vi >= 0 && vi < (long)vertices.size(), getName() + format(": Vertex index %ld out of bounds", vi));

    poly.addVertex(vertices[vi], vi);
  }

  long num_tris = poly.triangulate(triangulated_indices);

//   if ((long)triangulated_indices.size() != 3 * (num_vertices - 2))
//   {
//     THEA_ERROR << getName() << ": Triangulation of polygonal face with " << num_vertices << " vertices yielded "
//                << triangulated_indices.size() / 3.0f << " triangles, whereas " << num_vertices - 2 << " were expected";
//
//     for (int i = 0; i < num_vertices; ++i)
//       THEA_CONSOLE << "v[" << i << "] = " << vertices[face_vertex_indices_[i]];
//
//     throw FatalError("Triangulation error");
//   }

  if (num_tris <= 0)
    return Face();

  alwaysAssertM((source_face_index >= 0 && tri_source_face_indices.size() == tris.size())
             || (source_face_index < 0 && tri_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

  long starting_index = numTriangles();

  for (size_t i = 0; i < triangulated_indices.size(); ++i)
  {
    tris.push_back((uint32)triangulated_indices[i]);
    if (source_face_index >= 0) tri_source_face_indices.push_back(source_face_index);
  }

  invalidateGPUBuffers();

  return Face(this, num_vertices, true, starting_index, num_tris);
}

void
DisplayMesh::removeTriangle(long tri_index)
{
  debugAssertM(tri_index >= 0 && 3 * tri_index < (long)tris.size(), getNameStr() + ": Triangle index out of bounds");

  IndexArray::iterator ii = tris.begin() + 3 * tri_index;
  tris.erase(ii, ii + 3);

  invalidateGPUBuffers();
}

void
DisplayMesh::removeTriangles(long begin, long num_triangles)
{
  debugAssertM(begin >= 0 && 3 * begin < (long)tris.size(), getNameStr() + ": Triangle index out of bounds");

  IndexArray::iterator ii = tris.begin() + 3 * begin;
  tris.erase(ii, ii + 3 * num_triangles);

  invalidateGPUBuffers();
}

void
DisplayMesh::removeQuad(long quad_index)
{
  debugAssertM(quad_index >= 0 && 4 * quad_index < (long)quads.size(), getNameStr() + ": Quad index out of bounds");

  IndexArray::iterator ii = quads.begin() + 4 * quad_index;
  quads.erase(ii, ii + 4);

  invalidateGPUBuffers();
}

void
DisplayMesh::removeQuads(long begin, long num_quads)
{
  debugAssertM(begin >= 0 && 4 * begin < (long)quads.size(), getNameStr() + ": Quad index out of bounds");

  IndexArray::iterator ii = quads.begin() + 4 * begin;
  quads.erase(ii, ii + 4 * num_quads);

  invalidateGPUBuffers();
}

void
DisplayMesh::removeFace(Face const & face)
{
  if (!face)
    return;

  alwaysAssertM(face.getMesh() == this, getNameStr() + ": Face belongs to a different mesh");

  if (face.hasTriangles())
    removeTriangles(face.getFirstTriangle(), face.numTriangles());

  if (face.hasQuads())
    removeQuads(face.getFirstQuad(), face.numQuads());
}

void
DisplayMesh::computeAveragedVertexNormals()
{
  bool topo_change = (normals.size() != vertices.size());

  normals.resize(vertices.size());
  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = Vector3::zero();

  // TODO: weight normals by face area?
  Vector3 n;
  uint32 i0, i1, i2, i3;
  for (size_t i = 0; i < tris.size(); i += 3)
  {
    i0 = tris[i], i1 = tris[i + 1], i2 = tris[i + 2];
    Vector3 const & v0 = vertices[i0];
    Vector3 const & v1 = vertices[i1];
    Vector3 const & v2 = vertices[i2];

    n = (v2 - v1).cross(v0 - v1).unit();
    normals[i0] += n;
    normals[i1] += n;
    normals[i2] += n;
  }

  for (size_t i = 0; i < quads.size(); i += 4)
  {
    i0 = quads[i], i1 = quads[i + 1], i2 = quads[i + 2], i3 = quads[i + 3];
    Vector3 const & v0 = vertices[i0];
    Vector3 const & v1 = vertices[i1];
    Vector3 const & v2 = vertices[i2];

    n = (v2 - v1).cross(v0 - v1).unit();
    normals[i0] += n;
    normals[i1] += n;
    normals[i2] += n;
    normals[i3] += n;
  }

  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = normals[i].unit();

  invalidateGPUBuffers(topo_change ? BufferID::ALL : BufferID::NORMAL);
}

void
DisplayMesh::flipNormals()
{
  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = -normals[i];

  invalidateGPUBuffers(BufferID::NORMAL);
}

void
DisplayMesh::updateEdges()
{
  edges.clear();

  if (wireframe_enabled)
  {
    typedef std::pair<uint32, uint32> Edge;
    typedef TheaUnorderedSet<Edge> EdgeSet;

    EdgeSet added_edges;
    Edge edge;

    for (size_t i = 0; i < tris.size(); i += 3)
      for (size_t j = 0; j < 3; ++j)
      {
        uint32 ei0 = tris[i + j];
        uint32 ei1 = tris[j == 2 ? i : (i + j + 1)];

        // Order so lower index is first, since we're considering edges to be undirected
        if (ei0 < ei1)
          edge = Edge(ei0, ei1);
        else
          edge = Edge(ei1, ei0);

        EdgeSet::iterator existing = added_edges.find(edge);
        if (existing == added_edges.end())
        {
          edges.push_back(ei0);
          edges.push_back(ei1);
          added_edges.insert(edge);
        }
      }

    for (size_t i = 0; i < quads.size(); i += 4)
      for (size_t j = 0; j < 4; ++j)
      {
        uint32 ei0 = quads[i + j];
        uint32 ei1 = quads[j == 3 ? i : (i + j + 1)];

        // Order so lower index is first, since we're considering edges to be undirected
        if (ei0 < ei1)
          edge = Edge(ei0, ei1);
        else
          edge = Edge(ei1, ei0);

        EdgeSet::iterator existing = added_edges.find(edge);
        if (existing == added_edges.end())
        {
          edges.push_back(ei0);
          edges.push_back(ei1);
          added_edges.insert(edge);
        }
      }
  }

  THEA_DEBUG << getName() << ": Mesh has " << edges.size() / 2 << " edges";
}

void
DisplayMesh::isolateFaces()
{
  VertexArray    new_vertices;
  NormalArray    new_normals;
  ColorArray     new_colors;
  TexCoordArray  new_texcoords;

  for (size_t i = 0; i < tris.size(); i += 3)
  {
    size_t i0 = (size_t)tris[i    ];
    size_t i1 = (size_t)tris[i + 1];
    size_t i2 = (size_t)tris[i + 2];

    uint32 new_vindex = (uint32)new_vertices.size();

    new_vertices.push_back(vertices[i0]);
    new_vertices.push_back(vertices[i1]);
    new_vertices.push_back(vertices[i2]);

    if (!normals.empty())
    {
      new_normals.push_back(normals[i0]);
      new_normals.push_back(normals[i1]);
      new_normals.push_back(normals[i2]);
    }

    if (!colors.empty())
    {
      new_colors.push_back(colors[i0]);
      new_colors.push_back(colors[i1]);
      new_colors.push_back(colors[i2]);
    }

    if (!texcoords.empty())
    {
      new_texcoords.push_back(texcoords[i0]);
      new_texcoords.push_back(texcoords[i1]);
      new_texcoords.push_back(texcoords[i2]);
    }

    tris[i    ]  =  new_vindex;
    tris[i + 1]  =  new_vindex + 1;
    tris[i + 2]  =  new_vindex + 2;
  }

  for (size_t i = 0; i < quads.size(); i += 4)
  {
    size_t i0 = (size_t)quads[i    ];
    size_t i1 = (size_t)quads[i + 1];
    size_t i2 = (size_t)quads[i + 2];
    size_t i3 = (size_t)quads[i + 3];

    uint32 new_vindex = (uint32)new_vertices.size();

    new_vertices.push_back(vertices[i0]);
    new_vertices.push_back(vertices[i1]);
    new_vertices.push_back(vertices[i2]);
    new_vertices.push_back(vertices[i3]);

    if (!normals.empty())
    {
      new_normals.push_back(normals[i0]);
      new_normals.push_back(normals[i1]);
      new_normals.push_back(normals[i2]);
      new_normals.push_back(normals[i3]);
    }

    if (!colors.empty())
    {
      new_colors.push_back(colors[i0]);
      new_colors.push_back(colors[i1]);
      new_colors.push_back(colors[i2]);
      new_colors.push_back(colors[i3]);
    }

    if (!texcoords.empty())
    {
      new_texcoords.push_back(texcoords[i0]);
      new_texcoords.push_back(texcoords[i1]);
      new_texcoords.push_back(texcoords[i2]);
      new_texcoords.push_back(texcoords[i3]);
    }

    quads[i    ]  =  new_vindex;
    quads[i + 1]  =  new_vindex + 1;
    quads[i + 2]  =  new_vindex + 2;
    quads[i + 3]  =  new_vindex + 3;
  }

  vertices   =  new_vertices;
  normals    =  new_normals;
  colors     =  new_colors;
  texcoords  =  new_texcoords;

  invalidateGPUBuffers(BufferID::ALL);
}

void
DisplayMesh::updateBounds()
{
  if (valid_bounds) return;

  bounds = AxisAlignedBox3();
  for (size_t i = 0; i < vertices.size(); ++i)
    bounds.merge(vertices[i]);

  valid_bounds = true;
}

void
DisplayMesh::uploadToGraphicsSystem(RenderSystem & render_system)
{
  if (changed_buffers == 0) return;

  if (changed_buffers == BufferID::ALL)
  {
    if (var_area) var_area->reset();
    vertices_var = normals_var = colors_var = texcoords_var = tris_var = quads_var = edges_var = NULL;

    if (vertices.empty() || (tris.empty() && quads.empty()))
    {
      if (var_area)
      {
        render_system.destroyVARArea(var_area);
        var_area = NULL;
      }

      allGPUBuffersAreValid();
      return;
    }

    static int const PADDING = 32;
    long vertex_bytes    =  !vertices.empty()  ?  3 * 4 * (long)vertices.size()   +  PADDING : 0;  // 3 * float
    long normal_bytes    =  hasNormals()       ?  3 * 4 * (long)normals.size()    +  PADDING : 0;  // 3 * float
    long color_bytes     =  hasColors()        ?  4 * 4 * (long)colors.size()     +  PADDING : 0;  // 4 * float
    long texcoord_bytes  =  hasTexCoords()     ?  2 * 4 * (long)texcoords.size()  +  PADDING : 0;  // 2 * float

    updateEdges();

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    long num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + PADDING;
#else
    long tri_bytes   =  !tris.empty()   ?  4 * (long)tris.size()   +  PADDING : 0;  // uint32
    long quad_bytes  =  !quads.empty()  ?  4 * (long)quads.size()  +  PADDING : 0;  // uint32
    long edge_bytes  =  !edges.empty()  ?  4 * (long)edges.size()  +  PADDING : 0;  // uint32

    long num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + tri_bytes + quad_bytes + edge_bytes + PADDING;
#endif

    if (var_area)
    {
      if (var_area->getCapacity() <= num_bytes || var_area->getCapacity() > (long)(1.5 * num_bytes))
      {
        render_system.destroyVARArea(var_area);

        std::string vararea_name = getNameStr() + " VAR area";
        var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, VARArea::Usage::WRITE_OCCASIONALLY, true);
        if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
      }
      // Else no need to reset var_area, we've done it above
    }
    else
    {
      std::string vararea_name = getNameStr() + " VAR area";
      var_area = render_system.createVARArea(vararea_name.c_str(), num_bytes, VARArea::Usage::WRITE_OCCASIONALLY, true);
      if (!var_area) throw Error(getNameStr() + ": Couldn't create VAR area");
    }

    if (!vertices.empty())
    {
      vertices_var = var_area->createArray(vertex_bytes);
      if (!vertices_var) throw Error(getNameStr() + ": Couldn't create vertices VAR");
    }

    if (hasNormals())
    {
      normals_var = var_area->createArray(normal_bytes);
      if (!normals_var) throw Error(getNameStr() + ": Couldn't create normals VAR");
    }

    if (hasColors())
    {
      colors_var = var_area->createArray(color_bytes);
      if (!colors_var) throw Error(getNameStr() + ": Couldn't create colors VAR");
    }

    if (hasTexCoords())
    {
      texcoords_var = var_area->createArray(texcoord_bytes);
      if (!texcoords_var) throw Error(getNameStr() + ": Couldn't create texcoords VAR");
    }

#ifndef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    if (!tris.empty())
    {
      tris_var = var_area->createArray(tri_bytes);
      if (!tris_var) throw Error(getNameStr() + ": Couldn't create triangle indices VAR");
    }

    if (!quads.empty())
    {
      quads_var = var_area->createArray(quad_bytes);
      if (!quads_var) throw Error(getNameStr() + ": Couldn't create quad indices VAR");
    }

    if (!edges.empty())
    {
      edges_var = var_area->createArray(edge_bytes);
      if (!edges_var) throw Error(getNameStr() + ": Couldn't create edge indices VAR");
    }

    if (!tris.empty())  tris_var->updateIndices (0, (long)tris.size(),  &tris[0]);
    if (!quads.empty()) quads_var->updateIndices(0, (long)quads.size(), &quads[0]);
    if (!edges.empty()) edges_var->updateIndices(0, (long)edges.size(), &edges[0]);
#endif

    if (!vertices.empty()) vertices_var->updateVectors (0, (long)vertices.size(),  &vertices[0]);
    if (hasNormals())      normals_var->updateVectors  (0, (long)normals.size(),   &normals[0]);
    if (hasColors())       colors_var->updateColors    (0, (long)colors.size(),    &colors[0]);
    if (hasTexCoords())    texcoords_var->updateVectors(0, (long)texcoords.size(), &texcoords[0]);
  }
  else
  {
    if (!gpuBufferIsValid(BufferID::VERTEX) && !vertices.empty())
      vertices_var->updateVectors(0, (long)vertices.size(), &vertices[0]);

    if (!gpuBufferIsValid(BufferID::NORMAL) && hasNormals())
      normals_var->updateVectors (0, (long)normals.size(), &normals[0]);

    if (!gpuBufferIsValid(BufferID::COLOR) && hasColors())
      colors_var->updateColors(0, (long)colors.size(), &colors[0]);

    if (!gpuBufferIsValid(BufferID::TEXCOORD) && hasTexCoords())
      texcoords_var->updateVectors(0, (long)texcoords.size(), &texcoords[0]);
  }

  allGPUBuffersAreValid();
}

void
DisplayMesh::draw(RenderSystem & render_system, RenderOptions const & options) const
{
  if (options.drawEdges() && !wireframe_enabled)
    throw Error(getNameStr() + ": Can't draw mesh edges when wireframe mode is disabled");

  const_cast<DisplayMesh *>(this)->uploadToGraphicsSystem(render_system);

  if (!vertices_var) return;
  if (!options.drawFaces() && !options.drawEdges()) return;
  if (!options.drawFaces() && !edges_var) return;
  if (!options.drawEdges() && !tris_var && !quads_var) return;

  render_system.beginIndexedPrimitives();

    render_system.setVertexArray(vertices_var);
    if (options.sendNormals()    &&  normals_var)    render_system.setNormalArray(normals_var);
    if (options.sendColors()     &&  colors_var)     render_system.setColorArray(colors_var);
    if (options.sendTexCoords()  &&  texcoords_var)  render_system.setTexCoordArray(0, texcoords_var);

    if (options.drawFaces())
    {
      if (options.drawEdges())
      {
        render_system.pushShapeFlags();
        render_system.setPolygonOffset(true, 2);
      }

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
        if (!tris.empty()) render_system.sendIndices(RenderSystem::Primitive::TRIANGLES, (long)tris.size(), &tris[0]);
        if (!quads.empty()) render_system.sendIndices(RenderSystem::Primitive::QUADS, (long)quads.size(), &quads[0]);
#else
        if (!tris.empty())
        {
          render_system.setIndexArray(tris_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::TRIANGLES, 0, (long)tris.size());
        }

        if (!quads.empty())
        {
          render_system.setIndexArray(quads_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::QUADS, 0, (long)quads.size());
        }
#endif

      if (options.drawEdges())
        render_system.popShapeFlags();
    }

    if (options.drawEdges())
    {
      render_system.pushShader();
      render_system.pushColorFlags();
      render_system.pushTextures();

        render_system.setShader(NULL);
        render_system.setColorArray(NULL);
        render_system.setTexCoordArray(0, NULL);
        render_system.setNormalArray(NULL);
        render_system.setColor(options.edgeColor());  // set default edge color (TODO: handle per-edge colors)
        render_system.setTexture(0, NULL);

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
        if (!edges.empty()) render_system.sendIndices(RenderSystem::Primitive::LINES, (long)edges.size(), &edges[0]);
#else
        if (!edges.empty())
        {
          render_system.setIndexArray(edges_var);
          render_system.sendIndicesFromArray(RenderSystem::Primitive::LINES, 0, (long)edges.size());
        }
#endif

      render_system.popTextures();
      render_system.popColorFlags();
      render_system.popShader();
    }

  render_system.endIndexedPrimitives();
}

} // namespace Graphics
} // namespace Thea
