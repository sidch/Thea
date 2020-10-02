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

#include "DisplayMesh.hpp"
#include "../Polygon3.hpp"
#include "../UnorderedSet.hpp"

namespace Thea {
namespace Graphics {

DisplayMesh::DisplayMesh(std::string const & name)
: NamedObject(name),
  valid_bounds(true),
  wireframe_enabled(false),
  changed_buffers(BufferId::ALL),
  buf_pool(nullptr),
  vertices_buf(nullptr),
  tris_buf(nullptr),
  quads_buf(nullptr),
  normals_buf(nullptr),
  colors_buf(nullptr),
  texcoords_buf(nullptr),
  vertex_matrix(nullptr, 3, 0),
  tri_matrix(nullptr, 3, 0),
  quad_matrix(nullptr, 4, 0),
  vertex_wrapper(&vertex_matrix),
  tri_wrapper(&tri_matrix),
  quad_wrapper(&quad_matrix)
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
  changed_buffers(BufferId::ALL),
  buf_pool(nullptr),
  vertices_buf(nullptr),
  tris_buf(nullptr),
  quads_buf(nullptr),
  normals_buf(nullptr),
  colors_buf(nullptr),
  texcoords_buf(nullptr),
  vertex_matrix(nullptr, 3, 0),
  tri_matrix(nullptr, 3, 0),
  quad_matrix(nullptr, 4, 0),
  vertex_wrapper(&vertex_matrix),
  tri_wrapper(&tri_matrix),
  quad_wrapper(&quad_matrix)
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

  invalidateGpuBuffers();
}

IDenseMatrix<Real> const *
DisplayMesh::getVertexMatrix() const
{
  // Assume Vector3 is tightly packed and has no padding
  Vector3 const * buf = (vertices.empty() ? nullptr : &vertices[0]);
  new (&vertex_matrix) VertexMatrix(reinterpret_cast<Real *>(const_cast<Vector3 *>(buf)), 3, numVertices());
  return &vertex_wrapper;
}

IDenseMatrix<uint32> const *
DisplayMesh::getTriangleMatrix() const
{
  uint32 const * buf = (tris.empty() ? nullptr : &tris[0]);
  new (&tri_matrix) TriangleMatrix(const_cast<uint32 *>(buf), 3, numTriangles());
  return &tri_wrapper;
}

IDenseMatrix<uint32> const *
DisplayMesh::getQuadMatrix() const
{
  uint32 const * buf = (quads.empty() ? nullptr : &quads[0]);
  new (&quad_matrix) QuadMatrix(const_cast<uint32 *>(buf), 4, numQuads());
  return &quad_wrapper;
}

DisplayMesh::Vertex
DisplayMesh::getVertex(intx i)
{
  debugAssertM(i >= 0 && i < (intx)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  size_t si = (size_t)i;
  return Vertex(this, vertices[si],
                hasNormals()    ?  &normals[si]    :  nullptr,
                hasColors()     ?  &colors[si]     :  nullptr,
                hasTexCoords()  ?  &texcoords[si]  :  nullptr);
}

DisplayMesh::IndexTriple
DisplayMesh::getTriangle(intx tri_index) const
{
  debugAssertM(tri_index >= 0 && 3 * tri_index < (intx)tris.size(), getNameStr() + ": Triangle index out of bounds");

  size_t base_index = (size_t)(3 * tri_index);
  IndexTriple tri;
  tri[0] = (intx)tris[base_index];
  tri[1] = (intx)tris[base_index + 1];
  tri[2] = (intx)tris[base_index + 2];

  return tri;
}

DisplayMesh::IndexQuad
DisplayMesh::getQuad(intx quad_index) const
{
  debugAssertM(quad_index >= 0 && 4 * quad_index < (intx)quads.size(), getNameStr() + ": Quad index out of bounds");

  size_t base_index = (size_t)(4 * quad_index);
  IndexQuad quad;
  quad[0] = (intx)quads[base_index];
  quad[1] = (intx)quads[base_index + 1];
  quad[2] = (intx)quads[base_index + 2];
  quad[3] = (intx)quads[base_index + 3];

  return quad;
}

void
DisplayMesh::addColors()
{
  if (colors.size() < vertices.size())
  {
    colors.resize(vertices.size(), ColorRgba(0, 0, 0, 0));
    invalidateGpuBuffers();
  }
}

void
DisplayMesh::addNormals()
{
  if (normals.size() < vertices.size())
  {
    normals.resize(vertices.size(), Vector3::Zero());
    invalidateGpuBuffers();
  }
}

void
DisplayMesh::addTexCoords()
{
  if (texcoords.size() < vertices.size())
  {
    texcoords.resize(vertices.size(), Vector2::Zero());
    invalidateGpuBuffers();
  }
}

intx
DisplayMesh::addVertex(Vector3 const & point, intx source_index, Vector3 const * normal, ColorRgba const * color,
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

  intx index = (intx)vertices.size();

  if (valid_bounds)
    bounds.merge(point);

  vertices.push_back(point);
  if (source_index >= 0)  vertex_source_indices.push_back(source_index);
  if (normal)             normals.push_back(*normal);
  if (color)              colors.push_back(*color);
  if (texcoord)           texcoords.push_back(*texcoord);

  invalidateGpuBuffers();

  return index;
}

intx
DisplayMesh::addTriangle(intx vi0, intx vi1, intx vi2, intx source_face_index)
{
  debugAssertM(vi0 >= 0 && vi1 >= 0 && vi2 >= 0
            && vi0 < (intx)vertices.size()
            && vi1 < (intx)vertices.size()
            && vi2 < (intx)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  alwaysAssertM((source_face_index >= 0 && 3 * tri_source_face_indices.size() == tris.size())
             || (source_face_index < 0 && tri_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

  intx index = (intx)(tris.size() / 3);

  tris.push_back((uint32)vi0);
  tris.push_back((uint32)vi1);
  tris.push_back((uint32)vi2);

  if (source_face_index >= 0)
    tri_source_face_indices.push_back(source_face_index);

  invalidateGpuBuffers();

  return index;
}

intx
DisplayMesh::addQuad(intx vi0, intx vi1, intx vi2, intx vi3, intx source_face_index)
{
  debugAssertM(vi0 >= 0 && vi1 >= 0 && vi2 >= 0 && vi3 >= 0
            && vi0 < (intx)vertices.size()
            && vi1 < (intx)vertices.size()
            && vi2 < (intx)vertices.size()
            && vi3 < (intx)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  alwaysAssertM((source_face_index >= 0 && 4 * quad_source_face_indices.size() == quads.size())
             || (source_face_index < 0 && quad_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no quad face source indices");

  intx index = (intx)(quads.size() / 4);

  quads.push_back((uint32)vi0);
  quads.push_back((uint32)vi1);
  quads.push_back((uint32)vi2);
  quads.push_back((uint32)vi3);

  if (source_face_index >= 0)
    quad_source_face_indices.push_back(source_face_index);

  invalidateGpuBuffers();

  return index;
}

DisplayMesh::Face
DisplayMesh::addFace(int num_vertices, intx const * face_vertex_indices_, intx source_face_index)
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
    intx vi = face_vertex_indices_[i];
    debugAssertM(vi >= 0 && vi < (intx)vertices.size(), getName() + format(": Vertex index %ld out of bounds", vi));

    poly.addVertex(vertices[vi], vi);
  }

  intx num_tris = poly.triangulate(triangulated_indices);
  if (num_tris <= 0)
    return Face();

  // debugAssertM(num_tris == 3 * (num_vertices - 2),
  //              getName() + format(": Triangulation of polygonal face yielded %l triangles, whereas %l were expected",
  //                                 num_tris, num_vertices - 2));

  alwaysAssertM((source_face_index >= 0 && 3 * tri_source_face_indices.size() == tris.size())
             || (source_face_index < 0 && tri_source_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

  intx starting_index = numTriangles();
  size_t num_tri_verts = (size_t)(3 * num_tris);
  for (size_t i = 0; i < num_tri_verts; i += 3)
  {
    tris.push_back((uint32)triangulated_indices[i    ]);
    tris.push_back((uint32)triangulated_indices[i + 1]);
    tris.push_back((uint32)triangulated_indices[i + 2]);
    if (source_face_index >= 0) tri_source_face_indices.push_back(source_face_index);
  }

  invalidateGpuBuffers();

  return Face(this, num_vertices, true, starting_index, num_tris);
}

void
DisplayMesh::removeTriangle(intx tri_index)
{
  debugAssertM(tri_index >= 0 && 3 * tri_index < (intx)tris.size(), getNameStr() + ": Triangle index out of bounds");

  IndexArray::iterator ii = tris.begin() + 3 * tri_index;
  tris.erase(ii, ii + 3);

  invalidateGpuBuffers();
}

void
DisplayMesh::removeTriangles(intx begin, intx num_triangles)
{
  debugAssertM(begin >= 0 && 3 * begin < (intx)tris.size(), getNameStr() + ": Triangle index out of bounds");

  IndexArray::iterator ii = tris.begin() + 3 * begin;
  tris.erase(ii, ii + 3 * num_triangles);

  invalidateGpuBuffers();
}

void
DisplayMesh::removeQuad(intx quad_index)
{
  debugAssertM(quad_index >= 0 && 4 * quad_index < (intx)quads.size(), getNameStr() + ": Quad index out of bounds");

  IndexArray::iterator ii = quads.begin() + 4 * quad_index;
  quads.erase(ii, ii + 4);

  invalidateGpuBuffers();
}

void
DisplayMesh::removeQuads(intx begin, intx num_quads)
{
  debugAssertM(begin >= 0 && 4 * begin < (intx)quads.size(), getNameStr() + ": Quad index out of bounds");

  IndexArray::iterator ii = quads.begin() + 4 * begin;
  quads.erase(ii, ii + 4 * num_quads);

  invalidateGpuBuffers();
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
    normals[i] = Vector3::Zero();

  // TODO: weight normals by face area?
  Vector3 n;
  uint32 i0, i1, i2, i3;
  for (size_t i = 0; i < tris.size(); i += 3)
  {
    i0 = tris[i], i1 = tris[i + 1], i2 = tris[i + 2];
    Vector3 const & v0 = vertices[i0];
    Vector3 const & v1 = vertices[i1];
    Vector3 const & v2 = vertices[i2];

    n = (v2 - v1).cross(v0 - v1).normalized();
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

    n = (v2 - v1).cross(v0 - v1).normalized();
    normals[i0] += n;
    normals[i1] += n;
    normals[i2] += n;
    normals[i3] += n;
  }

  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = normals[i].normalized();

  invalidateGpuBuffers(topo_change ? BufferId::ALL : BufferId::NORMAL);
}

void
DisplayMesh::flipNormals()
{
  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = -normals[i];

  invalidateGpuBuffers(BufferId::NORMAL);
}

void
DisplayMesh::updateEdges()
{
  edges.clear();

  if (wireframe_enabled)
  {
    typedef std::pair<uint32, uint32> Edge;
    typedef UnorderedSet<Edge> EdgeSet;

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

  invalidateGpuBuffers(BufferId::ALL);
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
DisplayMesh::uploadToGraphicsSystem(IRenderSystem & render_system)
{
  if (changed_buffers == 0) return;

  if (changed_buffers == BufferId::ALL)
  {
    if (buf_pool) buf_pool->reset();
    vertices_buf = normals_buf = colors_buf = texcoords_buf = tris_buf = quads_buf = edges_buf = nullptr;

    if (vertices.empty() || (tris.empty() && quads.empty()))
    {
      if (buf_pool)
      {
        render_system.destroyBufferPool(buf_pool);
        buf_pool = nullptr;
      }

      allGpuBuffersAreValid();
      return;
    }

    static int const PADDING = 32;
    intx vertex_bytes    =  !vertices.empty()  ?  3 * 4 * (intx)vertices.size()   +  PADDING : 0;  // 3 * float
    intx normal_bytes    =  hasNormals()       ?  3 * 4 * (intx)normals.size()    +  PADDING : 0;  // 3 * float
    intx color_bytes     =  hasColors()        ?  4 * 4 * (intx)colors.size()     +  PADDING : 0;  // 4 * float
    intx texcoord_bytes  =  hasTexCoords()     ?  2 * 4 * (intx)texcoords.size()  +  PADDING : 0;  // 2 * float

    updateEdges();

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    intx num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + PADDING;
#else
    intx tri_bytes   =  !tris.empty()   ?  4 * (intx)tris.size()   +  PADDING : 0;  // uint32
    intx quad_bytes  =  !quads.empty()  ?  4 * (intx)quads.size()  +  PADDING : 0;  // uint32
    intx edge_bytes  =  !edges.empty()  ?  4 * (intx)edges.size()  +  PADDING : 0;  // uint32

    intx num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + tri_bytes + quad_bytes + edge_bytes + PADDING;
#endif

    if (buf_pool)
    {
      if (buf_pool->getCapacity() <= num_bytes || buf_pool->getCapacity() > (intx)(1.5 * num_bytes))
      {
        render_system.destroyBufferPool(buf_pool);

        std::string pool_name = getNameStr() + " buffer pool";
        buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
        if (!buf_pool) throw Error(getNameStr() + ": Couldn't create buffer pool");
      }
      // Else no need to reset buf_pool, we've done it above
    }
    else
    {
      std::string pool_name = getNameStr() + " buffer pool";
      buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
      if (!buf_pool) throw Error(getNameStr() + ": Couldn't create buffer pool");
    }

    if (!vertices.empty())
    {
      vertices_buf = buf_pool->createBuffer(vertex_bytes);
      if (!vertices_buf) throw Error(getNameStr() + ": Couldn't create vertices buffer");
    }

    if (hasNormals())
    {
      normals_buf = buf_pool->createBuffer(normal_bytes);
      if (!normals_buf) throw Error(getNameStr() + ": Couldn't create normals buffer");
    }

    if (hasColors())
    {
      colors_buf = buf_pool->createBuffer(color_bytes);
      if (!colors_buf) throw Error(getNameStr() + ": Couldn't create colors buffer");
    }

    if (hasTexCoords())
    {
      texcoords_buf = buf_pool->createBuffer(texcoord_bytes);
      if (!texcoords_buf) throw Error(getNameStr() + ": Couldn't create texcoords buffer");
    }

#ifndef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    if (!tris.empty())
    {
      tris_buf = buf_pool->createBuffer(tri_bytes);
      if (!tris_buf) throw Error(getNameStr() + ": Couldn't create triangle indices buffer");
    }

    if (!quads.empty())
    {
      quads_buf = buf_pool->createBuffer(quad_bytes);
      if (!quads_buf) throw Error(getNameStr() + ": Couldn't create quad indices buffer");
    }

    if (!edges.empty())
    {
      edges_buf = buf_pool->createBuffer(edge_bytes);
      if (!edges_buf) throw Error(getNameStr() + ": Couldn't create edge indices buffer");
    }

    if (!tris.empty())  tris_buf->updateIndices (0, (int64)tris.size(),  NumericType::UINT32, &tris[0]);
    if (!quads.empty()) quads_buf->updateIndices(0, (int64)quads.size(), NumericType::UINT32, &quads[0]);
    if (!edges.empty()) edges_buf->updateIndices(0, (int64)edges.size(), NumericType::UINT32, &edges[0]);
#endif

    if (!vertices.empty()) vertices_buf->updateAttributes (0, (int64)vertices.size(),  3, NumericType::REAL, &vertices[0]);
    if (hasNormals())      normals_buf->updateAttributes  (0, (int64)normals.size(),   3, NumericType::REAL, &normals[0]);
    if (hasColors())       colors_buf->updateAttributes   (0, (int64)colors.size(),    4, NumericType::REAL, &colors[0]);
    if (hasTexCoords())    texcoords_buf->updateAttributes(0, (int64)texcoords.size(), 2, NumericType::REAL, &texcoords[0]);
  }
  else
  {
    if (!gpuBufferIsValid(BufferId::VERTEX) && !vertices.empty())
      vertices_buf->updateAttributes(0, (int64)vertices.size(), 3, NumericType::REAL, &vertices[0]);

    if (!gpuBufferIsValid(BufferId::NORMAL) && hasNormals())
      normals_buf->updateAttributes(0, (int64)normals.size(), 3, NumericType::REAL, &normals[0]);

    if (!gpuBufferIsValid(BufferId::COLOR) && hasColors())
      colors_buf->updateAttributes(0, (int64)colors.size(), 4, NumericType::REAL, &colors[0]);

    if (!gpuBufferIsValid(BufferId::TEXCOORD) && hasTexCoords())
      texcoords_buf->updateAttributes(0, (int64)texcoords.size(), 2, NumericType::REAL, &texcoords[0]);
  }

  allGpuBuffersAreValid();
}

void
DisplayMesh::draw(IRenderSystem * render_system, IRenderOptions const * options) const
{
  if (!render_system) { THEA_ERROR << getName() << ": Can't display mesh on a null rendersystem"; return; }
  if (!options) options = RenderOptions::defaults();

  if (options->drawEdges() && !wireframe_enabled)
    throw Error(getNameStr() + ": Can't draw mesh edges when wireframe mode is disabled");

  const_cast<DisplayMesh *>(this)->uploadToGraphicsSystem(*render_system);

  if (!vertices_buf) return;
  if (!options->drawFaces() && !options->drawEdges()) return;
  if (!options->drawFaces() && !edges_buf) return;
  if (!options->drawEdges() && !tris_buf && !quads_buf) return;

  render_system->beginIndexedPrimitives();

    render_system->setVertexBuffer(vertices_buf);
    if (options->sendNormals()    &&  normals_buf)    render_system->setNormalBuffer(normals_buf);
    if (options->sendColors()     &&  colors_buf)     render_system->setColorBuffer(colors_buf);
    if (options->sendTexCoords()  &&  texcoords_buf)  render_system->setTexCoordBuffer(0, texcoords_buf);

    if (options->drawFaces())
    {
      if (options->drawEdges())
      {
        render_system->pushShapeFlags();
        render_system->setPolygonOffset(true, 2);
      }

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
        if (!tris.empty()) render_system->sendIndices(IRenderSystem::Primitive::TRIANGLES, (int64)tris.size(), &tris[0]);
        if (!quads.empty()) render_system->sendIndices(IRenderSystem::Primitive::QUADS, (int64)quads.size(), &quads[0]);
#else
        if (!tris.empty())
        {
          render_system->setIndexBuffer(tris_buf);
          render_system->sendIndicesFromBuffer(IRenderSystem::Primitive::TRIANGLES, 0, (int64)tris.size());
        }

        if (!quads.empty())
        {
          render_system->setIndexBuffer(quads_buf);
          render_system->sendIndicesFromBuffer(IRenderSystem::Primitive::QUADS, 0, (int64)quads.size());
        }
#endif

      if (options->drawEdges())
        render_system->popShapeFlags();
    }

    if (options->drawEdges())
    {
      render_system->pushShader();
      render_system->pushColorFlags();
      render_system->pushTextures();

        render_system->setShader(nullptr);
        render_system->setColorBuffer(nullptr);
        render_system->setTexCoordBuffer(0, nullptr);
        render_system->setNormalBuffer(nullptr);
        render_system->setColor(options->edgeColor());  // set default edge color (TODO: handle per-edge colors)
        render_system->setTexture(0, nullptr);

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
        if (!edges.empty()) render_system->sendIndices(IRenderSystem::Primitive::LINES, (int64)edges.size(), &edges[0]);
#else
        if (!edges.empty())
        {
          render_system->setIndexBuffer(edges_buf);
          render_system->sendIndicesFromBuffer(IRenderSystem::Primitive::LINES, 0, (int64)edges.size());
        }
#endif

      render_system->popTextures();
      render_system->popColorFlags();
      render_system->popShader();
    }

  render_system->endIndexedPrimitives();
}

} // namespace Graphics
} // namespace Thea
