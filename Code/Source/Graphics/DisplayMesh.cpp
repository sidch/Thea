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
  changed_buffers(BufferId::ALL),
  buf_pool(nullptr),
  vertices_buf(nullptr),
  tris_buf(nullptr),
  normals_buf(nullptr),
  colors_buf(nullptr),
  texcoords_buf(nullptr),
  vertex_matrix(nullptr, 3, 0),
  tri_matrix(nullptr, 3, 0),
  vertex_wrapper(&vertex_matrix),
  tri_wrapper(&tri_matrix)
{
  face_starting_tris.resize(1); face_starting_tris[0] = 0;
}

DisplayMesh::DisplayMesh(DisplayMesh const & src)
: NamedObject(src),
  vertices(src.vertices),
  normals(src.normals),
  colors(src.colors),
  texcoords(src.texcoords),
  tris(src.tris),
  face_starting_tris(src.face_starting_tris),
  large_poly_verts(src.large_poly_verts),
  valid_bounds(src.valid_bounds),
  bounds(src.bounds),
  changed_buffers(BufferId::ALL),
  buf_pool(nullptr),
  vertices_buf(nullptr),
  tris_buf(nullptr),
  normals_buf(nullptr),
  colors_buf(nullptr),
  texcoords_buf(nullptr),
  vertex_matrix(nullptr, 3, 0),
  tri_matrix(nullptr, 3, 0),
  vertex_wrapper(&vertex_matrix),
  tri_wrapper(&tri_matrix)
{}

void
DisplayMesh::clear()
{
  vertices.clear();
  normals.clear();
  colors.clear();
  texcoords.clear();
  tris.clear();
  edges.clear();

  vertex_src_indices.clear();
  tri_src_face_indices.clear();
  face_starting_tris.resize(1); face_starting_tris[0] = 0;
  large_poly_verts.clear();

  face_vertex_indices.clear();
  triangulated_indices.clear();

  valid_bounds = true;
  bounds = AxisAlignedBox3();

  invalidateGpuBuffers();
}

IDenseMatrix<Real> const *
DisplayMesh::getVertexMatrix() const
{
  return const_cast<DisplayMesh *>(this)->getVertexMatrix();
}

IDenseMatrix<uint32> const *
DisplayMesh::getTriangleMatrix() const
{
  return const_cast<DisplayMesh *>(this)->getTriangleMatrix();
}

IDenseMatrix<uint32> const *
DisplayMesh::getQuadMatrix() const
{
  // FIXME: Even if the mesh has quads, we pack everything as triangles (for compatibility with modern OpenGL) so we
  // currently have no efficient way of returning the quads.
  typedef Matrix<4, 0, uint32, MatrixLayout::COLUMN_MAJOR> EmptyQuadMatrix;
  static EmptyQuadMatrix EMPTY_QUAD_MATRIX;
  static MatrixWrapper<EmptyQuadMatrix> EMPTY_QUAD_WRAPPER(&EMPTY_QUAD_MATRIX);
  return &EMPTY_QUAD_WRAPPER;
}

IDenseMatrix<Real> *
DisplayMesh::getVertexMatrix()
{
  // Assume Vector3 is tightly packed and has no padding
  Vector3 * buf = (vertices.empty() ? nullptr : vertices.data());
  new (&vertex_matrix) VertexMatrix(reinterpret_cast<Real *>(buf), 3, numVertices());
  return &vertex_wrapper;
}

IDenseMatrix<uint32> *
DisplayMesh::getTriangleMatrix()
{
  uint32 * buf = (tris.empty() ? nullptr : tris.data());
  new (&tri_matrix) TriangleMatrix(buf, 3, numTriangles());
  return &tri_wrapper;
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
DisplayMesh::addVertex(Vector3 const & point, intx src_index, Vector3 const * normal, ColorRgba const * color,
                       Vector2 const * texcoord)
{
  alwaysAssertM((src_index >= 0 && vertex_src_indices.size() == vertices.size())
             || (src_index < 0 && vertex_src_indices.empty()),
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
  if (src_index >= 0)  vertex_src_indices.push_back(src_index);
  if (normal)          normals.push_back(*normal);
  if (color)           colors.push_back(*color);
  if (texcoord)        texcoords.push_back(*texcoord);

  invalidateGpuBuffers();

  return index;
}

intx
DisplayMesh::addTriangle(intx vi0, intx vi1, intx vi2, intx src_face_index)
{
  debugAssertM(vi0 >= 0 && vi1 >= 0 && vi2 >= 0
            && vi0 < (intx)vertices.size()
            && vi1 < (intx)vertices.size()
            && vi2 < (intx)vertices.size(), getNameStr() + ": Vertex index out of bounds");

  alwaysAssertM((src_face_index >= 0 && 3 * tri_src_face_indices.size() == tris.size())
             || (src_face_index < 0 && tri_src_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

  debugAssertM(face_starting_tris.back() == numTriangles(), getNameStr() + ": Last entry of face_starting_tris is invalid");

  tris.push_back((uint32)vi0);
  tris.push_back((uint32)vi1);
  tris.push_back((uint32)vi2);

  intx face_index = numFaces();
  face_starting_tris.push_back(numTriangles());  // old final entry of face_starting_tris is index of first new triangle

  if (src_face_index >= 0)
    tri_src_face_indices.push_back(src_face_index);

  addEdge((uint32)vi0, (uint32)vi1);
  addEdge((uint32)vi1, (uint32)vi2);
  addEdge((uint32)vi2, (uint32)vi0);

  invalidateGpuBuffers(BufferId::TOPOLOGY);

  return face_index;
}

intx
DisplayMesh::addFace(int num_vertices, intx const * vertex_indices, intx src_face_index)
{
  if (num_vertices < 3)
  {
    THEA_DEBUG << getName() << ": Skipping face -- too few vertices (" << num_vertices << ')';
    return -1;
  }

  if (num_vertices == 3)
    return addTriangle(vertex_indices[0], vertex_indices[1], vertex_indices[2], src_face_index);

#ifdef THEA_DEBUG
  for (int i = 0; i < num_vertices; ++i)
  {
    intx vi = vertex_indices[i];
    debugAssertM(vi >= 0 && vi < (intx)vertices.size(), getName() + format(": Vertex index %ld out of bounds", vi));
  }
#endif

  intx num_tris = 0;
  if (num_vertices == 4)
  {
    triangulated_indices.resize(6);  // std::vector spec minimizes reallocs by never shrinking storage
    num_tris = Polygon3::triangulateQuad(vertices[(size_t)vertex_indices[0]],
                                         vertices[(size_t)vertex_indices[1]],
                                         vertices[(size_t)vertex_indices[2]],
                                         vertices[(size_t)vertex_indices[3]],
                                         triangulated_indices[0], triangulated_indices[1], triangulated_indices[2],
                                         triangulated_indices[3], triangulated_indices[4], triangulated_indices[5]);
  }
  else
  {
    Polygon3 poly;
    for (int i = 0; i < num_vertices; ++i)
      poly.addVertex(vertices[(size_t)vertex_indices[i]]);

    num_tris = poly.triangulate(triangulated_indices);
  }

  if (num_tris <= 0)
  {
    THEA_DEBUG << getName() << ": Skipping degenerate face -- no triangles produced";
    return -1;
  }

  // debugAssertM(num_tris == 3 * (num_vertices - 2),
  //              getName() + format(": Triangulation of polygonal face yielded %l triangles, whereas %l were expected",
  //                                 num_tris, num_vertices - 2));

  alwaysAssertM((src_face_index >= 0 && 3 * tri_src_face_indices.size() == tris.size())
             || (src_face_index < 0 && tri_src_face_indices.empty()),
                getNameStr() + ": Mesh must have all or no triangle face source indices");

#ifdef THEA_DEBUG
  intx first_tri = numTriangles();
  debugAssertM(face_starting_tris.back() == first_tri, getNameStr() + ": Last entry of face_starting_tris is invalid");
#endif

  size_t num_tri_verts = (size_t)(3 * num_tris);
  for (size_t i = 0; i < num_tri_verts; i += 3)
  {
    tris.push_back((uint32)vertex_indices[triangulated_indices[i    ]]);
    tris.push_back((uint32)vertex_indices[triangulated_indices[i + 1]]);
    tris.push_back((uint32)vertex_indices[triangulated_indices[i + 2]]);
    if (src_face_index >= 0) tri_src_face_indices.push_back(src_face_index);
  }

  intx face_index = numFaces();
  face_starting_tris.push_back(numTriangles());  // old final entry of face_starting_tris is index of first new triangle

  if (num_vertices >= 5 && num_tris >= 2)  // bona fide higher-degree poly which does not reduce to a single triangle
    large_poly_verts[first_tri].assign(vertex_indices, vertex_indices + num_vertices);

  // Add edges of original polygon
  {
    uint32 i = (uint32)(num_vertices - 1);
    for (uint32 j = 0; j < (uint32)num_vertices; ++j)
    {
      addEdge(i, j);
      i = j;
    }
  }

  invalidateGpuBuffers(BufferId::TOPOLOGY);

  return face_index;
}

DisplayMesh::Face
DisplayMesh::getFace(intx face_index)
{
  if (face_index < 0 || face_index + 1 >= (intx)face_starting_tris.size())
    return Face();

  auto begin = face_starting_tris[(size_t)face_index], end = face_starting_tris[(size_t)face_index + 1];
  auto num_tris = end - begin;

  switch (num_tris)
  {
    case 1: return Face(this, 3, begin, 1);  // simple triangle

    case 2:  // actual quad or larger degenerate poly which only produced two valid triangles
    {
      auto loc = large_poly_verts.find(begin);
      return Face(this, (loc == large_poly_verts.end() ? 4 : (int)loc->second.size()), begin, 2);
    }

    case 0: return Face();  // totally degenerate, produced no triangles

    default:  // higher-degree polygon
    {
      auto loc = large_poly_verts.find(begin);
      debugAssertM(loc != large_poly_verts.end(), getNameStr() + format(": Vertices of face %ld not found", (long)face_index));
      return Face(this, (int)loc->second.size(), begin, num_tris);
    }
  }
}

void
DisplayMesh::computeAveragedVertexNormals()
{
  bool topo_change = (normals.size() != vertices.size());

  normals.resize(vertices.size());
  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = Vector3::Zero();

  // TODO: weight normals by face area?
  for (size_t i = 0; i < tris.size(); i += 3)
  {
    uint32 i0 = tris[i], i1 = tris[i + 1], i2 = tris[i + 2];
    Vector3 const & v0 = vertices[i0];
    Vector3 const & v1 = vertices[i1];
    Vector3 const & v2 = vertices[i2];

    Vector3 n = (v2 - v1).cross(v0 - v1).normalized();
    normals[i0] += n;
    normals[i1] += n;
    normals[i2] += n;
  }

  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = normals[i].normalized();

  invalidateGpuBuffers((topo_change ? BufferId::TOPOLOGY : 0) | BufferId::VERTEX_NORMAL);
}

void
DisplayMesh::flipNormals()
{
  for (size_t i = 0; i < normals.size(); ++i)
    normals[i] = -normals[i];

  invalidateGpuBuffers(BufferId::VERTEX_NORMAL);
}

void
DisplayMesh::recomputeEdges()
{
  edges.clear();

  intx num_faces = numFaces();
  for (intx i = 0; i < num_faces; ++i)
  {
    if (isTriangle(i))
    {
      intx tri_base = 3 * face_starting_tris[(size_t)i];
      for (uint32 j = 0; j < 3; ++j)
        addEdge(tris[tri_base + j], tris[tri_base + (j + 1) % 3]);
    }
    else
    {
      auto f = getFace(i);
      if (!f) continue;

      auto nv = f.numVertices();
      if (nv <= 0) continue;

      uint32 ei0 = f.getVertexIndex(nv - 1);
      for (int j = 0; j < nv; ++j)
      {
        uint32 ei1 = f.getVertexIndex(j);
        addEdge(ei0, ei1);
        ei0 = ei1;
      }
    }
  }

  invalidateGpuBuffers(BufferId::TOPOLOGY);

  THEA_DEBUG << getName() << ": Mesh has " << edges.size() / 2 << " edges";
}

void
DisplayMesh::packEdges() const
{
  packed_edges.resize(2 * edges.size());
  size_t i = 0;
  for (auto ei = edges.begin(); ei != edges.end(); ++ei, i += 2)
  {
    packed_edges[i    ] = ei->first;
    packed_edges[i + 1] = ei->second;
  }
}

void
DisplayMesh::isolateTriangles()
{
  // Note: does not change face_starting_tris

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

  vertices   =  new_vertices;
  normals    =  new_normals;
  colors     =  new_colors;
  texcoords  =  new_texcoords;

  recomputeEdges();

  invalidateGpuBuffers();
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

bool
DisplayMesh::uploadToGraphicsSystem(IRenderSystem & render_system, IRenderOptions const & options)
{
  if (changed_buffers == 0) return true;

  if (!isGpuBufferValid(BufferId::TOPOLOGY) || (options.drawEdges() == packed_edges.empty() && !edges.empty()))
    invalidateGpuBuffers(BufferId::ALL);  // need to reallocate pool

  if (changed_buffers == BufferId::ALL)
  {
    if (buf_pool) buf_pool->reset();
    vertices_buf = normals_buf = colors_buf = texcoords_buf = tris_buf = edges_buf = nullptr;

    if (vertices.empty() || (tris.empty() && edges.empty()))
    {
      if (buf_pool)
      {
        render_system.destroyBufferPool(buf_pool);
        buf_pool = nullptr;
      }

      setAllGpuBuffersValid();
      return true;
    }

    static int const PADDING = 32;
    intx vertex_bytes    =  !vertices.empty()  ?  3 * 4 * (intx)vertices.size()   +  PADDING : 0;  // 3 * float
    intx normal_bytes    =  hasNormals()       ?  3 * 4 * (intx)normals.size()    +  PADDING : 0;  // 3 * float
    intx color_bytes     =  hasColors()        ?  4 * 4 * (intx)colors.size()     +  PADDING : 0;  // 4 * float
    intx texcoord_bytes  =  hasTexCoords()     ?  2 * 4 * (intx)texcoords.size()  +  PADDING : 0;  // 2 * float

    // Create the edge buffer only if we need it, else save GPU memory. We assume toggling the edge display on/off will be a
    // relatively rare occurrence.
    if (options.drawEdges())  packEdges();
    else                      packed_edges.clear();

    // After this, options.drawEdges() != packed_edges.empty() except in the trivial case of edges.empty().

#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    intx num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + PADDING;
#else
    intx tri_bytes   =  !tris.empty()          ?  4 * (intx)tris.size()          +  PADDING : 0;  // uint32
    intx edge_bytes  =  !packed_edges.empty()  ?  4 * (intx)packed_edges.size()  +  PADDING : 0;  // uint32

    intx num_bytes = vertex_bytes + normal_bytes + color_bytes + texcoord_bytes + tri_bytes + edge_bytes + PADDING;
#endif

    if (buf_pool)
    {
      if (buf_pool->getCapacity() <= num_bytes || buf_pool->getCapacity() > (intx)(1.5 * num_bytes))
      {
        render_system.destroyBufferPool(buf_pool);

        std::string pool_name = getNameStr() + " buffer pool";
        buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
        if (!buf_pool) return false;
      }
      // Else no need to reset buf_pool, we've done it above
    }
    else
    {
      std::string pool_name = getNameStr() + " buffer pool";
      buf_pool = render_system.createBufferPool(pool_name.c_str(), num_bytes, IBufferPool::Usage::WRITE_OCCASIONALLY, true);
      if (!buf_pool) return false;
    }

    if (!vertices.empty() && !(vertices_buf = buf_pool->createBuffer(vertex_bytes))   ) return false;
    if (hasNormals()      && !(normals_buf = buf_pool->createBuffer(normal_bytes))    ) return false;
    if (hasColors()       && !(colors_buf = buf_pool->createBuffer(color_bytes))      ) return false;
    if (hasTexCoords()    && !(texcoords_buf = buf_pool->createBuffer(texcoord_bytes))) return false;

#ifndef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
    if (!tris.empty()  && !(tris_buf  = buf_pool->createBuffer(tri_bytes)))  return false;
    if (!packed_edges.empty() && !(edges_buf = buf_pool->createBuffer(edge_bytes))) return false;

    if (!tris.empty() && !tris_buf->updateIndices (0, (int64)tris.size(), NumericType::UINT32, &tris[0]))
      return false;

    if (!packed_edges.empty()
     && !edges_buf->updateIndices(0, (int64)packed_edges.size(), NumericType::UINT32, &packed_edges[0]))
      return false;
#endif

    if (!vertices.empty() && !vertices_buf->updateAttributes(0, (int64)vertices.size(), 3, NumericType::REAL, &vertices[0]))
      return false;

    if (hasNormals() && !normals_buf->updateAttributes(0, (int64)normals.size(), 3, NumericType::REAL, &normals[0]))
      return false;

    if (hasColors() && !colors_buf->updateAttributes(0, (int64)colors.size(), 4, NumericType::REAL, &colors[0]))
      return false;

    if (hasTexCoords() && !texcoords_buf->updateAttributes(0, (int64)texcoords.size(), 2, NumericType::REAL, &texcoords[0]))
      return false;
  }
  else
  {
    if (!isGpuBufferValid(BufferId::VERTEX_POSITION) && !vertices.empty()
     && !vertices_buf->updateAttributes(0, (int64)vertices.size(), 3, NumericType::REAL, &vertices[0])) return false;

    if (!isGpuBufferValid(BufferId::VERTEX_NORMAL) && hasNormals()
     && !normals_buf->updateAttributes(0, (int64)normals.size(), 3, NumericType::REAL, &normals[0])) return false;

    if (!isGpuBufferValid(BufferId::VERTEX_COLOR) && hasColors()
     && !colors_buf->updateAttributes(0, (int64)colors.size(), 4, NumericType::REAL, &colors[0])) return false;

    if (!isGpuBufferValid(BufferId::VERTEX_TEXCOORD) && hasTexCoords()
     && !texcoords_buf->updateAttributes(0, (int64)texcoords.size(), 2, NumericType::REAL, &texcoords[0])) return false;
  }

  setAllGpuBuffersValid();

  return true;
}

int8
DisplayMesh::draw(IRenderSystem * render_system, IRenderOptions const * options) const
{
  if (!render_system) { THEA_ERROR << getName() << ": Can't display mesh on a null rendersystem"; return false; }
  if (!options) options = RenderOptions::defaults();

  if (!const_cast<DisplayMesh *>(this)->uploadToGraphicsSystem(*render_system, *options)) return false;

  if (!vertices_buf) return true;
  if (!options->drawFaces() && !options->drawEdges()) return true;
  if (!options->drawFaces() && edges.empty()) return true;
  if (!options->drawEdges() && tris.empty()) return true;

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
        render_system->setPolygonOffset(true, 1);
      }

        if (!tris.empty())
        {
#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
          render_system->sendIndices(IRenderSystem::Primitive::TRIANGLES, (int64)tris.size(), &tris[0]);
#else
          render_system->setIndexBuffer(tris_buf);
          render_system->sendIndicesFromBuffer(IRenderSystem::Primitive::TRIANGLES, 0, (int64)tris.size());
#endif
        }

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

        if (!packed_edges.empty())
        {
#ifdef THEA_DISPLAY_MESH_NO_INDEX_ARRAY
          render_system->sendIndices(IRenderSystem::Primitive::LINES, (int64)packed_edges.size(), &packed_edges[0]);
#else
          render_system->setIndexBuffer(edges_buf);
          render_system->sendIndicesFromBuffer(IRenderSystem::Primitive::LINES, 0, (int64)packed_edges.size());
#endif
        }

      render_system->popTextures();
      render_system->popColorFlags();
      render_system->popShader();
    }

  render_system->endIndexedPrimitives();

  if (char const * err = render_system->getLastError())
  { THEA_ERROR << getName() << ": Rendering error (" << err << ')'; return false; }

  return true;
}

} // namespace Graphics
} // namespace Thea
