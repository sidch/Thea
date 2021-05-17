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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_MeshTriangles_hpp__
#define __Thea_Algorithms_MeshTriangles_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"
#include "../Polygon3.hpp"
#include "../Triangle3.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

/**
 * A set of three vertices of a single face of a mesh. If the mesh vertices change, so does this triplet. This base template
 * works for general and DCEL meshes.
 *
 * @see GeneralMesh, DcelMesh
 */
template <typename MeshT, typename Enable = void>
class MeshVertexTriple
{
  public:
    typedef MeshT                          Mesh;                   ///< The mesh type.
    typedef typename Mesh::Vertex       *  MeshVertexHandle;       ///< A handle to a vertex of the mesh.
    typedef typename Mesh::Vertex const *  MeshVertexConstHandle;  ///< A const handle to a vertex of the mesh.
    typedef typename Mesh::Face         *  MeshFaceHandle;         ///< A handle to a face of the mesh.
    typedef typename Mesh::Face const   *  MeshFaceConstHandle;    ///< A const handle to a face of the mesh.

    /** Default constructor. */
    MeshVertexTriple() {}

    /** Constructs the triple from three vertices of a mesh face. */
    MeshVertexTriple(MeshVertexHandle v0, MeshVertexHandle v1, MeshVertexHandle v2, MeshFaceHandle face_, Mesh * mesh_)
    : mesh(mesh_), vertices{ v0, v1, v2 }, face(face_)
    {
      debugAssertM(v0 && v1 && v2, "Mesh triangle: Null vertex provided");
    }

    /** Get the position of any one of the three vertices. */
    Vector3 const & getVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertices[i]->getPosition();
    }

    /** Get the normal at one of the three vertices. */
    Vector3 const & getVertexNormal(int i) const
    {
      return vertices[i]->getNormal();
    }

    /** Get a pointer to any one of the three mesh vertices. */
    MeshVertexConstHandle getMeshVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertices[i];
    }

    /** Get a pointer to any one of the three mesh vertices. */
    MeshVertexHandle getMeshVertex(int i)
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertices[i];
    }

    /** Get the associated mesh face from which the vertices were obtained. */
    MeshFaceConstHandle getMeshFace() const { return face; }

    /** Get the associated mesh face from which the vertices were obtained. */
    MeshFaceHandle getMeshFace() { return face; }

    /** Get the parent mesh. */
    Mesh const * getMesh() const { return mesh; }

    /** Get the parent mesh. */
    Mesh * getMesh() { return mesh; }

  private:
    Mesh * mesh;                   ///< The mesh containing the triangle.
    MeshVertexHandle vertices[3];  ///< The vertices of the triangle.
    MeshFaceHandle face;           ///< The mesh face containing the triangle.

}; // class MeshVertexTriple

/**
 * A set of three vertices of a single face of a display mesh. If the mesh vertices change, so does this triplet.
 *
 * @see DisplayMesh
 */
template <typename MeshT>
class MeshVertexTriple<MeshT, typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type>
{
  public:
    typedef MeshT Mesh;  ///< The mesh type.

    /** Default constructor. */
    MeshVertexTriple() {}

    /** Constructs the triple from three mesh vertices. */
    template <typename IntegerT>
    MeshVertexTriple(IntegerT vi0, IntegerT vi1, IntegerT vi2, Mesh * mesh_, intx tri_index_)
    : mesh(mesh_), mesh_vertices(mesh->getVertices().data()), vertex_indices{ (intx)vi0, (intx)vi1, (intx)vi2 },
      tri_index(tri_index_)
    {
      debugAssertM(vi0 >= 0 && vi0 < mesh->numVertices()
                && vi1 >= 0 && vi1 < mesh->numVertices()
                && vi2 >= 0 && vi2 < mesh->numVertices(), "Display mesh triangle: Mesh vertex index out of bounds");
    }

    /** Get the position of any one of the three vertices. */
    Vector3 const & getVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Display mesh triangle: Vertex index must be 0, 1 or 2");

      return mesh_vertices[(size_t)vertex_indices[i]];
    }

    /**
     * Get the normal at one of the three vertices. If the display mesh does not have explicit vertex normals, the face normal
     * is returned.
     */
    Vector3 getVertexNormal(int i) const
    {
      if (mesh->hasNormals())
        return mesh->getIndexedVertex(vertex_indices[i]).getNormal();
      else
      {
        return (mesh_vertices[(size_t)vertex_indices[1]] - mesh_vertices[(size_t)vertex_indices[0]])
               .cross(mesh_vertices[(size_t)vertex_indices[2]] - mesh_vertices[(size_t)vertex_indices[0]]).normalized();
      }
    }

    /** Get index in the mesh of any one of the three vertices. */
    intx getMeshVertexIndex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Display mesh triangle: Vertex index out of bounds");
      return vertex_indices[i];
    }

    /** Get the index, in the source mesh, of the mesh triangle. */
    intx getMeshTriangleIndex() const { return tri_index; }

    /** Get the parent mesh. */
    Mesh * getMesh() const { return mesh; }

  private:
    Mesh * mesh;                    ///< The mesh containing the triangle.
    Vector3 const * mesh_vertices;  ///< A direct reference to the vertex memory block of the mesh.
    intx vertex_indices[3];         ///< The indices of the vertices of the mesh triangle.
    intx tri_index;                 ///< The index of the triangle in the mesh.

}; // class MeshVertexTriple<DisplayMesh>

namespace MeshTrianglesInternal {

// Add a face of a general mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value >::type
addFace(MeshT & mesh, typename MeshT::Face & face, Array<TriangleT> & tris)
{
  if (face.isTriangle())
  {
    typename MeshT::Face::VertexIterator vi = face.verticesBegin();
    typename MeshT::Vertex * v0 = *(vi++);
    typename MeshT::Vertex * v1 = *(vi++);
    typename MeshT::Vertex * v2 = *vi;
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
  }
  else if (face.isQuad())
  {
    typename MeshT::Vertex * v[4];
    typename MeshT::Face::VertexIterator vi = face.verticesBegin();
    v[0] = *(vi++);
    v[1] = *(vi++);
    v[2] = *(vi++);
    v[3] = *vi;

    intx i0, j0, k0;
    intx i1, j1, k1;
    int num_tris = Polygon3::triangulateQuad(v[0]->getPosition(), v[1]->getPosition(), v[2]->getPosition(), v[3]->getPosition(),
                                             i0, j0, k0, i1, j1, k1);

    if (num_tris > 0)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(v[i0], v[j0], v[k0], &face, &mesh)));

      if (num_tris > 1)
        tris.push_back(TriangleT(typename TriangleT::VertexTriple(v[i1], v[j1], v[k1], &face, &mesh)));
    }
  }
  else
  {
    Array<typename MeshT::Vertex *> face_vertices;
    Polygon3 poly;
    Array<intx> tri_indices;

    intx i = 0;
    for (typename MeshT::Face::VertexIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi, ++i)
    {
      face_vertices.push_back(*vi);
      poly.addVertex((*vi)->getPosition(), i);
    }

    poly.triangulate(tri_indices);
    for (size_t j = 0; j < tri_indices.size(); j += 3)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(face_vertices[(size_t)tri_indices[j]],
                                                                face_vertices[(size_t)tri_indices[j + 1]],
                                                                face_vertices[(size_t)tri_indices[j + 2]],
                                                                &face, &mesh)));
    }
  }
}

// Convert the faces of a general mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value >::type
buildTriangleList(MeshT & mesh, Array<TriangleT> & tris)
{
  for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    addFace<MeshT>(mesh, *fi, tris);
}

// Add a face of a DCEL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value >::type
addFace(MeshT & mesh, typename MeshT::Face & face, Array<TriangleT> & tris)
{
  if (face.isTriangle())
  {
    typename MeshT::Halfedge * he = face.getHalfedge();
    typename MeshT::Vertex * v0 = he->getOrigin();
    typename MeshT::Vertex * v1 = he->next()->getOrigin();
    typename MeshT::Vertex * v2 = he->next()->next()->getOrigin();
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
  }
  else if (face.isQuad())
  {
    typename MeshT::Vertex * v[4];
    typename MeshT::Halfedge * he = face.getHalfedge();
    v[0] = he->getOrigin();
    v[1] = he->next()->getOrigin();
    v[2] = he->next()->next()->getOrigin();
    v[3] = he->next()->next()->next()->getOrigin();

    intx i0, j0, k0;
    intx i1, j1, k1;
    int num_tris = Polygon3::triangulateQuad(v[0]->getPosition(), v[1]->getPosition(), v[2]->getPosition(), v[3]->getPosition(),
                                             i0, j0, k0, i1, j1, k1);

    if (num_tris > 0)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(v[i0], v[j0], v[k0], &face, &mesh)));

      if (num_tris > 1)
        tris.push_back(TriangleT(typename TriangleT::VertexTriple(v[i1], v[j1], v[k1], &face, &mesh)));
    }
  }
  else
  {
    Array<typename MeshT::Vertex *> face_vertices;
    Polygon3 poly;
    Array<intx> tri_indices;

    typename MeshT::Halfedge * he = face.getHalfedge();
    intx num_verts = face.numVertices();
    for (intx i = 0; i < num_verts; ++i)
    {
      face_vertices.push_back(he->getOrigin());
      poly.addVertex(he->getOrigin()->getPosition(), i);
      he = he->next();
    }

    poly.triangulate(tri_indices);
    for (size_t j = 0; j < tri_indices.size(); j += 3)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(face_vertices[(size_t)tri_indices[j]],
                                                                face_vertices[(size_t)tri_indices[j + 1]],
                                                                face_vertices[(size_t)tri_indices[j + 2]],
                                                                &face, &mesh)));
    }
  }
}

// Convert the faces of a DCEL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value >::type
buildTriangleList(MeshT & mesh, Array<TriangleT> & tris)
{
  for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    addFace<MeshT>(mesh, *fi, tris);
}

// Add a face of a display mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type
addFace(MeshT & mesh, typename MeshT::Face & face, Array<TriangleT> & tris)
{
  typedef typename TriangleT::VertexTriple VertexTriple;

  typename MeshT::IndexArray const & tri_indices = mesh.getTriangleIndices();
  size_t beg = 3 * (size_t)face.getFirstTriangle();
  size_t end = beg + 3 * (size_t)face.numTriangles();
  for (size_t i = beg; i < end; i += 3)
  {
    tris.push_back(TriangleT(VertexTriple(tri_indices[i], tri_indices[i + 1], tri_indices[i + 2], &mesh, (intx)i / 3)));
  }
}

// Convert the faces of a display mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type
buildTriangleList(MeshT & mesh, Array<TriangleT> & tris)
{
  typedef MeshT Mesh;
  typedef TriangleT Triangle;
  typedef typename Triangle::VertexTriple VertexTriple;

  typename Mesh::IndexArray const & tri_indices = mesh.getTriangleIndices();
  for (size_t i = 0; i < tri_indices.size(); i += 3)
    tris.push_back(Triangle(VertexTriple(tri_indices[i], tri_indices[i + 1], tri_indices[i + 2], &mesh, (intx)i / 3)));
}

} // namespace MeshTrianglesInternal

/**
 * A set of triangles obtained by triangulating mesh faces. Implemented for general, DCEL and display meshes.
 *
 * @see GeneralMesh, DcelMesh, DisplayMesh
 */
template <typename MeshT>
class MeshTriangles
{
  public:
    THEA_DECL_SMART_POINTERS(MeshTriangles)

    typedef MeshT Mesh;                           ///< The mesh type.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;  ///< A group of meshes.
    typedef MeshVertexTriple<Mesh> VertexTriple;  ///< A triple of mesh vertices.
    typedef Triangle3< VertexTriple > Triangle;   ///< The triangle defined by a triple of mesh vertices.
    typedef Array<Triangle> TriangleArray;    ///< An array of mesh triangles.

    /** Triangulate the faces of a mesh and add them to the set. */
    void add(Mesh & mesh)
    {
      MeshTrianglesInternal::buildTriangleList(mesh, tris);
    }

    /** Triangulate the faces of a mesh group and add them to the set. */
    void add(MeshGroup & mg)
    {
      buildTriangleList(mg, tris);
    }

    /** Add a mesh face to the kd-tree. The tree is <b>not</b> updated until you call init(). */
    void addFace(Mesh & mesh, typename Mesh::Face & face)
    {
      MeshTrianglesInternal::addFace(mesh, face, tris);
    }

    /** Add a single triangle to the set. */
    void addTriangle(Triangle const & tri)
    {
      tris.push_back(tri);
    }

    /** Add a range of triangles to the set. */
    template <typename TriangleIterator> void addTriangles(TriangleIterator tris_begin, TriangleIterator tris_end)
    {
      for (TriangleIterator ti = tris_begin; ti != tris_end; ++ti)
        tris.push_back(*ti);
    }

    /** Check if the set is empty. */
    bool isEmpty() const { return tris.empty(); }

    /** Get the number of triangles in the set. */
    intx numTriangles() const { return (intx)tris.size(); }

    /** Get the triangles in the set. */
    TriangleArray const & getTriangles() const { return tris; }

    /** Get the triangles in the set. */
    TriangleArray & getTriangles() { return tris; }

    /** Clear the set of triangles. */
    void clear()
    {
      tris.clear();
    }

  private:
    /** Convert the faces of a group of meshes to a set of triangles. */
    static void buildTriangleList(MeshGroup & mg, TriangleArray & tris)
    {
      for (typename MeshGroup::MeshIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
        MeshTrianglesInternal::buildTriangleList(**mi, tris);

      for (typename MeshGroup::GroupIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
        buildTriangleList(**ci, tris);
    }

    TriangleArray tris;  ///< Set of mesh triangles.

}; // class MeshTriangles

} // namespace Algorithms
} // namespace Thea

#endif
