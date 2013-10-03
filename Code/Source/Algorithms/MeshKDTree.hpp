//============================================================================
//
// This file is part of the Thea project.
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

#ifndef __Thea_Algorithms_MeshKDTree_hpp__
#define __Thea_Algorithms_MeshKDTree_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"
#include "../Polygon3.hpp"
#include "../Triangle3.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include "KDTree3.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

/**
 * A set of three vertices of a single face of a mesh. Works for general and DCEL meshes.
 *
 * @see GeneralMesh, DCELMesh
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
    {
      debugAssertM(v0 && v1 && v2, "Mesh triangle: Null vertex provided");
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
      face = face_;
      mesh = mesh_;
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
    MeshVertexHandle vertices[3];  ///< The vertices of the triangle.
    MeshFaceHandle face;           ///< The mesh face containing the triangle.
    Mesh * mesh;                   ///< The mesh containing the triangle.

}; // class MeshVertexTriple

namespace MeshKDTreeInternal {

/** Convert a CGAL point to a Vector3. */
template <typename CGALPointT>
Vector3
cgalToVector3(CGALPointT const & p)
{
  return Vector3((Real)p.x(), (Real)p.y(), (Real)p.z());
}

} // namespace MeshKDTreeInternal

/**
 * A set of three vertices of a single face of a CGAL mesh.
 *
 * @see CGALMesh
 */
template <typename MeshT>
class MeshVertexTriple<MeshT, typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type>
{
  public:
    typedef MeshT                               Mesh;                   ///< The mesh type.
    typedef typename Mesh::Vertex_handle        MeshVertexHandle;       ///< A handle to a vertex of the mesh.
    typedef typename Mesh::Vertex_const_handle  MeshVertexConstHandle;  ///< A const handle to a vertex of the mesh.
    typedef typename Mesh::Facet_handle         MeshFaceHandle;         ///< A handle to a face of the mesh.
    typedef typename Mesh::Facet_const_handle   MeshFaceConstHandle;    ///< A const handle to a face of the mesh.

    /** Default constructor. */
    MeshVertexTriple() {}

    /** Constructs the triple from three vertices of a mesh face. */
    MeshVertexTriple(MeshVertexHandle v0, MeshVertexHandle v1, MeshVertexHandle v2, MeshFaceHandle face_, Mesh * mesh_)
    {
      debugAssertM(v0 != MeshVertexHandle()
                && v1 != MeshVertexHandle()
                && v2 != MeshVertexHandle(), "Mesh triangle: Null vertex provided");

      vertices[0] = MeshKDTreeInternal::cgalToVector3(v0->point());
      vertices[1] = MeshKDTreeInternal::cgalToVector3(v1->point());
      vertices[2] = MeshKDTreeInternal::cgalToVector3(v2->point());

      vertex_handles[0] = v0;
      vertex_handles[1] = v1;
      vertex_handles[2] = v2;

      face = face_;
      mesh = mesh_;
    }

    /** Get the position of any one of the three vertices. */
    Vector3 const & getVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertices[i];
    }

    /** Get the normal at one of the three vertices. For CGAL meshes this is just the uniform face normal. */
    Vector3 getVertexNormal(int i) const
    {
      if (Mesh::Facet::Supports_facet_plane)
        return MeshKDTreeInternal::cgalToVector3(face->plane().orthogonal_direction());
      else
        return (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).unit();
    }

    /** Get a handle to any one of the three mesh vertices. */
    MeshVertexConstHandle getMeshVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertex_handles[i];
    }

    /** Get a handle to any one of the three mesh vertices. */
    MeshVertexHandle getMeshVertex(int i)
    {
      debugAssertM(i >= 0 && i < 3, "Mesh triangle: Vertex index out of bounds");
      return vertex_handles[i];
    }

    /** Get the associated mesh face from which the vertices were obtained. */
    MeshFaceConstHandle getMeshFace() const { return face; }

    /** Get the associated mesh face from which the vertices were obtained. */
    MeshFaceHandle getMeshFace() { return face; }

    /** Get the parent mesh. */
    Mesh * getMesh() const { return mesh; }

  private:
    Vector3 vertices[3];                 ///< Positions of the vertices of the mesh triangle.
    MeshVertexHandle vertex_handles[3];  ///< Pointers to the vertices of the mesh triangle.
    MeshFaceHandle face;                 ///< The face containing the triangle.
    Mesh * mesh;                         ///< The mesh containing the triangle.

}; // class MeshVertexTriple<CGALMesh>

/**
 * A set of three vertices of a single face of a display mesh.
 *
 * @see DisplayMesh
 */
template <typename MeshT>
class MeshVertexTriple<MeshT, typename boost::enable_if< Graphics::IsDisplayMesh<MeshT> >::type>
{
  public:
    /** Type of face (enum class). */
    struct FaceType
    {
      /** Supported values. */
      enum Value
      {
        TRIANGLE,
        QUAD
      };

      THEA_ENUM_CLASS_BODY(FaceType)
    };

    /** A handle, in the form of an index/face type pair, to a mesh face. */
    class MeshFaceHandle
    {
      public:
        /** Default constructor. */
        MeshFaceHandle() {}

        /** Construct from an index/type pair. */
        MeshFaceHandle(long index_, FaceType type_) : index(index_), type(type_) {}

        /** Get the index of the face in the source mesh. */
        long getIndex() const { return index; }

        /** Get the type of the face (quad/triangle). */
        FaceType getType() const { return type; }

      private:
        long index;
        FaceType type;

    }; // class MeshFaceHandle

    typedef MeshT           Mesh;                   ///< The mesh type.
    typedef long            MeshVertexHandle;       ///< A handle to a vertex of the mesh.
    typedef long            MeshVertexConstHandle;  ///< A const handle to a vertex of the mesh.
    typedef MeshFaceHandle  MeshFaceConstHandle;    ///< A const handle to a face of the mesh.

    /** Default constructor. */
    MeshVertexTriple() {}

    /** Constructs the triple from three mesh vertices. */
    template <typename IntegerT>
    MeshVertexTriple(IntegerT vi0, IntegerT vi1, IntegerT vi2, Mesh * mesh_, long face_index_, FaceType face_type_)
    : mesh(mesh_), face_index(face_index_), face_type(face_type_)
    {
      typename Mesh::VertexArray const & mv = mesh->getVertices();
      vertices[0] = mv[(array_size_t)vi0];
      vertices[1] = mv[(array_size_t)vi1];
      vertices[2] = mv[(array_size_t)vi2];

      vertex_indices[0] = (long)vi0;
      vertex_indices[1] = (long)vi1;
      vertex_indices[2] = (long)vi2;
    }

    /** Get the position of any one of the three vertices. */
    Vector3 const & getVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Display mesh triangle: Vertex index must be 0, 1 or 2");
      return vertices[i];
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
        return (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).unit();
    }

    /** Get the index of any one of the three mesh vertices. */
    MeshVertexHandle getMeshVertex(int i) const
    {
      debugAssertM(i >= 0 && i < 3, "Display mesh triangle: Vertex index out of bounds");
      return vertex_indices[i];
    }

    /** Get the index, in the source mesh, of the associated mesh face from which the vertices were obtained. */
    long getMeshFaceIndex() const { return face_index; }

    /** Check if the associated mesh face is a triangle or a quad. */
    FaceType getMeshFaceType() const { return face_type; }

    /**
     * Get a handle, in the form of an index/face type pair, of the associated mesh face from which the vertices were obtained.
     */
    MeshFaceHandle getMeshFace() const { return MeshFaceHandle(face_index, face_type); }

    /** Get the parent mesh. */
    Mesh * getMesh() const { return mesh; }

  private:
    Vector3 vertices[3];     ///< The positions of the vertices of the mesh triangle.
    Mesh * mesh;             ///< The mesh containing the triangle.
    long vertex_indices[3];  ///< The indices of the vertices of the mesh triangle.
    long face_index;         ///< The index of the face containing the triangle.
    FaceType face_type;      ///< Type of face (triangle/quad).

}; // class MeshVertexTriple<DisplayMesh>

namespace MeshKDTreeInternal {

// Add a face of a general mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type
addFace(MeshT & mesh, typename MeshT::Face & face, TheaArray<TriangleT> & tris)
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
    typename MeshT::Face::VertexIterator vi = face.verticesBegin();
    typename MeshT::Vertex * v0 = *(vi++);
    typename MeshT::Vertex * v1 = *(vi++);
    typename MeshT::Vertex * v2 = *(vi++);
    typename MeshT::Vertex * v3 = *vi;
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v2, v3, &face, &mesh)));
  }
  else
  {
    TheaArray<typename MeshT::Vertex *> face_vertices;
    Polygon3 poly;
    TheaArray<long> tri_indices;

    long i = 0;
    for (typename MeshT::Face::VertexIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi, ++i)
    {
      face_vertices.push_back(*vi);
      poly.addVertex((*vi)->getPosition(), i);
    }

    poly.triangulate(tri_indices);
    for (array_size_t j = 0; j < tri_indices.size(); j += 3)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(face_vertices[(array_size_t)tri_indices[j]],
                                                                face_vertices[(array_size_t)tri_indices[j + 1]],
                                                                face_vertices[(array_size_t)tri_indices[j + 2]],
                                                                &face, &mesh)));
    }
  }
}

// Convert the faces of a general mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type
buildTriangleList(MeshT & mesh, TheaArray<TriangleT> & tris)
{
  for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    addFace<MeshT>(mesh, *fi, tris);
}

// Add a face of a DCEL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type
addFace(MeshT & mesh, typename MeshT::Face & face, TheaArray<TriangleT> & tris)
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
    typename MeshT::Halfedge * he = face.getHalfedge();
    typename MeshT::Vertex * v0 = he->getOrigin();
    typename MeshT::Vertex * v1 = he->next()->getOrigin();
    typename MeshT::Vertex * v2 = he->next()->next()->getOrigin();
    typename MeshT::Vertex * v3 = he->next()->next()->next()->getOrigin();
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v2, v3, &face, &mesh)));
  }
  else
  {
    TheaArray<typename MeshT::Vertex *> face_vertices;
    Polygon3 poly;
    TheaArray<long> tri_indices;

    typename MeshT::Halfedge * he = face.getHalfedge();
    long num_verts = face.numVertices();
    for (long i = 0; i < num_verts; ++i)
    {
      face_vertices.push_back(he->getOrigin());
      poly.addVertex(he->getOrigin()->getPosition(), i);
      he = he->next();
    }

    poly.triangulate(tri_indices);
    for (array_size_t j = 0; j < tri_indices.size(); j += 3)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(face_vertices[(array_size_t)tri_indices[j]],
                                                                face_vertices[(array_size_t)tri_indices[j + 1]],
                                                                face_vertices[(array_size_t)tri_indices[j + 2]],
                                                                &face, &mesh)));
    }
  }
}

// Convert the faces of a DCEL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type
buildTriangleList(MeshT & mesh, TheaArray<TriangleT> & tris)
{
  for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    addFace<MeshT>(mesh, **fi, tris);
}

// Add a face of a CGAL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type
addFace(MeshT & mesh, typename MeshT::Face & face, TheaArray<TriangleT> & tris)
{
  if (face.is_triangle())
  {
    typename MeshT::Halfedge_handle he = face.halfedge();
    typename MeshT::Vertex_handle v0 = he->vertex();
    typename MeshT::Vertex_handle v1 = he->next()->vertex();
    typename MeshT::Vertex_handle v2 = he->next()->next()->vertex();
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
  }
  else if (face.is_quad())
  {
    typename MeshT::Halfedge_handle he = face.halfedge();
    typename MeshT::Vertex_handle v0 = he->vertex();
    typename MeshT::Vertex_handle v1 = he->next()->vertex();
    typename MeshT::Vertex_handle v2 = he->next()->next()->vertex();
    typename MeshT::Vertex_handle v3 = he->next()->next()->next()->vertex();
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v1, v2, &face, &mesh)));
    tris.push_back(TriangleT(typename TriangleT::VertexTriple(v0, v2, v3, &face, &mesh)));
  }
  else
  {
    TheaArray<typename MeshT::Vertex_handle> face_vertices;
    Polygon3 poly;
    TheaArray<long> tri_indices;

    typename MeshT::Halfedge_handle he = face.halfedge();
    long num_verts = (long)face.facet_degree();
    for (long i = 0; i < num_verts; ++i)
    {
      face_vertices.push_back(he->vertex());
      poly.addVertex(cgalToVector3(he->vertex()->point()), i);
      he = he->next();
    }

    poly.triangulate(tri_indices);
    for (array_size_t j = 0; j < tri_indices.size(); j += 3)
    {
      tris.push_back(TriangleT(typename TriangleT::VertexTriple(face_vertices[(array_size_t)tri_indices[j]],
                                                                face_vertices[(array_size_t)tri_indices[j + 1]],
                                                                face_vertices[(array_size_t)tri_indices[j + 2]],
                                                                &face, &mesh)));
    }
  }
}

// Convert the faces of a CGAL mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type
buildTriangleList(MeshT & mesh, TheaArray<TriangleT> & tris)
{
  for (typename MeshT::Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
    addFace<MeshT>(mesh, *fi, tris);
}

// Add a face of a display mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsDisplayMesh<MeshT> >::type
addFace(MeshT & mesh, typename MeshT::Face & face, TheaArray<TriangleT> & tris)
{
  typedef typename TriangleT::VertexTriple VertexTriple;

  if (face.hasTriangles())
  {
    typename MeshT::IndexArray const & tri_indices = mesh.getTriangleIndices();
    array_size_t beg = 3 * (array_size_t)face.getFirstTriangle();
    array_size_t end = beg + 3 * (array_size_t)face.numTriangles();
    for (array_size_t i = beg; i < end; i += 3)
    {
      tris.push_back(TriangleT(VertexTriple(tri_indices[i], tri_indices[i + 1], tri_indices[i + 2], &mesh, (long)i / 3,
                                            VertexTriple::FaceType::TRIANGLE)));
    }
  }

  if (face.hasQuads())
  {
    typename MeshT::IndexArray const & quad_indices = mesh.getQuadIndices();
    array_size_t beg = 4 * (array_size_t)face.getFirstQuad();
    array_size_t end = beg + 4 * (array_size_t)face.numQuads();
    for (array_size_t i = beg; i < end; i += 4)
    {
      long quad = (long)i / 4;
      tris.push_back(TriangleT(VertexTriple(quad_indices[i], quad_indices[i + 1], quad_indices[i + 2], &mesh, quad,
                                            VertexTriple::FaceType::QUAD)));
      tris.push_back(TriangleT(VertexTriple(quad_indices[i], quad_indices[i + 2], quad_indices[i + 3], &mesh, quad,
                                            VertexTriple::FaceType::QUAD)));
    }
  }
}

// Convert the faces of a display mesh to a set of triangles.
template <typename MeshT, typename TriangleT>
typename boost::enable_if< Graphics::IsDisplayMesh<MeshT> >::type
buildTriangleList(MeshT & mesh, TheaArray<TriangleT> & tris)
{
  typedef MeshT Mesh;
  typedef TriangleT Triangle;
  typedef typename Triangle::VertexTriple VertexTriple;

  typename Mesh::IndexArray const & tri_indices = mesh.getTriangleIndices();
  for (array_size_t i = 0; i < tri_indices.size(); i += 3)
    tris.push_back(Triangle(VertexTriple(tri_indices[i], tri_indices[i + 1], tri_indices[i + 2], &mesh, (long)i / 3,
                            VertexTriple::FaceType::TRIANGLE)));

  // FIXME: Handle non-convex faces
  typename Mesh::IndexArray const & quad_indices = mesh.getQuadIndices();
  for (array_size_t i = 0; i < quad_indices.size(); i += 4)
  {
    long quad = (long)i / 4;
    tris.push_back(Triangle(VertexTriple(quad_indices[i], quad_indices[i + 1], quad_indices[i + 2], &mesh, quad,
                            VertexTriple::FaceType::QUAD)));
    tris.push_back(Triangle(VertexTriple(quad_indices[i], quad_indices[i + 2], quad_indices[i + 3], &mesh, quad,
                            VertexTriple::FaceType::QUAD)));
  }
}

} // namespace MeshKDTreeInternal

/**
 * A kd-tree on mesh triangles. Implemented for general, DCEL, CGAL and display meshes.
 *
 * @see GeneralMesh, DCELMesh, CGALMesh, DisplayMesh
 */
template <typename MeshT, typename NodeAttributeT = NullAttribute>
class MeshKDTree : public Algorithms::KDTree3< Triangle3< MeshVertexTriple<MeshT> >, NodeAttributeT >
{
  private:
    typedef Algorithms::KDTree3< Triangle3< MeshVertexTriple<MeshT> >, NodeAttributeT > BaseT;

  public:
    THEA_DEF_POINTER_TYPES(MeshKDTree, shared_ptr, weak_ptr)

    typedef MeshT Mesh;                           ///< The mesh type.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;  ///< A group of meshes.
    typedef MeshVertexTriple<Mesh> VertexTriple;  ///< A triple of mesh vertices.
    typedef Triangle3< VertexTriple > Triangle;   ///< The triangle defined by a triple of mesh vertices.
    typedef TheaArray<Triangle> TriangleArray;    ///< An array of mesh triangles.

    /**
     * Add a mesh to the kd-tree. The mesh is converted to triangles which are cached internally. The tree is <b>not</b>
     * actually constructed until you call init().
     */
    void add(Mesh & mesh)
    {
      MeshKDTreeInternal::buildTriangleList(mesh, tris);
    }

    /**
     * Add a group of meshes to the kd-tree. The meshes are converted to triangles which are cached internally. The tree is
     * <b>not</b> actually constructed until you call init().
     */
    void add(MeshGroup & mg)
    {
      buildTriangleList(mg, tris);
    }

    /** Add a mesh face to the kd-tree. The tree is <b>not</b> updated until you call init(). */
    void addFace(Mesh & mesh, typename Mesh::Face & face)
    {
      MeshKDTreeInternal::addFace(mesh, face, tris);
    }

    /** Add a single triangle to the kd-tree. The tree is <b>not</b> updated until you call init(). */
    void addTriangle(Triangle const & tri)
    {
      tris.push_back(tri);
    }

    /**
     * Add a set of triangles to the kd-tree. TriangleIterator must dereference to Triangle. The tree is <b>not</b> updated
     * until you call init().
     */
    template <typename TriangleIterator> void addTriangles(TriangleIterator tris_begin, TriangleIterator tris_end)
    {
      for (TriangleIterator ti = tris_begin; ti != tris_end; ++ti)
        tris.push_back(*ti);
    }

    /**
     * Get direct access to the triangles cached by the tree, which will be used to build the tree on the next call to init().
     */
    TriangleArray const & getTriangles() const { return tris; }

    /**
     * Get direct access to the triangles cached by the tree, which will be used to build the tree on the next call to init().
     */
    TriangleArray & getTriangles() { return tris; }

    /**
     * Compute the kd-tree from the added meshes. You <b>must</b> call this function to construct (or recompute) the tree after
     * any addMesh() or addMeshGroup() calls. Clears the triangle cache.
     */
    void init(int max_depth = -1, int max_elems_in_leaf = -1, bool save_memory = false, bool deallocate_previous_memory = true)
    {
      BaseT::init(tris.begin(), tris.end(), max_depth, max_elems_in_leaf, save_memory, deallocate_previous_memory);
      tris.clear();
    }

    /**
     * Clear the tree. If \a deallocate_all_memory is false, memory allocated in pools is held to be reused if possible by the
     * next init() operation.
     */
    void clear(bool deallocate_all_memory = true)
    {
      BaseT::clear(deallocate_all_memory);
      tris.clear();
    }

  private:
    /** Convert the faces of a group of meshes to a set of triangles. */
    static void buildTriangleList(MeshGroup & mg, TriangleArray & tris)
    {
      for (typename MeshGroup::MeshIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
        MeshKDTreeInternal::buildTriangleList(**mi, tris);

      for (typename MeshGroup::GroupIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
        buildTriangleList(**ci, tris);
    }

    TriangleArray tris;  ///< Internal cache of triangles used to initialize the tree.

}; // class MeshKDTree

} // namespace Algorithms
} // namespace Thea

#endif
