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

#ifndef __Thea_Graphics_DisplayMesh_hpp__
#define __Thea_Graphics_DisplayMesh_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Colors.hpp"
#include "../MatrixWrapper.hpp"
#include "../UnorderedMap.hpp"
#include "IMesh.hpp"
#include "DefaultMeshCodecs.hpp"
#include "IncrementalDisplayMeshBuilder.hpp"
#include <array>

namespace Thea {
namespace Graphics {

// Forward declaration
class DisplayMesh;

/**
 * A vertex of of a DisplayMesh. This is created on the fly when a vertex is accessed, vertex data is not actually stored in
 * this format.
 */
class THEA_API DisplayMeshVertex
{
  public:
    typedef DisplayMesh Mesh;  ///< Parent mesh class.

  private:
    DisplayMesh  *  mesh;
    Vector3      *  point;
    Vector3      *  normal;
    ColorRgba    *  color;
    Vector2      *  texcoord;

  public:
    /** Default constructor. Creates an invalid vertex reference. Only for compatibility with standard containers. */
    DisplayMeshVertex() : mesh(nullptr), point(nullptr) {}

    /** Constructor. */
    DisplayMeshVertex(DisplayMesh * mesh_, Vector3 & point_, Vector3 * normal_ = nullptr, ColorRgba * color_ = nullptr,
                      Vector2 * texcoord_ = nullptr)
    : mesh(mesh_), point(&point_), normal(normal_), color(color_), texcoord(texcoord_)
    {}

    /** Get the parent mesh. */
    DisplayMesh const * getMesh() const { return mesh; }

    /** Get the parent mesh. */
    DisplayMesh * getMesh() { return mesh; }

    /** Check if the vertex is valid. */
    operator bool() const { return mesh && point; }

    /** Get the position of the vertex. */
    Vector3 const & getPosition() const { return *point; }

    /** Set the position of the vertex. */
    void setPosition(Vector3 const & point_);

    /** Check if the vertex has a normal. */
    bool hasNormal() const { return normal != nullptr; }

    /** Get the normal at the vertex. Call only if hasNormal() returns true. */
    Vector3 const & getNormal() const
    {
      debugAssertM(hasNormal(), "DisplayMeshVertex: Vertex does not have a normal");
      return *normal;
    }

    /** Set the normal at the vertex. Call only if hasNormal() returns true. */
    void setNormal(Vector3 const & normal_);

    /** Check if the vertex has a color. */
    bool hasColor() const { return color != nullptr; }

    /** Get the color at the vertex. Call only if hasColor() returns true. */
    ColorRgba const & getColor() const
    {
      debugAssertM(hasColor(), "DisplayMeshVertex: Vertex does not have a color");
      return *color;
    }

    /** Set the color at the vertex. Call only if hasColor() returns true. */
    void setColor(ColorRgba const & color_);

    /** Check if the vertex has a texture coordinate. */
    bool hasTexCoord() const { return texcoord != nullptr; }

    /**
     * Get the texture coordinates at the vertex, or null if no such texture coordinates exist. Call only if hasTexCoord()
     * returns true.
     */
    Vector2 const & getTexCoord() const
    {
      debugAssertM(hasTexCoord(), "DisplayMeshVertex: Vertex does not have texture coordinates");
      return *texcoord;
    }

    /** Set the texture coordinates at the vertex. Call only if hasTexCoord() returns true. */
    void setTexCoord(Vector2 const & texcoord_);

}; // class DisplayMeshVertex

/**
 * A reference to a vertex of a display mesh, via its index. This is a convenience class that is better in some cases than
 * DisplayMeshVertex (e.g. when references are created on the fly as vertices are added to the mesh (adding vertices invalidates
 * previously created DisplayMeshVertex objects), or when memory conservation is important.).
 */
class THEA_API DisplayMeshIndexedVertex
{
  public:
    typedef DisplayMesh Mesh;  ///< Parent mesh class.

    /** Constructor. */
    DisplayMeshIndexedVertex(DisplayMesh * mesh_ = nullptr, intx index_ = -1) : mesh(mesh_), index(index_) {}

    /** Check if the vertex reference is valid. */
    operator bool() const { return mesh && index >= 0; }

    /** Get the parent mesh. */
    DisplayMesh const * getMesh() const { return mesh; }

    /** Get the parent mesh. */
    DisplayMesh * getMesh() { return mesh; }

    /** Get the index of the vertex in the parent mesh. */
    intx getIndex() const { return index; }

    /** Get a direct (pointer-based) reference to the vertex in the parent mesh. */
    DisplayMeshVertex getVertex();

    /** Get the position of the vertex. */
    Vector3 const & getPosition() const;

    /** Set the position of the vertex. */
    void setPosition(Vector3 const & point_);

    /** Check if the vertex has a normal. */
    bool hasNormal() const;

    /** Get the normal at the vertex. Call only if hasNormal() returns true. */
    Vector3 const & getNormal() const;

    /** Set the normal at the vertex. Call only if hasNormal() returns true. */
    void setNormal(Vector3 const & normal_);

    /** Check if the vertex has a color. */
    bool hasColor() const;

    /** Get the color at the vertex. Call only if hasColor() returns true. */
    ColorRgba const & getColor() const;

    /** Set the color at the vertex. Call only if hasColor() returns true. */
    void setColor(ColorRgba const & color_);

    /** Check if the vertex has a texture coordinate. */
    bool hasTexCoord() const;

    /**
     * Get the texture coordinates at the vertex, or null if no such texture coordinates exist. Call only if hasTexCoord()
     * returns true.
     */
    Vector2 const & getTexCoord() const;

    /** Set the texture coordinates at the vertex. Call only if hasTexCoord() returns true. */
    void setTexCoord(Vector2 const & texcoord_);

  private:
    DisplayMesh * mesh;
    intx index;

}; // class DisplayMeshIndexedVertex

/**
 * A face of a DisplayMesh. This is created on the fly when a face is accessed, face data is not actually stored in this format.
 * It can be freely copied and used as a handle to a face.
 */
class THEA_API DisplayMeshFace
{
  public:
    typedef DisplayMesh Mesh;  ///< Parent mesh class.

    /** Default constructor. Creates an invalid face. */
    DisplayMeshFace() : mesh(nullptr), num_vertices(0), first_tri(-1), num_triangles(0) {}

    /** Constructor. */
    DisplayMeshFace(DisplayMesh * mesh_, int num_vertices_, intx first_tri_, int num_triangles_)
    : mesh(mesh_), num_vertices(num_vertices_), first_tri(first_tri_), num_triangles(num_triangles_)
    {
      debugAssertM(mesh_, "DisplayMeshFace: Parent mesh must exist");
      debugAssertM(first_tri_ >= 0, "DisplayMeshFace: Index of first triangle must be non-negative");
    }

    /** Get the parent mesh. */
    DisplayMesh const * getMesh() const { return mesh; }

    /** Get the parent mesh. */
    DisplayMesh * getMesh() { return mesh; }

    /** Check if the face is valid. */
    operator bool() const { return mesh && first_tri >= 0; }

    /** Get the number of vertices of the face. */
    int numVertices() const { return num_vertices; }

    /** Check if the face is a single triangle. */
    bool isTriangle() const { return num_vertices == 3; }

    /** Check if the face is a quad. */
    bool isQuad() const { return num_vertices == 4; }

    /**
     * Get the <tt>i</tt>'th vertex of the face, in the original winding order but not necessarily starting at the same first
     * vertex.
     */
    uint32 getVertexIndex(int i) const;

    /** Get the number of triangles in the face. */
    int numTriangles() const { return num_triangles; }

    /** Get the index of the first triangle of the face. */
    intx getFirstTriangle() const { return first_tri; }

  private:
    DisplayMesh * mesh;
    int num_vertices;
    intx first_tri;
    int num_triangles;

}; // class DisplayMeshFace

/** A class for storing meshes for display, without detailed topology information. */
class THEA_API DisplayMesh : public NamedObject, public virtual IMesh
{
  public:
    THEA_DECL_SMART_POINTERS(DisplayMesh)

    /** Mesh type tag. */
    struct DISPLAY_MESH_TAG {};

    typedef Array<Vector3>    VertexArray;    ///< Array of vertex positions.
    typedef Array<Vector3>    NormalArray;    ///< Array of normals.
    typedef Array<Vector2>    TexCoordArray;  ///< Array of texture coordinates.
    typedef Array<ColorRgba>  ColorArray;     ///< Array of colors.
    typedef Array<uint32>     IndexArray;     ///< Array of indices.

    typedef DisplayMeshVertex Vertex;  ///< A convenience wrapper for accessing a vertex's properties.
    typedef DisplayMeshIndexedVertex IndexedVertex;  ///< A reference to a vertex's properties via its index.
    typedef DisplayMeshFace Face;  ///< A convenience wrapper for accessing a face's properties.

    typedef std::array<intx, 3> IndexTriple;  ///< Vertex indices of a single triangle.

    // Generic typedefs, each mesh class must define these for builder and codec compatibility
    typedef intx VertexHandle;       ///< Handle to a mesh vertex.
    typedef intx VertexConstHandle;  ///< Handle to an immutable mesh vertex.
    typedef intx FaceHandle;         ///< Handle to a mesh face.
    typedef intx FaceConstHandle;    ///< Handle to an immutable mesh face.

    /** Identifiers for various GPU buffers (enum class). */
    struct BufferId
    {
      /** Supported values. */
      enum Value
      {
        ALL              =  0xFFFF,  // handles upto 31 buffers, should be enough
        VERTEX_POSITION  =  0x0001,
        VERTEX_NORMAL    =  0x0002,
        VERTEX_COLOR     =  0x0004,
        VERTEX_TEXCOORD  =  0x0008,
        TRIANGLE         =  0x0010,
      };

      THEA_ENUM_CLASS_BODY(BufferId)

    }; // struct BufferId

  private:
    // Vertex data
    VertexArray    vertices;   ///< Vertex positions.
    NormalArray    normals;    ///< Vertex normals.
    ColorArray     colors;     ///< Vertex colors.
    TexCoordArray  texcoords;  ///< Vertex texture coordinates.

    // Face data
    IndexArray tris;   ///< Triangle indices (in triplets).

    // Edge data
    IndexArray edges;  ///< Edge indices (in pairs).

    // Element source indices (typically from source files)
    Array<intx> vertex_src_indices;     ///< Index in source file of each vertex.
    Array<intx> tri_src_face_indices;   ///< Index in source file of the parent face for each triangle.

    /**
     * The index of the first triangle in each polygonal face. Has one additional entry to store the total number of triangles
     * so that the triangle count of face <tt>i</tt> is always <tt>face_starting_tris[i + 1] - face_starting_tris[i]</tt>.
     */
    Array<intx> face_starting_tris;

    typedef UnorderedMap< intx, Array<intx> > PolyVertexMap;  ///< Map from first triangle of a polygon to its full vertex list.
    PolyVertexMap large_poly_verts;  ///< Map the index of the first triangle of a polygon to the latter's full vertex list.

    bool valid_bounds;  ///< Is the bounding box valid?
    AxisAlignedBox3 bounds;  ///< Bounding box.

    bool wireframe_enabled;  ///< If true, mesh edges will be packed in buffers and sent to the GPU.

    int changed_buffers;        ///< A bitwise OR of the flags of the buffers that have changed.
    IBufferPool * buf_pool;     ///< GPU buffer pool.
    IBuffer * vertices_buf;     ///< GPU buffer for vertex positions.
    IBuffer * tris_buf;         ///< GPU buffer for triangle indices.
    IBuffer * normals_buf;      ///< GPU buffer for vertex normals.
    IBuffer * colors_buf;       ///< GPU buffer for vertex colors.
    IBuffer * texcoords_buf;    ///< GPU buffer for texture coordinates.
    IBuffer * edges_buf;        ///< GPU buffer for edges.

    // Map the vertex and index buffers as matrices for the IMesh API
    typedef MatrixMap<3, Eigen::Dynamic, Real,   MatrixLayout::COLUMN_MAJOR>  VertexMatrix;    ///< Wraps vertices as a matrix.
    typedef MatrixMap<3, Eigen::Dynamic, uint32, MatrixLayout::COLUMN_MAJOR>  TriangleMatrix;  ///< Wraps triangles as a matrix.

    mutable VertexMatrix    vertex_matrix;  ///< Vertex data as a dense 3xN column-major matrix.
    mutable TriangleMatrix  tri_matrix;     ///< Triangle indices as a dense 3xN column-major matrix.

    mutable MatrixWrapper<VertexMatrix>    vertex_wrapper;
    mutable MatrixWrapper<TriangleMatrix>  tri_wrapper;

    // Friend classes
    friend class DisplayMeshVertex;
    friend class DisplayMeshIndexedVertex;
    friend class DisplayMeshFace;

  public:
    /** Constructor. */
    DisplayMesh(std::string const & name = "AnonymousMesh");

    /** Copy constructor. */
    DisplayMesh(DisplayMesh const & src);

    /** Destructor. */
    virtual ~DisplayMesh() {}

    // Abstract mesh interface
    IDenseMatrix<Real> const * THEA_ICALL getVertexMatrix() const;
    IDenseMatrix<uint32> const * THEA_ICALL getTriangleMatrix() const;
    IDenseMatrix<uint32> const * THEA_ICALL getQuadMatrix() const;

    /**
     * Get a read-write handle to the 3xN column-major matrix of vertex positions. The matrix is not resizable. Remember to
     * call invalidateBounds() and invalidateGpuBuffers(BufferId::VERTEX_POSITION) after making any changes to the vertices.
     */
    IDenseMatrix<Real> * getVertexMatrix();

    /**
     * Get a read-write handle to the 3 x &#35;triangles column-major matrix of triangle vertex indices. The matrix is not
     * resizable. Remember to call invalidateGpuBuffers(BufferId::TRIANGLE) after making any changes to the indices.
     */
    IDenseMatrix<uint32> * getTriangleMatrix();

    /** Get the set of vertex positions. */
    VertexArray const & getVertices() const { return vertices; }

    /** Get the set of vertex normals. */
    NormalArray const & getNormals() const { return normals; }

    /** Get the set of vertex colors. */
    ColorArray const & getColors() const { return colors; }

    /** Get the set of vertex texture coordinates. */
    TexCoordArray const & getTexCoords() const { return texcoords; }

    /**
     * Get a structure referencing the position, normal, color and texture coordinates of a single vertex. This may be used to
     * access and modify the properties of a vertex.
     */
    Vertex getVertex(intx i);

    /**
     * Get a structure referencing the position, normal, color and texture coordinates of a single vertex, via its index. This
     * may be used to access and modify the properties of a vertex.
     */
    IndexedVertex getIndexedVertex(intx i) { return IndexedVertex(this, i); }

    /** Get the vertex indices of the triangular faces. Each successive triplet defines a triangle. */
    IndexArray const & getTriangleIndices() const { return tris; }

    /**
     * Get a structure referencing the three indices of a single triangle. This is a convenience function -- the same
     * information can be obtained from getTriangleIndices().
     */
    IndexTriple getTriangle(intx tri_index) const;

    /**
     * Get the source index of a given vertex, or negative if none exists. This is typically the index of the vertex in the
     * source mesh file.
     */
    intx getVertexSourceIndex(intx i) const
    {
      return i >= 0 && i < (intx)vertex_src_indices.size() ? vertex_src_indices[(size_t)i] : -1;
    }

    /**
     * Get the index of the source face of a given triangle, or negative if none exists. This is typically the index of the face
     * in the source mesh file.
     */
    intx getTriangleSourceFaceIndex(intx tri_index) const
    {
      return tri_index >= 0 && tri_index < (intx)tri_src_face_indices.size() ? tri_src_face_indices[(size_t)tri_index] : -1;
    }

    /** Deletes all data in the mesh. */
    virtual void clear();

    /** True if and only if the mesh contains no objects. */
    bool isEmpty() const { return vertices.empty() && tris.empty(); }

    /** Get the number of vertices of the mesh. */
    intx numVertices() const { return (intx)vertices.size(); };

    /**
     * Get the number of mesh triangles, obtained <i>after</i> triangulating any higher-degree faces.
     *
     * @see numFaces()
     */
    intx numTriangles() const { return (intx)(tris.size() / 3); };

    /**
     * Get the number of polygonal faces in the mesh.
     *
     * @see numTriangles()
     */
    intx numFaces() const { return (intx)face_starting_tris.size() - 1; };

    /** Recompute and cache the bounding box for the mesh. Make sure this has been called before calling getBounds(). */
    void updateBounds();

    /**
     * Get the cached bounding box of the mesh. Will be out-of-date unless updateBounds() has been called after all
     * modifications.
     */
    AxisAlignedBox3 const & getBounds() const
    {
      const_cast<DisplayMesh *>(this)->updateBounds();
      return bounds;
    }

    /** Check if the vertices have attached normal information. */
    bool hasNormals() const { return !normals.empty(); }

    /** Check if the vertices have attached color information. */
    bool hasColors() const { return !colors.empty(); }

    /** Check if the vertices have attached texture coordinate information. */
    bool hasTexCoords() const { return !texcoords.empty(); }

    /**
     * Adds a color property to each vertex, initially initialized to ColorRgba(0, 0, 0, 0). If a color property already exists,
     * no action is taken.
     */
    virtual void addColors();

    /**
     * Adds a normal property to each vertex, initially initialized to Vector3(0, 0, 0). If a normal property already exists, no
     * action is taken.
     */
    virtual void addNormals();

    /**
     * Adds a texture coordinate property to each vertex, initially initialized to Vector2(0, 0). If a texture coordinate
     * property already exists, no action is taken.
     */
    virtual void addTexCoords();

    /**
     * Add a vertex to the mesh, with optional normal, color and texture coordinate attributes, as well as an optional source
     * index (typically the index of the vertex in the mesh source file).
     *
     * Each of these attributes is (independently) an all or nothing choice for the mesh: it must be specified for all vertices,
     * or no vertices. E.g. if the normal is specified, it must be specified for all existing and future vertices of the mesh.
     * Similarly, if the source index is non-negative, it must be so for all vertices.
     *
     * @return The index of the new vertex in the mesh (distinct from the source index input to this function). Indices are
     *   guaranteed to be sequentially generated, starting from 0.
     */
    virtual intx addVertex(Vector3 const & point, intx src_index = -1, Vector3 const * normal = nullptr,
                           ColorRgba const * color = nullptr, Vector2 const * texcoord = nullptr);

    /**
     * Add a triangular face to the mesh, specified by three vertex indices and an optional source face index (typically the
     * index of the face in the mesh source file).
     *
     * @return The index of the new face in the mesh, i.e. a valid argument to getFace(), or a negative number on error. Note,
     *   this is <b>NOT</b> necessarily the same as the index of the triangle in the triangle list, which can be computed as
     *   numTriangles() <i>before</i> the addition. Nor is it necessarily the same as \a src_face_index.
     */
    virtual intx addTriangle(intx vi0, intx vi1, intx vi2, intx src_face_index = -1);

    /**
     * Add a polygonal face to the mesh, specified as a sequence of vertex indices and an optional source face index (typically
     * the index of the face in the mesh source file). Polygons with less than 3 vertices are ignored. If the polygon has 3
     * vertices, it is added to the triangle list. If it has more than 3 vertices, it is triangulated and added to the triangle
     * list.
     *
     * @return The index of the new face in the mesh, i.e. a valid argument to getFace(), or a negative number on error. This is
     *   not necessarily the same as \a src_face_index.
     */
    virtual intx addFace(int num_vertices, intx const * vertex_indices, intx src_face_index = -1);

    /**
     * Add a polygonal face to the mesh, specified as a sequence of vertex indices obtained by dereferencing [vbegin, vend), and
     * an optional source face index (typically the index of the face in the mesh source file). Polygons with less than 3
     * vertices are ignored. If the polygon has 3 vertices, it is added to the triangle list. If it has more than 3 vertices, it
     * is triangulated and added to the triangle list.
     *
     * @return The index of the new face in the mesh, i.e. a valid argument to getFace(), or a negative number on error. This is
     *   not necessarily the same as \a src_face_index.
     */
    template <typename IndexIterator> intx addFace(IndexIterator vi_begin, IndexIterator vi_end, intx src_face_index = -1)
    {
      face_vertex_indices.clear();
      for (IndexIterator vi = vi_begin; vi != vi_end; ++vi)
        face_vertex_indices.push_back((intx)*vi);

      return addFace((int)face_vertex_indices.size(), &face_vertex_indices[0], src_face_index);
    }

    /**
     * Get a handle to a polygonal face which may comprise several triangles.
     *
     * @param face_index The index of the face in the mesh, as returned by the call to addFace() which added the face.
     */
    Face getFace(intx face_index);

    /** Check if a particular face consists of a single triangle or not. */
    bool isTriangle(intx face_index) const
    {
      debugAssertM(face_index >= 0 && face_index + 1 < (intx)face_starting_tris.size(),
                   getNameStr() + ": Face index out of bounds");

      return face_starting_tris[(size_t)face_index + 1] - face_starting_tris[(size_t)face_index] == 1;
    }

    /**
     * Check if a particular face is a non-degenerate quad or not. This requires that it was created (addFace()) with 4
     * vertices, and that it produced two triangles.
     */
    bool isQuad(intx face_index) const
    {
      debugAssertM(face_index >= 0 && face_index + 1 < (intx)face_starting_tris.size(),
                   getNameStr() + ": Face index out of bounds");

      auto begin = face_starting_tris[(size_t)face_index], end = face_starting_tris[(size_t)face_index + 1];
      return end - begin == 2 && large_poly_verts.find(begin) == large_poly_verts.end();
    }

    /** Set the position of a mesh vertex. */
    virtual void setVertex(intx vertex_index, Vector3 const & position)
    {
      alwaysAssertM(vertex_index >= 0 && vertex_index < (intx)vertices.size(),
                    getNameStr() + ": Vertex index out of bounds");

      vertices[(size_t)vertex_index] = position;
      invalidateBounds();
      invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_POSITION);
    }

    /** Set the normal of a mesh vertex. */
    virtual void setNormal(intx vertex_index, Vector3 const & normal)
    {
      alwaysAssertM(vertex_index >= 0 && vertex_index < (intx)normals.size(),
                    getNameStr() + ": Vertex index out of bounds, or vertex does not have associated normal field");

      normals[(size_t)vertex_index] = normal;
      invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_NORMAL);
    }

    /** Set the color of a mesh vertex. */
    virtual void setColor(intx vertex_index, ColorRgba const & color)
    {
      alwaysAssertM(vertex_index >= 0 && vertex_index < (intx)colors.size(),
                    getNameStr() + ": Vertex index out of bounds, or vertex does not have associated color field");

      colors[(size_t)vertex_index] = color;
      invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_COLOR);
    }

    /** Set the texture coordinates of a mesh vertex. */
    virtual void setTexCoord(intx vertex_index, Vector2 const & texcoord)
    {
      alwaysAssertM(vertex_index >= 0 && vertex_index < (intx)texcoords.size(),
                    getNameStr()
                  + ": Vertex index out of bounds, or vertex does not have associated texture coordinates field");

      texcoords[(size_t)vertex_index] = texcoord;
      invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_TEXCOORD);
    }

    /** Set the normal at each vertex as the average of the normals of all faces incident at the vertex. */
    virtual void computeAveragedVertexNormals();

    /** Flip the direction of all mesh normals. Has no effect if the mesh vertices do not have normal information. */
    virtual void flipNormals();

    /**
     * Duplicate vertices (and normals, colors and texture coordinates) as necessary to ensure no two triangles share
     * vertex-level data. Useful for drawing faces with per-face information. The faces retain their sequence, but reference new
     * sets of vertex data.
     */
    virtual void isolateTriangles();

    /**
     * Enable/disable drawing the edges of the mesh. Enabling this function will <b>not</b> draw any edges unless you turn on
     * the appropriate RenderOptions flag. The edges will be uploaded to the graphics system on the next call to
     * uploadToGraphicsSystem().
     *
     * Wireframe mode is initially disabled to save video memory.
     *
     * @see wireframeIsEnabled()
     */
    void setWireframeEnabled(bool value)
    {
      if (value && !wireframe_enabled)
        wireframe_enabled = true;
      else if (!value && wireframe_enabled)
        wireframe_enabled = false;

      invalidateGpuBuffers();
    }

    /**
     * Check if wireframe drawing is enabled. Wireframe mode is initially disabled.
     *
     * @see setWireframeEnabled()
     */
    bool wireframeIsEnabled() const { return wireframe_enabled; }

    /** Invalidate part or all of the current GPU data for the mesh. */
    void invalidateGpuBuffers(int changed_buffers_ = BufferId::ALL) { changed_buffers |= changed_buffers_; }

    /** Invalidate the current bounding box of the mesh. */
    void invalidateBounds() { valid_bounds = false; }

    int8 THEA_ICALL draw(IRenderSystem * render_system, IRenderOptions const * options = nullptr) const;

  protected:
    /** Check if a buffer has is synchronized with the mesh or not. */
    bool gpuBufferIsValid(BufferId buffer) const { return (changed_buffers & (int)buffer) == 0; }

    /** Clear the set of changed buffers. */
    void allGpuBuffersAreValid() { changed_buffers = 0; }

    /** Upload GPU resources to the graphics system. */
    bool uploadToGraphicsSystem(IRenderSystem & render_system);

  private:
    /** Update the set of edge indices from face data. */
    void updateEdges();

    // Temporary storage for triangulating polygons with more than 4 vertices
    Array<intx> face_vertex_indices;
    Array<intx> triangulated_indices;

}; // class DisplayMesh

// Inline functions that require the full declaration of DisplayMesh.

//============================================================================================================================
// DisplayMeshVertex
//============================================================================================================================

inline void
DisplayMeshVertex::setPosition(Vector3 const & point_)
{
  *point = point_;
  mesh->invalidateBounds();
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_POSITION);
}

inline void
DisplayMeshVertex::setNormal(Vector3 const & normal_)
{
  debugAssertM(hasNormal(), "DisplayMeshVertex: Vertex does not have a normal");
  *normal = normal_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_NORMAL);
}

inline void
DisplayMeshVertex::setColor(ColorRgba const & color_)
{
  debugAssertM(hasColor(), "DisplayMeshVertex: Vertex does not have a color");
  *color = color_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_COLOR);
}

inline void
DisplayMeshVertex::setTexCoord(Vector2 const & texcoord_)
{
  debugAssertM(hasTexCoord(), "DisplayMeshVertex: Vertex does not have texture coordinates");
  *texcoord = texcoord_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_TEXCOORD);
}

//============================================================================================================================
// DisplayMeshIndexedVertex
//============================================================================================================================

inline DisplayMeshVertex
DisplayMeshIndexedVertex::getVertex()
{
  debugAssertM(*this, "DisplayMeshIndexedVertex: Cannot get vertex from invalid reference");
  return mesh->getVertex(index);
}

inline Vector3 const &
DisplayMeshIndexedVertex::getPosition() const
{
  return mesh->vertices[(size_t)index];
}

inline void
DisplayMeshIndexedVertex::setPosition(Vector3 const & point_)
{
  mesh->vertices[(size_t)index] = point_;
  mesh->invalidateBounds();
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_POSITION);
}

inline bool
DisplayMeshIndexedVertex::hasNormal() const
{
  return mesh->hasNormals();
}

inline Vector3 const &
DisplayMeshIndexedVertex::getNormal() const
{
  debugAssertM(hasNormal(), "DisplayMeshIndexedVertex: Vertex does not have a normal");
  return mesh->normals[(size_t)index];
}

inline void
DisplayMeshIndexedVertex::setNormal(Vector3 const & normal_)
{
  debugAssertM(hasNormal(), "DisplayMeshIndexedVertex: Vertex does not have a normal");
  mesh->normals[(size_t)index] = normal_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_NORMAL);
}

inline bool
DisplayMeshIndexedVertex::hasColor() const
{
  return mesh->hasColors();
}

inline ColorRgba const &
DisplayMeshIndexedVertex::getColor() const
{
  debugAssertM(hasColor(), "DisplayMeshIndexedVertex: Vertex does not have a color");
  return mesh->colors[(size_t)index];
}

inline void
DisplayMeshIndexedVertex::setColor(ColorRgba const & color_)
{
  debugAssertM(hasColor(), "DisplayMeshIndexedVertex: Vertex does not have a color");
  mesh->colors[(size_t)index] = color_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_COLOR);
}

inline bool
DisplayMeshIndexedVertex::hasTexCoord() const
{
  return mesh->hasTexCoords();
}

inline Vector2 const &
DisplayMeshIndexedVertex::getTexCoord() const
{
  debugAssertM(hasTexCoord(), "DisplayMeshIndexedVertex: Vertex does not have texture coordinates");
  return mesh->texcoords[(size_t)index];
}

inline void
DisplayMeshIndexedVertex::setTexCoord(Vector2 const & texcoord_)
{
  debugAssertM(hasTexCoord(), "DisplayMeshIndexedVertex: Vertex does not have texture coordinates");
  mesh->texcoords[(size_t)index] = texcoord_;
  mesh->invalidateGpuBuffers(DisplayMesh::BufferId::VERTEX_TEXCOORD);
}

//============================================================================================================================
// DisplayMeshFace
//============================================================================================================================

inline uint32
DisplayMeshFace::getVertexIndex(int i) const
{
  debugAssertM(i >= 0 && i < num_vertices, "DisplayMeshFace: Vertex index out of bounds");

  switch (num_vertices)
  {
    case 3:
    case 4: return mesh->tris[3 * first_tri + (i < 3 ? i : 4)];  // exploit vertex ordering guarantee of triangulateQuad()

    default:
    {
      auto loc = mesh->large_poly_verts.find(first_tri);
      debugAssertM(loc != mesh->large_poly_verts.end(), "DisplayMeshFace: Vertices of high-degree polygon not found");
      return loc->second[(size_t)i];
    }
  }
}

} // namespace Graphics
} // namespace Thea

#endif
