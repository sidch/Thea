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

#ifndef __Thea_Algorithms_ImlsSurface_hpp__
#define __Thea_Algorithms_ImlsSurface_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Ball3.hpp"
#include "../MatVec.hpp"
#include "../Noncopyable.hpp"
#include "../Polygon3.hpp"
#include "../Triangle3.hpp"
#include "../UnorderedMap.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include "KdTreeN.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

namespace ImlsSurfaceInternal {

// Three indexed vertices plus base position array.
class IndexedVertexTriple
{
  private:
    Array<Vector3> const * array;
    size_t indices[3];

  public:
    /** Default constructor. */
    IndexedVertexTriple() {}

    /** Initializing constructor. */
    IndexedVertexTriple(Array<Vector3> const * array_, size_t i0, size_t i1, size_t i2)
    {
      array = array_;
      indices[0] = i0;
      indices[1] = i1;
      indices[2] = i2;
    }

    /** Get the i'th index. */
    size_t getIndex(int i) const { return indices[i]; }

    /** Get the i'th vertex. */
    Vector3 const & getVertex(int i) const { return (*array)[indices[i]]; }

}; // IndexedVertexTriple

typedef Triangle3<IndexedVertexTriple> IndexedTriangle;  ///< Triangle with indexed vertices.

// KD-tree node attributes.
struct NodeAttribute
{
  double    area;        ///< Surface area of the contents of the node
  double    unweighted;  ///< Unweighted integrals of phi over the surface contained by the node.
  Vector3d  centroid;    ///< Area weighted centroid of the node.
  Vector3d  normal;      ///< Area weighted normal of the node.
  double    normal_len;  ///< Cached magnitude of the normal.
};

typedef KdTreeN<IndexedTriangle, 3, Real, NodeAttribute> TriangleKdTree;

// A data structure to hold the integral terms.
struct IntegralData
{
  IntegralData() : I(0), I_phi(0), dI(0), dI_phi(0) {}

  IntegralData(double i, double i_phi, Vector3d const & di, Vector3d const & di_phi)
  : I(i), I_phi(i_phi), dI(di), dI_phi(di_phi) {}

  IntegralData operator+(IntegralData const & other) const
  { return IntegralData(I + other.I, I_phi + other.I_phi, dI + other.dI, dI_phi + other.dI_phi); }

  IntegralData & operator+=(IntegralData const & other)
  {
    I  += other.I;  I_phi  += other.I_phi;
    dI += other.dI; dI_phi += other.dI_phi;

    return *this;
  }

  IntegralData & operator*=(double s)
  {
    I  *= s; I_phi  *= s;
    dI *= s; dI_phi *= s;

    return *this;
  }

  double I, I_phi;
  Vector3d dI, dI_phi;

}; // struct IntegralData

} // namespace ImlsSurfaceInternal

/**
 * Use the Implicit Moving Least Squares (IMLS) method to compute a smooth approximation to an input mesh. We use the algorithm
 * of:
 *   Chen Shen, James O'Brien and Jonathan R. Shewchuk, "Interpolating and Approximating Implicit Surfaces from Polygon Soup".
 *   Proc. SIGGRAPH, pp. 896-904 (2004).
 *
 * The result is an implicit function whose zero level set (i.e. the set of points at which the function assumes the value zero)
 * gives the desired surface.
 *
 * Our code is derived from the Brown Indexed Mesh Library
 * (http://public.kitware.com/vxl/doc/release/contrib/brl/bbas/imesh/html/index.html), part of the VXL project
 * (http://vxl.sourceforge.net).
 */
class THEA_API ImlsSurface : private Noncopyable
{
  public:
    /**
     * Construct an IMLS surface from an input mesh.
     *
     * @param mesh The input mesh, treated as a polygon soup.
     * @param eps_ Controls smoothness of the generated isosurface (negative for default value).
     * @param lambda_ Controls accuracy of computations (negative for default value).
     * @param enforce_bounded If true, a loose bounding box of the input mesh is added to the polygon soup to ensure the
     *   generated surface is bounded.
     */
    template <typename MeshT>
    ImlsSurface(MeshT const & mesh, double eps_ = -1, double lambda_ = -1, bool enforce_bounded = false);

    /**
     * Construct an IMLS surface from an input mesh group.
     *
     * @param mg The input mesh group, treated as a polygon soup.
     * @param eps_ Controls smoothness of the generated isosurface.
     * @param lambda_ Controls accuracy of computations.
     * @param enforce_bounded If true, a loose bounding box of the input mesh is added to the polygon soup to ensure the
     *   generated surface is bounded.
     */
    template <typename MeshT>
    ImlsSurface(Graphics::MeshGroup<MeshT> const & mg, double eps_ = -1, double lambda_ = -1, bool enforce_bounded = false);

    /** Evaluate the implicit surface at a point represented by a Vector3. */
    double operator()(Vector3 const & p) const;

    /** Evaluate the implicit surface at a point. PointT must have public member functions x(), y() and z(). */
    template <typename PointT>
    double operator()(PointT const & p) const { return operator()(Vector3((Real)p.x(), (Real)p.y(), (Real)p.z())); }

    /** Get a bounding box for the original input mesh. */
    AxisAlignedBox3 const & getInputBounds() const;

    /** Evaluate the function and its derivative (returned by reference). */
    double deriv(Vector3 const & p, Vector3d & dp) const;

    /** Evaluate the function and its first and second derivatives (returned by reference). */
    double deriv2(Vector3 const & p, Vector3d & dp, Matrix3d & ddp) const;

    /** Get the smoothness of the surface (epsilon). */
    double getSmoothness() const { return eps; }

    /** Set the smoothness of the surface (epsilon). Pass a negative argument to choose the default value. */
    void setSmoothness(double eps_);

    /** Get the accuracy parameter (lambda). A larger value increases speed but reduces computation accuracy. */
    double getAccuracy() const { return lambda; }

    /**
     * Set the accuracy parameter (lambda). A larger value increases speed but reduces computation accuracy. Pass a negative
     * argument to choose the default value.
     */
    void setAccuracy(double lambda_);

    /**
     * Get a bounding ball for the surface. The function will have a negative value at its center. An exception will be thrown
     * if no such ball can be found.
     *
     * @todo Implement actual search for negative center
     */
    Ball3 getBoundingBallWithNegativeCenter() const
    {
      Vector3 center = getInputBounds().getCenter();  // will do for now, but implement actual search later
      if ((*this)(center) >= 0)
        throw Error("ImlsSurface: Couldn't find negative value center for bounding ball");

      Real radius = getInputBounds().getExtent().norm();  // purposely leave some slack in case the surface deviates
                                                          // significantly from the input mesh
      return Ball3(center, radius);
    }

  private:
    typedef ImlsSurfaceInternal::IndexedVertexTriple  IndexedVertexTriple;
    typedef ImlsSurfaceInternal::IndexedTriangle      IndexedTriangle;
    typedef ImlsSurfaceInternal::NodeAttribute        NodeAttribute;
    typedef ImlsSurfaceInternal::TriangleKdTree       TriangleKdTree;
    typedef ImlsSurfaceInternal::IntegralData         IntegralData;

    /** Functor prototype. */
    struct Functor
    {
      /** Destructor. */
      virtual ~Functor() {}

      /** Process a triangle. */
      virtual void evalTri(Vector3 const & p, IndexedTriangle const & tri) = 0;

      /** Process a node. */
      virtual void evalNode(Vector3 const & p, TriangleKdTree::Node const & node) = 0;

      /** Check if the error in a sum is within acceptable limits. */
      virtual bool acceptableError(double err) const = 0;

    }; // struct Functor

    /** Functor for evaluating the function. */
    struct EvalFunctor : public Functor
    {
      ImlsSurface const & surf;
      double sum, sum_phi;

      EvalFunctor(ImlsSurface const & surf_) : surf(surf_), sum(0), sum_phi(0) {}

      void evalTri(Vector3 const & p, IndexedTriangle const & tri);
      void evalNode(Vector3 const & p, TriangleKdTree::Node const & node);
      bool acceptableError(double err) const;

    }; // struct EvalFunctor

    /** Functor for evaluating the function and its derivative. */
    struct DerivFunctor : public Functor
    {
      ImlsSurface const & surf;
      IntegralData sums;

      DerivFunctor(ImlsSurface const & surf_) : surf(surf_) {}

      void evalTri(Vector3 const & p, IndexedTriangle const & tri);
      void evalNode(Vector3 const & p, TriangleKdTree::Node const & node);
      bool acceptableError(double err) const;

    }; // struct DerivFunctor

    friend struct Functor;
    friend struct EvalFunctor;
    friend struct DerivFunctor;

    /** Construction template. */
    template <typename MeshT> void constructFromMesh(MeshT const & mesh);

    /** Add a general mesh to the polygon soup. \a tris is used to store the generated triangles. */
    template < typename MeshT, typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value, int >::type = 0 >
    void addMesh(MeshT const & mesh, Array<IndexedTriangle> & tris);

    /** Add a DCEL mesh to the polygon soup. \a tris is used to store the generated triangles. */
    template < typename MeshT, typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value, int >::type = 0 >
    void addMesh(MeshT const & mesh, Array<IndexedTriangle> & tris);

    /** Add a mesh group to the polygon soup. \a tris is used to store the generated triangles. */
    template <typename MeshT>
    void addMeshGroup(Graphics::MeshGroup<MeshT> const & mg, Array<IndexedTriangle> & tris);

    /** Add a loose bounding box of all input meshes to the polygon soup, ensuring the isosurface is bounded. */
    void enforceBounds(Array<IndexedTriangle> & tris);

    /** Evaluate the weighting function given a squared distance. */
    double weight(double dist2) const
    {
      return 1.0 / (dist2 + eps2);
    }

    /** Evaluate the squared weighting function given a squared distance. */
    double weight2(double dist2) const
    {
      double w = weight(dist2);
      return w * w;
    }

    /** Recursively compute the area weighted centroids. */
    void computeCentroidsRec(TriangleKdTree::Node const * start);

    /** Recursively compute the unweighted integrals. */
    void computeUnweightedRec(TriangleKdTree::Node const * start);

    /** Compute the iso value such that the mean value at the vertices is zero. */
    void computeIsolevel();

    /** Adjust the phi values until all vertices are within the isosurface. Also computes the isolevel. */
    void computeEnclosingPhi();

    /** Apply a functor to kd-tree nodes, to approximate the implicit function or its derivative. */
    void eval(Vector3 const & p, Functor & functor) const;

    /**
     * Apply a functor to kd-tree nodes, to approximate the implicit function or its derivative. The tree is traversed
     * recursively, without using a priority queue.
     */
    void evalRec(Vector3 const & p, TriangleKdTree::Node const * start, Functor & functor) const;

    Array<Vector3>  verts;   ///< Mesh vertices.
    Array<double>   phi;     ///< Phi values assigned per vertex.
    TriangleKdTree  kdtree;  ///< KD-tree of mesh triangles.

    double  mesh_size;  ///< Size of the input mesh, measured as the diagonal of its bounding box.
    double  eps;        ///< Smoothness.
    double  eps2;       ///< Square of smoothness parameter.
    double  lambda;     ///< Accuracy parameter.
    double  isolevel;   ///< Mean value of function at vertices
    bool    bounded;    ///< Ensure the isosurface is bounded?

    static double const DEFAULT_SMOOTHNESS;
    static double const DEFAULT_ACCURACY;
    static intx   const DEFAULT_MAX_TRIS_PER_LEAF;

}; // class ImlsSurface

template <typename MeshT>
ImlsSurface::ImlsSurface(MeshT const & mesh, double eps_, double lambda_, bool enforce_bounded)
: isolevel(0), bounded(enforce_bounded)
{
  setSmoothness(eps_);
  setAccuracy(lambda_);
  constructFromMesh(mesh);
}

template <typename MeshT>
ImlsSurface::ImlsSurface(Graphics::MeshGroup<MeshT> const & mg, double eps_, double lambda_, bool enforce_bounded)
: isolevel(0), bounded(enforce_bounded)
{
  setSmoothness(eps_);
  setAccuracy(lambda_);

  Array<IndexedTriangle> tris;
  addMeshGroup(mg, tris);

  if (bounded)
    enforceBounds(tris);

  kdtree.init(tris.begin(), tris.end(), -1, DEFAULT_MAX_TRIS_PER_LEAF, true);
  mesh_size = kdtree.getBounds().getExtent().norm();

  computeCentroidsRec(kdtree.getRoot());
  computeUnweightedRec(kdtree.getRoot());
  computeEnclosingPhi();
}

template <typename MeshT>
void
ImlsSurface::constructFromMesh(MeshT const & mesh)
{
  Array<IndexedTriangle> tris;
  addMesh(mesh, tris);

  if (bounded)
    enforceBounds(tris);

  kdtree.init(tris.begin(), tris.end(), -1, -1, true);
  mesh_size = kdtree.getBounds().getExtent().norm();

  computeCentroidsRec(kdtree.getRoot());
  computeUnweightedRec(kdtree.getRoot());
  computeEnclosingPhi();
}

template < typename MeshT, typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value, int >::type >
void
ImlsSurface::addMesh(MeshT const & mesh, Array<IndexedTriangle> & tris)
{
  typedef UnorderedMap<typename MeshT::Vertex const *, int> VertexIndexMap;
  VertexIndexMap vertex_indices;

  size_t base = verts.size();
  size_t new_size = verts.size() + static_cast<size_t>(mesh.numVertices());
  verts.resize(new_size);
  phi.resize(new_size);

  size_t i = base;
  for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++i)
  {
    verts[i]  =  vi->getPosition();
    phi[i]    =  0;
    vertex_indices[&(*vi)] = i;
  }

  size_t i0, i1, i2, i3;
  Polygon3 poly;
  Array<intx> tri_indices;
  for (typename MeshT::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    if (fi->isTriangle())
    {
      typename MeshT::Face::VertexConstIterator vi = fi->verticesBegin();
      i0 = vertex_indices[&(**vi)]; ++vi;
      i1 = vertex_indices[&(**vi)]; ++vi;
      i2 = vertex_indices[&(**vi)];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
    }
    else if (fi->isQuad())
    {
      // FIXME: Handle non-convex quads

      typename MeshT::Face::VertexConstIterator vi = fi->verticesBegin();
      i0 = vertex_indices[&(**vi)]; ++vi;
      i1 = vertex_indices[&(**vi)]; ++vi;
      i2 = vertex_indices[&(**vi)]; ++vi;
      i3 = vertex_indices[&(**vi)];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i2, i3)));
    }
    else
    {
      poly.clear();
      for (typename MeshT::Face::VertexConstIterator vi = fi->verticesBegin(); vi != fi->verticesEnd(); ++vi)
        poly.addVertex((*vi)->getPosition(), static_cast<intx>(vertex_indices[&(**vi)]));

#ifdef THEA_DEBUG_BUILD
      intx num_tris = poly.triangulate(tri_indices);
      debugAssertM(tri_indices.size() == static_cast<size_t>(3 * num_tris),
                   "ImlsSurface: MeshT face triangulation error");
#else
      poly.triangulate(tri_indices);
#endif

      for (size_t i = 0; i < tri_indices.size(); i += 3)
        tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts,
                                       static_cast<size_t>(tri_indices[i]),
                                       static_cast<size_t>(tri_indices[i + 1]),
                                       static_cast<size_t>(tri_indices[i + 2]))));
    }
  }
}

template < typename MeshT, typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value, int >::type >
void
ImlsSurface::addMesh(MeshT const & mesh, Array<IndexedTriangle> & tris)
{
  typedef UnorderedMap<typename MeshT::Vertex const *, int> VertexIndexMap;
  VertexIndexMap vertex_indices;

  size_t base = verts.size();
  size_t new_size = verts.size() + static_cast<size_t>(mesh.numVertices());
  verts.resize(new_size);
  phi.resize(new_size);

  size_t i = base;
  for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++i)
  {
    verts[i]  =  vi->getPosition();
    phi[i]    =  0;
    vertex_indices[&(*vi)] = i;
  }

  size_t i0, i1, i2, i3;
  Polygon3 poly;
  Array<intx> tri_indices;
  for (typename MeshT::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    if (fi->isTriangle())
    {
      typename MeshT::Halfedge const * he = fi->getHalfedge();
      i0 = vertex_indices[he->getOrigin()]; he = he->next();
      i1 = vertex_indices[he->getOrigin()]; he = he->next();
      i2 = vertex_indices[he->getOrigin()];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
    }
    else if (fi->isQuad())
    {
      // FIXME: Handle non-convex quads

      typename MeshT::Halfedge const * he = fi->getHalfedge();
      i0 = vertex_indices[he->getOrigin()]; he = he->next();
      i1 = vertex_indices[he->getOrigin()]; he = he->next();
      i2 = vertex_indices[he->getOrigin()]; he = he->next();
      i3 = vertex_indices[he->getOrigin()];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i2, i3)));
    }
    else if (fi->numVertices() > 0)
    {
      poly.clear();
      typename MeshT::Halfedge const * he = fi->getHalfedge();
      do
      {
        poly.addVertex(he->getOrigin()->getPosition(), static_cast<intx>(vertex_indices[he->getOrigin()]));
        he = he->next();

      } while (he != fi->getHalfedge());

#ifdef THEA_DEBUG_BUILD
      intx num_tris = poly.triangulate(tri_indices);
      debugAssertM(tri_indices.size() == static_cast<size_t>(3 * num_tris), "ImlsSurface: Mesh face triangulation error");
#else
      poly.triangulate(tri_indices);
#endif

      for (size_t i = 0; i < tri_indices.size(); i += 3)
        tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts,
                                       static_cast<size_t>(tri_indices[i]),
                                       static_cast<size_t>(tri_indices[i + 1]),
                                       static_cast<size_t>(tri_indices[i + 2]))));
    }
  }
}

template <typename MeshT>
void
ImlsSurface::addMeshGroup(Graphics::MeshGroup<MeshT> const & mg, Array<IndexedTriangle> & tris)
{
  for (typename Graphics::MeshGroup<MeshT>::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
    addMesh(**mi, tris);

  for (typename Graphics::MeshGroup<MeshT>::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
    addMeshGroup(**ci, tris);
}

} // namespace Algorithms
} // namespace Thea

#endif
