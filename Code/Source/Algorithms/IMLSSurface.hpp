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

#ifndef __Thea_Algorithms_IMLSSurface_hpp__
#define __Thea_Algorithms_IMLSSurface_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Ball3.hpp"
#include "../MatrixMN.hpp"
#include "../Noncopyable.hpp"
#include "../Polygon3.hpp"
#include "../Triangle3.hpp"
#include "../UnorderedMap.hpp"
#include "../VectorN.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include "KDTreeN.hpp"
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

namespace IMLSSurfaceInternal {

typedef VectorN<2, double>      DoubleVector2;  // Double precision 2-vector.
typedef VectorN<3, double>      DoubleVector3;  // Double precision 3-vector.
typedef MatrixMN<3, 3, double>  DoubleMatrix3;  // Double precision 3x3 matrix.

// Three indexed vertices plus base position array.
class IndexedVertexTriple
{
  private:
    TheaArray<Vector3> const * array;
    array_size_t indices[3];

  public:
    /** Default constructor. */
    IndexedVertexTriple() {}

    /** Initializing constructor. */
    IndexedVertexTriple(TheaArray<Vector3> const * array_, array_size_t i0, array_size_t i1, array_size_t i2)
    {
      array = array_;
      indices[0] = i0;
      indices[1] = i1;
      indices[2] = i2;
    }

    /** Get the i'th index. */
    array_size_t getIndex(int i) const { return indices[i]; }

    /** Get the i'th vertex. */
    Vector3 const & getVertex(int i) const { return (*array)[indices[i]]; }

}; // IndexedVertexTriple

typedef Triangle3<IndexedVertexTriple> IndexedTriangle;  ///< Triangle with indexed vertices.

// KD-tree node attributes.
struct NodeAttribute
{
  double         area;        ///< Surface area of the contents of the node
  double         unweighted;  ///< Unweighted integrals of phi over the surface contained by the node.
  DoubleVector3  centroid;    ///< Area weighted centroid of the node.
  DoubleVector3  normal;      ///< Area weighted normal of the node.
  double         normal_len;  ///< Cached magnitude of the normal.
};

typedef KDTreeN<IndexedTriangle, 3, Real, NodeAttribute> TriangleKDTree;

// A data structure to hold the integral terms.
struct IntegralData
{
  IntegralData() : I(0), I_phi(0), dI(0), dI_phi(0) {}

  IntegralData(double i, double i_phi, DoubleVector3 const & di, DoubleVector3 const & di_phi)
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
  DoubleVector3 dI, dI_phi;

}; // struct IntegralData

} // namespace IMLSSurfaceInternal

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
class THEA_API IMLSSurface : private Noncopyable
{
  public:
    typedef IMLSSurfaceInternal::DoubleVector2  DoubleVector2;  ///< Double precision 2-vector.
    typedef IMLSSurfaceInternal::DoubleVector3  DoubleVector3;  ///< Double precision 3-vector.
    typedef IMLSSurfaceInternal::DoubleMatrix3  DoubleMatrix3;  ///< Double precision 3x3 matrix.

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
    IMLSSurface(MeshT const & mesh, double eps_ = -1, double lambda_ = -1, bool enforce_bounded = false);

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
    IMLSSurface(Graphics::MeshGroup<MeshT> const & mg, double eps_ = -1, double lambda_ = -1, bool enforce_bounded = false);

    /** Evaluate the implicit surface at a point represented by a Vector3. */
    double operator()(Vector3 const & p) const;

    /** Evaluate the implicit surface at a point. PointT must have public member functions x(), y() and z(). */
    template <typename PointT>
    double operator()(PointT const & p) const { return operator()(Vector3((Real)p.x(), (Real)p.y(), (Real)p.z())); }

    /** Get a bounding box for the original input mesh. */
    AxisAlignedBox3 const & getInputBounds() const;

    /** Evaluate the function and its derivative (returned by reference). */
    double deriv(Vector3 const & p, DoubleVector3 & dp) const;

    /** Evaluate the function and its first and second derivatives (returned by reference). */
    double deriv2(Vector3 const & p, DoubleVector3 & dp, DoubleMatrix3 & ddp) const;

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
        throw Error("IMLSSurface: Couldn't find negative value center for bounding ball");

      Real radius = getInputBounds().getExtent().length();  // purposely leave some slack in case the surface deviates
                                                         // significantly from the input mesh
      return Ball3(center, radius);
    }

  private:
    typedef IMLSSurfaceInternal::IndexedVertexTriple  IndexedVertexTriple;
    typedef IMLSSurfaceInternal::IndexedTriangle      IndexedTriangle;
    typedef IMLSSurfaceInternal::NodeAttribute        NodeAttribute;
    typedef IMLSSurfaceInternal::TriangleKDTree       TriangleKDTree;
    typedef IMLSSurfaceInternal::IntegralData         IntegralData;

    /** Functor prototype. */
    struct Functor
    {
      /** Destructor. */
      virtual ~Functor() {}

      /** Process a triangle. */
      virtual void evalTri(Vector3 const & p, IndexedTriangle const & tri) = 0;

      /** Process a node. */
      virtual void evalNode(Vector3 const & p, TriangleKDTree::Node const & node) = 0;

      /** Check if the error in a sum is within acceptable limits. */
      virtual bool acceptableError(double err) const = 0;

    }; // struct Functor

    /** Functor for evaluating the function. */
    struct EvalFunctor : public Functor
    {
      IMLSSurface const & surf;
      double sum, sum_phi;

      EvalFunctor(IMLSSurface const & surf_) : surf(surf_), sum(0), sum_phi(0) {}

      void evalTri(Vector3 const & p, IndexedTriangle const & tri);
      void evalNode(Vector3 const & p, TriangleKDTree::Node const & node);
      bool acceptableError(double err) const;

    }; // struct EvalFunctor

    /** Functor for evaluating the function and its derivative. */
    struct DerivFunctor : public Functor
    {
      IMLSSurface const & surf;
      IntegralData sums;

      DerivFunctor(IMLSSurface const & surf_) : surf(surf_) {}

      void evalTri(Vector3 const & p, IndexedTriangle const & tri);
      void evalNode(Vector3 const & p, TriangleKDTree::Node const & node);
      bool acceptableError(double err) const;

    }; // struct DerivFunctor

    friend struct Functor;
    friend struct EvalFunctor;
    friend struct DerivFunctor;

    /** Construction template. */
    template <typename MeshT> void constructFromMesh(MeshT const & mesh);

    /** Add a general mesh to the polygon soup. \a tris is used to store the generated triangles. */
    template <typename MeshT>
    void addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                 typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy = NULL);

    /** Add a DCEL mesh to the polygon soup. \a tris is used to store the generated triangles. */
    template <typename MeshT>
    void addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                 typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy = NULL);

    /** Add a CGAL mesh to the polygon soup. \a tris is used to store the generated triangles. */
    template <typename MeshT>
    void addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                 typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type * dummy = NULL);

    /** Add a mesh group to the polygon soup. \a tris is used to store the generated triangles. */
    template <typename MeshT>
    void addMeshGroup(Graphics::MeshGroup<MeshT> const & mg, TheaArray<IndexedTriangle> & tris);

    /** Add a loose bounding box of all input meshes to the polygon soup, ensuring the isosurface is bounded. */
    void enforceBounds(TheaArray<IndexedTriangle> & tris);

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
    void computeCentroidsRec(TriangleKDTree::Node const * start);

    /** Recursively compute the unweighted integrals. */
    void computeUnweightedRec(TriangleKDTree::Node const * start);

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
    void evalRec(Vector3 const & p, TriangleKDTree::Node const * start, Functor & functor) const;

    TheaArray<Vector3>  verts;   ///< Mesh vertices.
    TheaArray<double>   phi;     ///< Phi values assigned per vertex.
    TriangleKDTree      kdtree;  ///< KD-tree of mesh triangles.

    double  mesh_size;  ///< Size of the input mesh, measured as the diagonal of its bounding box.
    double  eps;        ///< Smoothness.
    double  eps2;       ///< Square of smoothness parameter.
    double  lambda;     ///< Accuracy parameter.
    double  isolevel;   ///< Mean value of function at vertices
    bool    bounded;    ///< Ensure the isosurface is bounded?

    static double const DEFAULT_SMOOTHNESS;
    static double const DEFAULT_ACCURACY;
    static long   const DEFAULT_MAX_TRIS_PER_LEAF;

}; // class IMLSSurface

template <typename MeshT>
IMLSSurface::IMLSSurface(MeshT const & mesh, double eps_, double lambda_, bool enforce_bounded)
: isolevel(0), bounded(enforce_bounded)
{
  setSmoothness(eps_);
  setAccuracy(lambda_);
  constructFromMesh(mesh);
}

template <typename MeshT>
IMLSSurface::IMLSSurface(Graphics::MeshGroup<MeshT> const & mg, double eps_, double lambda_, bool enforce_bounded)
: isolevel(0), bounded(enforce_bounded)
{
  setSmoothness(eps_);
  setAccuracy(lambda_);

  TheaArray<IndexedTriangle> tris;
  addMeshGroup(mg, tris);

  if (bounded)
    enforceBounds(tris);

  kdtree.init(tris.begin(), tris.end(), -1, DEFAULT_MAX_TRIS_PER_LEAF, true);
  mesh_size = kdtree.getBounds().getExtent().length();

  computeCentroidsRec(kdtree.getRoot());
  computeUnweightedRec(kdtree.getRoot());
  computeEnclosingPhi();
}

template <typename MeshT>
void
IMLSSurface::constructFromMesh(MeshT const & mesh)
{
  TheaArray<IndexedTriangle> tris;
  addMesh(mesh, tris);

  if (bounded)
    enforceBounds(tris);

  kdtree.init(tris.begin(), tris.end(), -1, -1, true);
  mesh_size = kdtree.getBounds().getExtent().length();

  computeCentroidsRec(kdtree.getRoot());
  computeUnweightedRec(kdtree.getRoot());
  computeEnclosingPhi();
}

template <typename MeshT>
void
IMLSSurface::addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                     typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy)
{
  typedef TheaUnorderedMap<typename MeshT::Vertex const *, int> VertexIndexMap;
  VertexIndexMap vertex_indices;

  array_size_t base = verts.size();
  array_size_t new_size = verts.size() + static_cast<array_size_t>(mesh.numVertices());
  verts.resize(new_size);
  phi.resize(new_size);

  array_size_t i = base;
  for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++i)
  {
    verts[i]  =  vi->getPosition();
    phi[i]    =  0;
    vertex_indices[&(*vi)] = i;
  }

  array_size_t i0, i1, i2, i3;
  Polygon3 poly;
  TheaArray<long> tri_indices;
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
        poly.addVertex((*vi)->getPosition(), static_cast<long>(vertex_indices[&(**vi)]));

#ifdef THEA_DEBUG_BUILD
      long num_tris = poly.triangulate(tri_indices);
      debugAssertM(tri_indices.size() == static_cast<array_size_t>(3 * num_tris),
                   "IMLSSurface: MeshT face triangulation error");
#else
      poly.triangulate(tri_indices);
#endif

      for (array_size_t i = 0; i < tri_indices.size(); i += 3)
        tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts,
                                       static_cast<array_size_t>(tri_indices[i]),
                                       static_cast<array_size_t>(tri_indices[i + 1]),
                                       static_cast<array_size_t>(tri_indices[i + 2]))));
    }
  }
}

template <typename MeshT>
void
IMLSSurface::addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                     typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy)
{
  typedef TheaUnorderedMap<typename MeshT::Vertex const *, int> VertexIndexMap;
  VertexIndexMap vertex_indices;

  array_size_t base = verts.size();
  array_size_t new_size = verts.size() + static_cast<array_size_t>(mesh.numVertices());
  verts.resize(new_size);
  phi.resize(new_size);

  array_size_t i = base;
  for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++i)
  {
    verts[i]  =  vi->getPosition();
    phi[i]    =  0;
    vertex_indices[&(*vi)] = i;
  }

  array_size_t i0, i1, i2, i3;
  Polygon3 poly;
  TheaArray<long> tri_indices;
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
        poly.addVertex(he->getOrigin()->getPosition(), static_cast<long>(vertex_indices[he->getOrigin()]));
        he = he->next();

      } while (he != fi->getHalfedge());

#ifdef THEA_DEBUG_BUILD
      long num_tris = poly.triangulate(tri_indices);
      debugAssertM(tri_indices.size() == static_cast<array_size_t>(3 * num_tris), "IMLSSurface: Mesh face triangulation error");
#else
      poly.triangulate(tri_indices);
#endif

      for (array_size_t i = 0; i < tri_indices.size(); i += 3)
        tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts,
                                       static_cast<array_size_t>(tri_indices[i]),
                                       static_cast<array_size_t>(tri_indices[i + 1]),
                                       static_cast<array_size_t>(tri_indices[i + 2]))));
    }
  }
}

template <typename MeshT>
void
IMLSSurface::addMesh(MeshT const & mesh, TheaArray<IndexedTriangle> & tris,
                     typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type * dummy)
{
  typedef TheaUnorderedMap<typename MeshT::Vertex const *, int> VertexIndexMap;
  VertexIndexMap vertex_indices;

  array_size_t base = verts.size();
  array_size_t new_size = verts.size() + static_cast<array_size_t>(mesh.size_of_vertices());
  verts.resize(new_size);
  phi.resize(new_size);

  array_size_t i = base;
  for (typename MeshT::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi, ++i)
  {
    verts[i]  =  vi->point();
    phi[i]    =  0;
    vertex_indices[&(*vi)] = i;
  }

  array_size_t i0, i1, i2, i3;
  Polygon3 poly;
  TheaArray<long> tri_indices;
  for (typename MeshT::Facet_const_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi)
  {
    if (fi->is_triangle())
    {
      typename MeshT::Facet::Halfedge_around_facet_const_circulator hc = fi->facet_begin();
      i0 = vertex_indices[&(*(hc->vertex()))]; ++hc;
      i1 = vertex_indices[&(*(hc->vertex()))]; ++hc;
      i2 = vertex_indices[&(*(hc->vertex()))];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
    }
    else if (fi->is_quad())
    {
      // FIXME: Handle non-convex quads

      typename MeshT::Facet::Halfedge_around_facet_const_circulator hc = fi->facet_begin();
      i0 = vertex_indices[&(*(hc->vertex()))]; ++hc;
      i1 = vertex_indices[&(*(hc->vertex()))]; ++hc;
      i2 = vertex_indices[&(*(hc->vertex()))]; ++hc;
      i3 = vertex_indices[&(*(hc->vertex()))];

      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i1, i2)));
      tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, i0, i2, i3)));
    }
    else if (fi->facet_degree() > 0)
    {
      poly.clear();
      typename MeshT::Facet::Halfedge_around_facet_const_circulator hc = fi->facet_begin();
      do
      {
        poly.addVertex(hc->vertex()->point(), static_cast<long>(vertex_indices[&(*(hc->vertex()))]));

      } while (++hc != fi->facet_begin());

#ifdef THEA_DEBUG_BUILD
      long num_tris = poly.triangulate(tri_indices);
      debugAssertM(tri_indices.size() == static_cast<array_size_t>(3 * num_tris), "IMLSSurface: Mesh face triangulation error");
#else
      poly.triangulate(tri_indices);
#endif

      for (array_size_t i = 0; i < tri_indices.size(); i += 3)
        tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts,
                                       static_cast<array_size_t>(tri_indices[i]),
                                       static_cast<array_size_t>(tri_indices[i + 1]),
                                       static_cast<array_size_t>(tri_indices[i + 2]))));
    }
  }
}

template <typename MeshT>
void
IMLSSurface::addMeshGroup(Graphics::MeshGroup<MeshT> const & mg, TheaArray<IndexedTriangle> & tris)
{
  for (typename Graphics::MeshGroup<MeshT>::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
    addMesh(**mi, tris);

  for (typename Graphics::MeshGroup<MeshT>::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
    addMeshGroup(**ci, tris);
}

} // namespace Algorithms
} // namespace Thea

#endif
