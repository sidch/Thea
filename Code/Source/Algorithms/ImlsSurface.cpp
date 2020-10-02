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

#include "ImlsSurface.hpp"
#include "../Math.hpp"
#include <algorithm>
#include <functional>
#include <limits>

#define IMLS_SURFACE_EVAL_RECURSIVE
// #define IMLS_SURFACE_EVAL_TIMER

#ifdef IMLS_SURFACE_EVAL_TIMER
#  include "../Stopwatch.hpp"
#endif

namespace Thea {
namespace Algorithms {

double const ImlsSurface::DEFAULT_SMOOTHNESS         =  0.01;
double const ImlsSurface::DEFAULT_ACCURACY           =  0.1;
intx   const ImlsSurface::DEFAULT_MAX_TRIS_PER_LEAF  =  3;

namespace ImlsSurfaceInternal {

// A kd-tree node and an associated scalar value.
struct NodeValuePair
{
  NodeValuePair() {}
  NodeValuePair(TriangleKdTree::Node const * node_, double value_) : node(node_), value(value_) {}

  // Susceptible to floating point error?
  bool operator<(NodeValuePair const & rhs) const { return value < rhs.value; }
  bool operator>(NodeValuePair const & rhs) const { return value > rhs.value; }

  TriangleKdTree::Node const * node;
  double value;

}; // struct NodeValuePair

// Traverse the tree looking for leaf nodes that contain the query. Returns a vector of leaf nodes paired with min square
// distance to the node's bounding box. Also returns a vector of the unexplored internal node-distance^2 pairs. Resulting
// vectors are unsorted.
void traverseKdTree(TriangleKdTree const & kdtree, Vector3 const & query, Array<NodeValuePair> & leaves,
                    Array<NodeValuePair> & internal_nodes);

// Compute the minimum squared distance between a query point and the triangles in a kd-tree node.
double triSquaredDistance(TriangleKdTree const & kdtree, TriangleKdTree::Node const & node, Vector3 const & query);

// Find the leaf node closest to a query point (w.r.t. the triangles in the node), and queue all explored nodes.
TriangleKdTree::Node const * closestLeaf(TriangleKdTree const & kdtree, Vector3 const & query,
                                         Array<NodeValuePair> & explored);

// Integrals of f(x) dx and x * f(x) dx over [0, 1] where f(x) = 1 / ((x + k1)^2 + k2)^2.
void lineIntegrals(double k1, double k2, double & I1, double & Ix);

// Integrals of f(x) dx and x * f(x) dx over [0, 1] where f(x) = 1 / ((x + k1)^2 + k2)^2. Also compute the integrals when
// f(x) = 1 / ((x + k1)^2 + k2)^3 (for use in derivatives).
void lineIntegrals(double k1, double k2, double & I1, double & Ix, double & dI1, double & dIx, double & dIx2);

// Line integral of the squared weight function times a linear value on the line from p0 to p1 (value at p0 is v0 and at p1
// is v1).
//
// @param eps2 epsilon^2.
double lineIntegral(Vector3 const & x, Vector3 const & p0, Vector3 const & p1, double v0, double v1, double eps2_);

// The derivative of the line integral with respect to \a x.
//
// @param eps2 epsilon^2.
Vector3d lineIntegralDeriv(Vector3 const & x, Vector3 const & p0, Vector3 const & p1, double v0, double v1, double eps2_);


// Area integral of the squared weight function times a linearly interpolated value. Call triangleQuadrature() to first
// split an arbitrary triangle.
//
// @param x Sample point.
// @param m Closest point on the triangle to sample point \a x.
// @param pp Second closest vertex.
// @param pm Furthest vertex.
// @param eps2 epsilon^2.
Vector2d splitTriangleQuadrature(Vector3 const & x, Vector3 const & pm, Vector3 const & p1, Vector3 const & p2, double vm,
                                      double v1, double v2, double eps2_);

// Area integral of the squared weight function times a linearly interpolated value. Also computes vector term used in the
// derivative
//
// @param x Sample point.
// @param m Closest point on the triangle to sample point \a x.
// @param pp Second closest vertex.
// @param pm Furthest vertex.
// @param eps2 epsilon^2.
IntegralData splitTriangleQuadratureWithDeriv(Vector3 const & x, Vector3 const & pm, Vector3 const & p1, Vector3 const & p2,
                                              double vm, double v1, double v2, double eps2_);

// Find the closest point on the triangle a, b, c to point p. The un-normalized normal vector (b - a) x (c - a) is
// precomputed and also passed in.
//
// @param sqdist is the square of the distance to the triangle (returned by reference)
// @param u Barycentric coordinate of the closest point (with v)
// @param v Barycentric coordinate of the closest point (with u)
//
// @return A code indicating that the closest point:
//  - 0 does not exist (should not occur)
//  - 1 is \a a
//  - 2 is \a b
//  - 3 is on the edge from \a a to \a b
//  - 4 is \a c
//  - 5 is on the edge from \a a to \a c
//  - 6 is on the edge from \a b to \a c
//  - 7 is on the face of the triangle
//
char triangleClosestPoint(Vector3 const & p, Vector3 const & a, Vector3 const & b, Vector3 const & c, Vector3 const & n,
                          double & sqdist, double & u, double & v);

// Area integral of the squared weight function times a linearly interpolated value.
//
// @param eps2 epsilon^2
template <typename T, typename F> T triangleQuadrature(F quad_func, Vector3 const & x, Vector3 const & p0, Vector3 const & p1,
                                                       Vector3 const & p2, Vector3 const & n, double v0, double v1, double v2,
                                                       double eps2_);

} // namespace ImlsSurfaceInternal

void
ImlsSurface::enforceBounds(Array<IndexedTriangle> & tris)
{
  // build enclosure
  AxisAlignedBox3 box;
  for (size_t i = 0; i < verts.size(); ++i)
    box.merge(verts[i]);

  box.scaleCentered(2);

  size_t base = verts.size();
  size_t new_size = verts.size() + 8;
  verts.resize(new_size);
  phi.resize(new_size);

  Vector3 const & lo = box.getLow();
  Vector3 const & hi = box.getHigh();

  verts[base    ] = Vector3(lo.x(), lo.y(), lo.z());
  verts[base + 1] = Vector3(lo.x(), lo.y(), hi.z());
  verts[base + 2] = Vector3(lo.x(), hi.y(), lo.z());
  verts[base + 3] = Vector3(lo.x(), hi.y(), hi.z());
  verts[base + 4] = Vector3(hi.x(), lo.y(), lo.z());
  verts[base + 5] = Vector3(hi.x(), lo.y(), hi.z());
  verts[base + 6] = Vector3(hi.x(), hi.y(), lo.z());
  verts[base + 7] = Vector3(hi.x(), hi.y(), hi.z());

  std::fill(&phi[base], &phi[base + 8], 1.0);

  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base    , base + 1, base + 2)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 1, base + 3, base + 2)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base    , base + 6, base + 4)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base    , base + 2, base + 6)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 2, base + 7, base + 6)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 2, base + 3, base + 7)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 3, base + 5, base + 7)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 3, base + 1, base + 5)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 1, base + 4, base + 5)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 1, base    , base + 4)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 4, base + 7, base + 5)));
  tris.push_back(IndexedTriangle(IndexedVertexTriple(&verts, base + 4, base + 6, base + 7)));
}

void
ImlsSurface::computeIsolevel()
{
  double mean = 0;
  size_t const num_verts = bounded ? verts.size() - 8: verts.size();

  if (num_verts > 0)
  {
    for (size_t i = 0; i < num_verts; ++i)
      mean += (*this)(verts[i]);

    isolevel = mean / num_verts;
  }
  else
    isolevel = 0;

  THEA_DEBUG << "ImlsSurface: Isolevel = " << isolevel;
}

void
ImlsSurface::computeEnclosingPhi()
{
  typedef std::pair<size_t, double> PairId;
  Array<PairId> outside;

  double mean = 0;
  size_t const num_verts = bounded ? verts.size() - 8: verts.size();

  if (num_verts > 0)
  {
    for (size_t i = 0; i < num_verts; ++i)
    {
      double val = (*this)(verts[i]);
      mean += val;
      if (val > 0)
        outside.push_back(PairId(i, val));
    }

    isolevel = mean / num_verts;
  }
  else
    isolevel = 0;

  THEA_DEBUG << "ImlsSurface: Isolevel = " << isolevel;

  for (double s = 1; !outside.empty(); s *= 2)
  {
    THEA_DEBUG << "ImlsSurface: " << outside.size() << " vertices outside";

    for (Array<PairId>::const_iterator oi = outside.begin(); oi != outside.end(); ++oi)
      phi[oi->first] -= s * oi->second;

    computeUnweightedRec(kdtree.getRoot());

    Array<PairId> next_outside;
    for (Array<PairId>::const_iterator oi = outside.begin(); oi != outside.end(); ++oi)
    {
      double val = (*this)(verts[oi->first]);
      if (val > std::fabs(std::numeric_limits<double>::epsilon() * phi[oi->first]))
        next_outside.push_back(PairId(oi->first, val));
    }

    outside.swap(next_outside);
  }
}

void
ImlsSurface::computeCentroidsRec(TriangleKdTree::Node const * start)
{
  using namespace ImlsSurfaceInternal;

  if (!start) return;

  NodeAttribute & my_attrib = const_cast<NodeAttribute &>(start->attr());  // yes, yes, it's ugly...

  if (start->isLeaf())
  {
    my_attrib.area      =  0;
    my_attrib.centroid  =  Vector3d::Zero();
    my_attrib.normal    =  Vector3d::Zero();

    for (TriangleKdTree::Node::ElementIndexConstIterator ei = start->elementIndicesBegin(); ei != start->elementIndicesEnd();
         ++ei)
    {
      IndexedTriangle const & tri = kdtree.getElements()[*ei];

      my_attrib.area      +=  tri.getArea();
      my_attrib.centroid  +=  tri.getArea() * tri.getCentroid().cast<double>();
      my_attrib.normal    +=  tri.getNormal().cast<double>();
    }

    if (my_attrib.area > 0) my_attrib.centroid /= my_attrib.area;
    my_attrib.normal_len = my_attrib.normal.norm();
  }
  else
  {
    TriangleKdTree::Node const * lo = start->getLowChild();
    TriangleKdTree::Node const * hi = start->getHighChild();

    computeCentroidsRec(lo);
    computeCentroidsRec(hi);

    NodeAttribute const & lo_attrib = lo->attr();
    NodeAttribute const & hi_attrib = hi->attr();

    my_attrib.area = lo_attrib.area + hi_attrib.area;
    if (my_attrib.area > 0)
      my_attrib.centroid = (lo_attrib.area * lo_attrib.centroid + hi_attrib.area * hi_attrib.centroid) / my_attrib.area;
    else
      my_attrib.centroid = Vector3d::Zero();

    my_attrib.normal      =  lo_attrib.normal + hi_attrib.normal;
    my_attrib.normal_len  =  my_attrib.normal.norm();
  }
}

void
ImlsSurface::computeUnweightedRec(TriangleKdTree::Node const * start)
{
  if (!start) return;

  NodeAttribute & my_attrib = const_cast<NodeAttribute &>(start->attr());  // yes, yes, it's ugly...

  if (start->isLeaf())
  {
    my_attrib.unweighted = 0;

    if (my_attrib.area > 0)
    {
      for (TriangleKdTree::Node::ElementIndexConstIterator ei = start->elementIndicesBegin(); ei != start->elementIndicesEnd();
           ++ei)
      {
        IndexedTriangle const & tri = kdtree.getElements()[*ei];
        my_attrib.unweighted += ((phi[tri.getVertices().getIndex(0)]
                                + phi[tri.getVertices().getIndex(1)]
                                + phi[tri.getVertices().getIndex(2)]) / 3) * tri.getArea();

        debugAssertM(!Math::isNaN(my_attrib.unweighted), "ImlsSurface: Unweighted value is NaN");
      }
    }
  }
  else
  {
    TriangleKdTree::Node const * lo = start->getLowChild();
    TriangleKdTree::Node const * hi = start->getHighChild();

    computeUnweightedRec(lo);
    computeUnweightedRec(hi);

    my_attrib.unweighted = lo->attr().unweighted + hi->attr().unweighted;
  }
}

AxisAlignedBox3 const &
ImlsSurface::getInputBounds() const
{
  return kdtree.getBounds();
}

void
ImlsSurface::setSmoothness(double eps_)
{
  eps = eps_ < 0 ? DEFAULT_SMOOTHNESS : eps_;
  eps2 = eps * eps;
  isolevel = 0;

  size_t const num_verts = bounded ? verts.size() - 8: verts.size();
  for (size_t i = 0;         i < num_verts;  ++i) phi[i] = 0.0;
  for (size_t i = num_verts; i < phi.size(); ++i) phi[i] = 1.0;

  computeUnweightedRec(kdtree.getRoot());
  computeEnclosingPhi();
}

void
ImlsSurface::setAccuracy(double lambda_)
{
  lambda = lambda_ < 0 ? DEFAULT_ACCURACY : lambda_;
}

void
ImlsSurface::eval(Vector3 const & p, Functor & functor) const
{
  using namespace ImlsSurfaceInternal;

  // Find the nodes that need to be evaluated
  Array<NodeValuePair> remain;
  closestLeaf(kdtree, p, remain);

  // Compute the (negative) maximum error of integration -- stored negative so that the max error is first when sorted by <
  for (Array<NodeValuePair>::iterator ri = remain.begin(); ri != remain.end(); ++ri)
  {
    double max_w2 = weight2(ri->value);
    double min_w2 = weight2(ri->node->getBounds().squaredMaxDistance(p));
    ri->value = (max_w2 - min_w2) * ri->node->attr().area;
  }

  // Heapify the set of nodes, by decreasing order of approximation error
  std::make_heap(remain.begin(), remain.end());

  // Now add up the contributions from the queued nodes
  IndexedTriangle const * tris = kdtree.getElements();
  while (!remain.empty() && !functor.acceptableError(remain.front().value))
  {
    TriangleKdTree::Node const * current = remain.front().node;
    std::pop_heap(remain.begin(), remain.end());
    remain.pop_back();

    if (current->isLeaf())
    {
      for (TriangleKdTree::Node::ElementIndexConstIterator ii = current->elementIndicesBegin();
           ii != current->elementIndicesEnd(); ++ii)
        functor.evalTri(p, tris[*ii]);
    }
    else
    {
      TriangleKdTree::Node const * lo = current->getLowChild();
      double max_w2 = weight2(lo->getBounds().squaredDistance(p));
      double min_w2 = weight2(lo->getBounds().squaredMaxDistance(p));
      remain.push_back(NodeValuePair(lo, (max_w2 - min_w2) * lo->attr().area));
      std::push_heap(remain.begin(), remain.end());

      TriangleKdTree::Node const * hi = current->getHighChild();
      max_w2 = weight2(hi->getBounds().squaredDistance(p));
      min_w2 = weight2(hi->getBounds().squaredMaxDistance(p));
      remain.push_back(NodeValuePair(hi, (max_w2 - min_w2) * hi->attr().area));
      std::push_heap(remain.begin(), remain.end());
    }
  }

  // Approximate the contribution from the remaining nodes
  for (Array<NodeValuePair>::iterator ri = remain.begin(); ri != remain.end(); ++ri)
    functor.evalNode(p, *ri->node);
}

void
ImlsSurface::evalRec(Vector3 const & p, TriangleKdTree::Node const * start, Functor & functor) const
{
  if (!start) return;

  if (start->isLeaf())
  {
    IndexedTriangle const * tris = kdtree.getElements();
    for (TriangleKdTree::Node::ElementIndexConstIterator ii = start->elementIndicesBegin(); ii != start->elementIndicesEnd();
         ++ii)
      functor.evalTri(p, tris[*ii]);
  }
  else
  {
    TriangleKdTree::Node const * c[2] = { start->getLowChild(), start->getHighChild() };

    // Make sure c[1] is the further child
    double d[2] = { c[0]->getBounds().squaredDistance(p), c[1]->getBounds().squaredDistance(p) };
    if (d[0] > d[1])
    {
      std::swap(c[0], c[1]);
      std::swap(d[0], d[1]);
    }

    // Nearer child first, since it's likely to have larger error
    for (int i = 0; i < 2; ++i)
    {
      double max_w2 = weight2(d[i]);
      double min_w2 = weight2(c[i]->getBounds().squaredMaxDistance(p));
      double err = (max_w2 - min_w2) * c[i]->attr().area;

      if (functor.acceptableError(err))
        functor.evalNode(p, *c[i]);
      else
        evalRec(p, c[i], functor);
    }
  }
}

double
ImlsSurface::operator()(Vector3 const & p) const
{
#ifdef IMLS_SURFACE_EVAL_TIMER
  static double total_time = 0;
  static intx num_calls = 1;
  Stopwatch timer;
  timer.tick();
#endif

  EvalFunctor efunc(*this);

#ifdef IMLS_SURFACE_EVAL_RECURSIVE
  evalRec(p, kdtree.getRoot(), efunc);
#else
  eval(p, efunc);
#endif

#if IMLS_SURFACE_EVAL_TIMER
  timer.tock();
  total_time += timer.elapsedTime();
  if (num_calls % 1000 == 0)
    THEA_CONSOLE << num_calls << " calls to operator(), average time " << 1000.0 * total_time / num_calls << " milliseconds";
  ++num_calls;
#endif

  debugAssertM(efunc.sum != 0, "ImlsSurface: Sum is zero");

  return efunc.sum_phi / efunc.sum - isolevel;
}

double
ImlsSurface::deriv(Vector3 const & p, Vector3d & dp) const
{
  DerivFunctor dfunc(*this);

#ifdef IMLS_SURFACE_EVAL_RECURSIVE
  evalRec(p, kdtree.getRoot(), dfunc);
#else
  eval(p, dfunc);
#endif

  debugAssertM(dfunc.sums.I != 0, "ImlsSurface: Sum is zero");

  dp = (-dfunc.sums.I_phi / (dfunc.sums.I * dfunc.sums.I)) * dfunc.sums.dI + dfunc.sums.dI_phi / dfunc.sums.I;
  return dfunc.sums.I_phi / dfunc.sums.I - isolevel;
}

void
ImlsSurface::EvalFunctor::evalTri(Vector3 const & p, IndexedTriangle const & tri)
{
  using namespace ImlsSurfaceInternal;

  // static intx num_calls = 1;
  // if ((num_calls++) % 1000 == 0) THEA_CONSOLE << num_calls << " calls to EvalFunctor::evalTri";

  size_t i0 = tri.getVertices().getIndex(0);
  size_t i1 = tri.getVertices().getIndex(1);
  size_t i2 = tri.getVertices().getIndex(2);

  Vector2d I = triangleQuadrature<Vector2d>(splitTriangleQuadrature, p, surf.verts[i0], surf.verts[i1],
                                                      surf.verts[i2], (Real)2 * tri.getNormal(), surf.phi[i0], surf.phi[i1],
                                                      surf.phi[i2], surf.eps2);

  debugAssertM(!Math::isNaN(I[0]) && !Math::isNaN(I[1]), "ImlsSurface: Integral is NaN");

  sum      +=  I[1];
  sum_phi  +=  (I[0] + tri.getNormal().dot(p - tri.getCentroid()) * I[1]);  // Triangle3 has unit normal
}

void
ImlsSurface::EvalFunctor::evalNode(Vector3 const & p, TriangleKdTree::Node const & node)
{
  using namespace ImlsSurfaceInternal;

  // static intx num_calls = 1;
  // if ((num_calls++) % 1000 == 0) THEA_CONSOLE << num_calls << " calls to EvalFunctor::evalNode";

  NodeAttribute const & attrib = node.attr();
  Vector3d v = p.cast<double>() - attrib.centroid;
  double w2 = surf.weight2(v.squaredNorm());

  sum      +=  w2 * attrib.area;
  sum_phi  +=  w2 * (attrib.unweighted + attrib.normal.dot(v));
}

bool
ImlsSurface::EvalFunctor::acceptableError(double err) const
{
  return err <= surf.lambda * sum;
}

void
ImlsSurface::DerivFunctor::evalTri(Vector3 const & p, IndexedTriangle const & tri)
{
  using namespace ImlsSurfaceInternal;

  size_t i0 = tri.getVertices().getIndex(0);
  size_t i1 = tri.getVertices().getIndex(1);
  size_t i2 = tri.getVertices().getIndex(2);

  IntegralData Id = triangleQuadrature<IntegralData>(splitTriangleQuadratureWithDeriv, p, surf.verts[i0], surf.verts[i1],
                                                     surf.verts[i2], (Real)2 * tri.getNormal(), surf.phi[i0], surf.phi[i1],
                                                     surf.phi[i2], surf.eps2);

  /* debugAssertM */ alwaysAssertM(!Math::isNaN(Id.I) && !Math::isNaN(Id.I_phi), "ImlsSurface: Integral is NaN");

  double plane_dist = tri.getNormal().dot(p - tri.getCentroid());

  sums.I       +=  Id.I;
  sums.I_phi   +=  Id.I_phi + plane_dist * Id.I;
  sums.dI      +=  Id.dI;
  sums.dI_phi  +=  Id.dI_phi + (Id.I * tri.getNormal().cast<double>() + plane_dist * Id.dI);
}

void
ImlsSurface::DerivFunctor::evalNode(Vector3 const & p, TriangleKdTree::Node const & node)
{
  using namespace ImlsSurfaceInternal;

  NodeAttribute const & attrib = node.attr();
  Vector3d v(p.cast<double>() - attrib.centroid);
  double w = surf.weight(v.squaredNorm());
  Vector3d dw(-4 * w * v);
  w   *=  w;
  dw  *=  w;

  double plane_dist = attrib.normal.dot(v);
  sums.I       +=  w * attrib.area;
  sums.I_phi   +=  w * attrib.unweighted + plane_dist * w;
  sums.dI      +=  dw * attrib.area;
  sums.dI_phi  +=  attrib.normal * w + dw * (attrib.unweighted + plane_dist);
}

bool
ImlsSurface::DerivFunctor::acceptableError(double err) const
{
  return err <= surf.lambda * sums.I;
}

double
ImlsSurface::deriv2(Vector3 const & p, Vector3d & dp, Matrix3d & ddp) const
{
  Real e = 1.0e-5f * mesh_size;
  double val = deriv(p, dp);

  Vector3d dpx, dpy, dpz;
  deriv(p.cwiseProduct(Vector3(1 + e, 1,     1    )), dpx);
  deriv(p.cwiseProduct(Vector3(1,     1 + e, 1    )), dpy);
  deriv(p.cwiseProduct(Vector3(1,     1,     1 + e)), dpz);

  dpx -= dp;
  dpx /= e;
  dpy -= dp;
  dpy /= e;
  dpz -= dp;
  dpz /= e;

  ddp(0, 0) = dpx[0]; ddp(0, 1) = dpy[0]; ddp(0, 2) = dpz[0];
  ddp(1, 0) = dpx[1]; ddp(1, 1) = dpy[1]; ddp(1, 2) = dpz[1];
  ddp(2, 0) = dpx[2]; ddp(2, 1) = dpy[2]; ddp(2, 2) = dpz[2];

  return val;
}

namespace ImlsSurfaceInternal {

void
traverseKdTree(TriangleKdTree const & kdtree, Vector3 const & query, Array<NodeValuePair> & leaves,
               Array<NodeValuePair> & internal_nodes)
{
  internal_nodes.clear();
  leaves.clear();

  if (kdtree.isEmpty())
    return;

  Array<TriangleKdTree::Node const *> to_examine;
  to_examine.push_back(kdtree.getRoot());

  for (size_t i = 0; i < to_examine.size(); ++i)
  {
    TriangleKdTree::Node const * current = to_examine[i];
    if (current->isLeaf())
      leaves.push_back(NodeValuePair(current, current->getBounds().squaredDistance(query)));
    else
    {
      if (current->getLowChild()->getBounds().contains(query))
        to_examine.push_back(current->getLowChild());
      else
        internal_nodes.push_back(NodeValuePair(current, current->getLowChild()->getBounds().squaredDistance(query)));

      if (current->getHighChild()->getBounds().contains(query))
        to_examine.push_back(current->getHighChild());
      else
        internal_nodes.push_back(NodeValuePair(current, current->getHighChild()->getBounds().squaredDistance(query)));
    }
  }
}

double
triSquaredDistance(TriangleKdTree const & kdtree, TriangleKdTree::Node const & node, Vector3 const & query)
{
  IndexedTriangle const * tris = kdtree.getElements();
  double min_dist2 = std::numeric_limits<double>::infinity();
  for (TriangleKdTree::Node::ElementIndexConstIterator ii = node.elementIndicesBegin(); ii != node.elementIndicesEnd(); ++ii)
  {
    double dist2 = tris[*ii].squaredDistance(query);
    if (dist2 < min_dist2)
      min_dist2 = dist2;
  }

  return min_dist2;
}

TriangleKdTree::Node const *
closestLeaf(TriangleKdTree const & kdtree, Vector3 const & query, Array<NodeValuePair> & explored)
{
  // Find the leaves containing the query point, if any
  Array<NodeValuePair> leaf_queue, internal_queue;
  traverseKdTree(kdtree, query, leaf_queue, internal_queue);

  explored.clear();

  // A comparator to heapify the queues in increasing order of distance from the query
  std::greater<NodeValuePair> cmp;

  double closest_dist2 = std::numeric_limits<double>::infinity();
  TriangleKdTree::Node const * closest_node = nullptr;

  // Find the closest leaf (at triangle level) in the leaf queue
  if (!leaf_queue.empty())
  {
    // Create a minheap on the leaf queue
    std::make_heap(leaf_queue.begin(), leaf_queue.end(), cmp);

    size_t num_remaining = leaf_queue.size();
    while (num_remaining > 0 && leaf_queue.front().value < closest_dist2)
    {
      NodeValuePair const & first = leaf_queue.front();
      double new_dist2 = triSquaredDistance(kdtree, *first.node, query);
      if (new_dist2 < closest_dist2)
      {
        closest_node = first.node;
        closest_dist2 = new_dist2;
      }

      explored.push_back(NodeValuePair(first.node, new_dist2));

      std::pop_heap(leaf_queue.begin(), leaf_queue.end(), cmp);
      num_remaining--;
    }

    // Add remaining leaves
    explored.insert(explored.end(), leaf_queue.begin(), leaf_queue.begin() + num_remaining);
  }

  // Now look for closer nodes by searching the unexplored internal nodes
  std::make_heap(internal_queue.begin(), internal_queue.end(), cmp);
  while (!internal_queue.empty() && internal_queue.front().value < closest_dist2)
  {
    NodeValuePair first = internal_queue.front();
    std::pop_heap(internal_queue.begin(), internal_queue.end(), cmp);
    internal_queue.pop_back();

    if (first.node->isLeaf())  // queued node is a leaf, check if it is the closest by examining triangles
    {
      double new_dist2 = triSquaredDistance(kdtree, *first.node, query);
      if (new_dist2 < closest_dist2)
      {
        closest_node = first.node;
        closest_dist2 = new_dist2;
      }

      explored.push_back(NodeValuePair(first.node, new_dist2));
    }
    else  // queued node is internal, check if its children could be closer than the current minimum
    {
      double lo_dist2 = first.node->getLowChild()->getBounds().squaredDistance(query);
      if (lo_dist2 < closest_dist2)
      {
        internal_queue.push_back(NodeValuePair(first.node->getLowChild(), lo_dist2));
        std::push_heap(internal_queue.begin(), internal_queue.end(), cmp);
      }
      else
        explored.push_back(NodeValuePair(first.node->getLowChild(), lo_dist2));

      double hi_dist2 = first.node->getHighChild()->getBounds().squaredDistance(query);
      if (hi_dist2 < closest_dist2)
      {
        internal_queue.push_back(NodeValuePair(first.node->getHighChild(), hi_dist2));
        std::push_heap(internal_queue.begin(), internal_queue.end(), cmp);
      }
      else
        explored.push_back(NodeValuePair(first.node->getHighChild(), hi_dist2));
    }
  }

  // Add remaining nodes
  explored.insert(explored.end(), internal_queue.begin(), internal_queue.end());

  return closest_node;
}

void
lineIntegrals(double k1, double k2, double & I1, double & Ix)
{
  // These equations are wrong in the paper, they should be (for a = 1):
  // Beta = atan( k1/sqrt(k2) ) - atan( (k1+1)/sqrt(k2) )
  // I1 = -Beta * 1/(2*k2^(3/2))  + (k2 - k1*(k1+1)) / (2*k2*(k1^2+k2)*((k1+1)^2+k2))
  // Ix = Beta * k1/(2*k2^(3/2))  + (k1 + 1) / (2*k2*((k1+1)^2+k2))

  static double const K2_EPSILON = 1.0e-200;
  if (k2 < K2_EPSILON)
  {
    I1 = Ix = 0;
    return;
  }

  double sqrt_k2  =  std::sqrt(k2);
  double k1_p1    =  k1 + 1;

  I1 = Ix = (std::atan(k1_p1 / sqrt_k2) - std::atan(k1 / sqrt_k2) ) / (2 * k2 * sqrt_k2);
  Ix *= -k1;

  double denom = 1.0 / (2 * k2 * (k1_p1 * k1_p1 + k2));
  Ix += k1_p1 * denom;
  I1 += (k2 - k1 * k1_p1) * denom / (k1 * k1 + k2);
}

void
lineIntegrals(double k1, double k2, double & I1, double & Ix, double & dI1, double & dIx, double & dIx2)
{
  // Beta = atan( k1/sqrt(k2) ) - atan( (k1+1)/sqrt(k2) )
  // I1 = -Beta * 1/(2*k2^(3/2))  + (k2 - k1*(k1+1)) / (2*k2*(k1^2+k2)*((k1+1)^2+k2))
  // Ix = Beta * k1/(2*k2^(3/2))  + (k1 + 1) / (2*k2*((k1+1)^2+k2))
  // dI1 = 1/8 * ( -Beta*3/k2^(5/2) +
  //               (5*k2^3 - (k1*(k1+1)-3)*k2^2 - k1*(k1+1)*(9*k1*(k1+1)+5)*k2 - 3*k1^3*(k1+1)^3)
  //               / (k2^2*(k1^2+k2)^2*((k1+1)^2+k2)^2) )
  // dIx = 1/8 * ( Beta*3*k1/k2^(5/2) +
  //               ((k1^2+k2)*((3*(k1+1)+1)*k2^2 + (k1+1)*(6*k1^2+3*k1+2)*k2 + 3*k1^2*(k1+1)^3))
  //               /(k2^2*(k1^2+k2)^2*((k1+1)^2+k2)^2)  )
  // dIx2 = 1/8 * ( -Beta*(3*k1^2+k2)/k2^(5/2) -
  //                ((k1^2+k2)^2*(k2^2 + (k1+1)*(4*k1-1)*k2 + 3*k1*(k1+1)^3))
  //                /(k2^2*(k1^2+k2)^2*((k1+1)^2+k2)^2)  )

  static double const K2_EPSILON = 1.0e-200;
  if (k2 < K2_EPSILON)
  {
    I1 = Ix = dI1 = dIx = dIx2 = 0;
    return;
  }

  double sqrt_k2  =  std::sqrt(k2);
  double k1_p1    =  k1 + 1;
  double k1_2     =  k1 * k1;
  double k2_2     =  k2 * k2;
  double k1_p1_2  =  k1_p1 * k1_p1;
  double t1       =  k2 + k1_p1_2;
  double t2       =  k2 + k1_2;
  double t3       =  3 * k1_p1_2 * k1_p1 * k1;

  I1    =  dI1  =  (std::atan(k1_p1 / sqrt_k2) - std::atan(k1 / sqrt_k2)) / (2 * k2 * sqrt_k2);
  dI1  *=  0.75 / k2;
  Ix    =  -k1 * I1;
  dIx   =  -k1 * dI1;
  dIx2  =  0.25 * I1 * (3 * k1_2 + k2) / k2;

  double denom   =  0.5 / (k2 * t1);
  double ddenom  =  0.125 / (k2_2 * t2 * t2 * t1 * t1);

  I1 += (k2 - k1 * k1_p1) * denom / t2;
  Ix += k1_p1 * denom;

  dI1   +=  (5 * k2 * k2_2 - (k1 * k1_p1 - 3) * k2_2 - k1 * k1_p1 * (9 * k1 * k1_p1 + 5) * k2 - k1_2 * t3) * ddenom;
  dIx   +=  t2 * ((3 * k1_p1 + 1) * k2_2 + k1_p1 * (6 * k1_2 + 3 * k1 + 2) * k2 + k1 * t3) * ddenom;
  dIx2  -=  t2 * t2 * (k2_2 + k1_p1 * (4 * k1 - 1) * k2 + t3) * ddenom;
}

double
lineIntegral(Vector3 const & x, Vector3 const & p0, Vector3 const & p1, double v0, double v1, double eps2_)
{
  using namespace ImlsSurfaceInternal;

  Vector3d ab = (p1 - p0).cast<double>();
  Vector3d xa = (p0 - x).cast<double>();
  double denom = 1 / ab.squaredNorm();
  double k1 = ab.dot(xa) * denom;
  double k2 = (xa.squaredNorm() + eps2_) * denom - k1 * k1;
  double I1, Ix;

  lineIntegrals(k1, k2, I1, Ix);

  return (v0 * I1 + (v1 - v0) * Ix) * std::sqrt(denom) * denom;
}

Vector3d
lineIntegralDeriv(Vector3 const & x, Vector3 const & p0, Vector3 const & p1, double v0, double v1, double eps2_)
{
  Vector3d ab = (p1 - p0).cast<double>();
  Vector3d xa = (p0 - x).cast<double>();
  double denom = 1 / ab.squaredNorm();
  double k1 = ab.dot(xa) * denom;
  double k2 = (xa.squaredNorm() + eps2_) * denom - k1 * k1;
  double I1, Ix, dI1, dIx, dIx2;

  lineIntegrals(k1, k2, I1, Ix, dI1, dIx, dIx2);

  denom = 4 * std::sqrt(denom) * denom * denom;
  return (xa * (v0 * dI1 + (v1 - v0) * dIx)
        + ab * (v0 * dIx + (v1 - v0) * dIx2)) * denom;
}

Vector2d
splitTriangleQuadrature(Vector3 const & x, Vector3 const & pm, Vector3 const & p1, Vector3 const & p2, double vm, double v1,
                        double v2, double eps2_)
{
  Vector3 pp(p1), pn(p2);
  double vp(v1), vn(v2);
  if ((pp - x).squaredNorm() > (pn - x).squaredNorm())
  {
    // Swap so that pp is closest to x
    pp = p2; pn = p1;
    vp = v2; vn = v1;
  }

  Vector3d d1 = (pp - pm).cast<double>();
  Vector3d d2 = (pn - pm).cast<double>();
  Vector3d d3 = (pm - x).cast<double>();

  static double const T1_EPSILON = 1.0e-100;
  double t1 = d1.squaredNorm();
  if (std::fabs(t1) < T1_EPSILON) return Vector2d::Zero();  // early exit if t1 is 0

  double t2 = d1.dot(d2);
  double t3 = d1.dot(d3);
  double t4 = d2.squaredNorm();
  double t5 = d2.dot(d3) * 2;
  double t6 = d3.squaredNorm() + eps2_;

  // Compute height (divided by 2 * sqrt(t1)); early exit if triangle flat
  double height = t4 / t1 - t2 * t2 / (t1 * t1);
  if (height <= 0) return Vector2d::Zero();
  height = std::sqrt(height) / 2.0;

  double vt1 = vn - vm;
  double vt2 = vp - vm;

  double alpha = 2.0 / 3.0;
  double sum1 = 0, sum2 = 0;
  double weight = (1.0 / alpha - alpha);
  double I1, Ix, u_1, denom, k1, k2;
  double u = alpha;
  double last_li1 = 0, last_li2 = 0;

  // Integrate using the trapezoid rule with non-uniform sampling
  for (; u > 0.01f; u *= alpha)
  {
    sum1 += last_li1;
    sum2 += last_li2;

    u_1 = 1.0 - u;
    denom = 1.0 / (u_1 * u_1 * t1);
    k1 = (t3 + u * t2) * u_1 * denom;
    k2 = (t6 + u * t5 + u * u * t4) * denom - k1 * k1;

    lineIntegrals(k1, k2, I1, Ix);

    // THEA_CONSOLE << "k1 = " << k1 << ", k2 = " << k2 << ", I1 = " << I1 << ", Ix = " << Ix
    //              << ", vm = " << vm << ", vt1 = " << vt1 << ", u_1 = " << u_1 << ", vt2 = " << vt2;

    denom *= (u / u_1);
    last_li1 = ((vm + u * vt1) * I1 + u_1 * vt2 * Ix) * denom;
    last_li2 = I1 * denom;
  }

  sum1 *= weight;
  sum1 += last_li1 / alpha;
  sum2 *= weight;
  sum2 += last_li2 / alpha;

  // Add the last trapezoid covering the remaining area
  denom = 1.0 / t1;
  k1 = t3 * denom;
  k2 = t6 * denom - k1 * k1;

  lineIntegrals(k1, k2, I1, Ix);

  denom *= u / alpha;
  sum1 += (vm * I1 + vt2 * Ix) * denom;
  sum2 += I1 * denom;

  sum1 *= height;
  sum2 *= height;

  Vector2d result;
  result[0] = sum1;
  result[1] = sum2;
  return result;
}

IntegralData
splitTriangleQuadratureWithDeriv(Vector3 const & x, Vector3 const & pm, Vector3 const & p1, Vector3 const & p2, double vm,
                                 double v1, double v2, double eps2_)
{
  Vector3 pp(p1), pn(p2);
  double vp(v1), vn(v2);
  if ((pp - x).squaredNorm() > (pn - x).squaredNorm())
  {
    // Swap so that pp is closest to x
    pp = p2; pn = p1;
    vp = v2; vn = v1;
  }

  Vector3d d1 = (pp - pm).cast<double>();
  Vector3d d2 = (pn - pm).cast<double>();
  Vector3d d3 = (pm - x).cast<double>();

  static double const T1_EPSILON = 1.0e-100;
  double t1 = d1.squaredNorm();
  if (std::fabs(t1) < T1_EPSILON) return IntegralData();  // early exit if t1 is 0

  double t2 = d1.dot(d2);
  double t3 = d1.dot(d3);
  double t4 = d2.squaredNorm();
  double t5 = d2.dot(d3) * 2;
  double t6 = d3.squaredNorm() + eps2_;

  // Compute height (divided by 2 * sqrt(t1)); early exit if triangle flat
  double height = t4 / t1 - t2 * t2 / (t1 * t1);
  if (height <= 0.0) return IntegralData();
  height = std::sqrt(height) / 2;

  double vt1 = vn - vm;
  double vt2 = vp - vm;

  double alpha = 2.0 / 3.0;
  IntegralData i_data, last_i_data;
  double weight = (1.0 / alpha - alpha);
  double I1, Ix, dI1, dIx, dIx2, u_1, denom, k1, k2;
  double u = alpha;

  // Integrate using the trapezoid rule with non-uniform sampling
  double const lower_bound = 0.01;  // ((t6<t4)?(t6/t4):1.0) * 0.01;
  for (; u > lower_bound; u *= alpha)
  {
    i_data += last_i_data;

    u_1 = 1.0 - u;
    denom = 1.0 / (u_1 * u_1 * t1);
    k1 = (t3 + u * t2) * u_1 * denom;
    k2 = (t6 + u * t5 + u * u * t4) * denom - k1 * k1;

    lineIntegrals(k1, k2, I1, Ix, dI1, dIx, dIx2);

    double phi_c = vm + u * vt1;
    double phi_x = u_1 * vt2;
    last_i_data.I      =  I1;
    last_i_data.I_phi  =  (phi_c * I1 + phi_x * Ix);

    Vector3d d_c(d3 + u * d2), d_x(d1 * u_1);
    last_i_data.dI      =  (d_c * dI1 + d_x * dIx) * denom;
    last_i_data.dI_phi  =  (d_c * (phi_c * dI1 + phi_x * dIx) + d_x * (phi_c * dIx + phi_x * dIx2) ) * denom;

    denom *= (u / u_1);
    last_i_data *= denom;
  }

  i_data *= weight;
  last_i_data *= 1.0 / alpha;
  i_data += last_i_data;

  // Add the last trapezoid covering the remaining area
  denom = 1.0 / t1;
  k1 = t3 * denom;
  k2 = t6 * denom - k1 * k1;
  lineIntegrals(k1, k2, I1, Ix, dI1, dIx, dIx2);
  denom *= u / alpha;

  i_data.I      +=  I1 * denom;
  i_data.I_phi  +=  (vm * I1 + vt2 * Ix) * denom;

  denom /= t1;
  i_data.dI      +=  (d3 * dI1 + d1 * dIx) * denom;
  i_data.dI_phi  +=  (d3 * (vm * dI1 + vt2 * dIx) + d1 * (vm * dIx + vt2 * dIx2)) * denom;

  i_data         *=  height;
  i_data.dI      *=  4;
  i_data.dI_phi  *=  4;

  return i_data;
}

template <typename T, typename F>
T
triangleQuadrature(F quad_func, Vector3 const & x, Vector3 const & p0, Vector3 const & p1, Vector3 const & p2,
                   Vector3 const & n, double v0, double v1, double v2, double eps2_)
{
  double sqdist, u, v;
  char flag = triangleClosestPoint(x, p0, p1, p2, n, sqdist, u, v);
  switch (flag)
  {
    case 1: return quad_func(x, p0, p1, p2, v0, v1, v2, eps2_);
    case 2: return quad_func(x, p1, p2, p0, v1, v2, v0, eps2_);
    case 4: return quad_func(x, p2, p0, p1, v2, v0, v1, eps2_);

    case 3:
    {
      double t = 1 - u;
      Vector3 pi = static_cast<Real>(t) * p0 + static_cast<Real>(u) * p1;
      double vi = t * v0 + u * v1;
      return quad_func(x, pi, p2, p0, vi, v2, v0, eps2_)
           + quad_func(x, pi, p1, p2, vi, v1, v2, eps2_);
    }

    case 6:
    {
      Vector3 pi = static_cast<Real>(u) * p1 + static_cast<Real>(v) * p2;
      double vi = u * v1 + v * v2;
      return quad_func(x, pi, p0, p1, vi, v0, v1, eps2_)
           + quad_func(x, pi, p2, p0, vi, v2, v0, eps2_);
    }

    case 5:
    {
      double t = 1 - v;
      Vector3 pi = static_cast<Real>(t) * p0 + static_cast<Real>(v) * p2;
      double vi = t * v0 + v * v2;
      return quad_func(x, pi, p0, p1, vi, v0, v1, eps2_)
           + quad_func(x, pi, p1, p2, vi, v1, v2, eps2_);
    }

    case 7:
    {
      double t = 1 - u - v;
      Vector3 pi = static_cast<Real>(t) * p0 + static_cast<Real>(u) * p1 + static_cast<Real>(v) * p2;
      double vi = t * v0 + u * v1 + v * v2;
      return quad_func(x, pi, p0, p1, vi, v0, v1, eps2_)
           + quad_func(x, pi, p1, p2, vi, v1, v2, eps2_)
           + quad_func(x, pi, p2, p0, vi, v2, v0, eps2_);
    }

    default: debugAssertM(false, "ImlsSurface: Invalid flag");  // should never be reached
  }

  return T();  // dummy, to avoid possible compiler warnings
}

char
triangleClosestPoint(Vector3 const & p, Vector3 const & a, Vector3 const & b, Vector3 const & c, Vector3 const & n,
                     double & sqdist, double & u, double & v)
{
  double denom = 1.0 / n.squaredNorm();

  Vector3d ap = (p - a).cast<double>();
  Vector3d bp = (p - b).cast<double>();
  Vector3d cp = (p - c).cast<double>();

  Vector3d t = n.cross(p - a).cast<double>();
  v = bp.dot(t) * denom;
  u = -cp.dot(t) * denom;

  Vector3d ab = (b - a).cast<double>();
  Vector3d bc = (c - b).cast<double>();
  Vector3d ca = (a - c).cast<double>();

  double const eps = std::numeric_limits<double>::epsilon();
  char state = 0;

  if (u <= eps)
  {
    double p_v = v - u * ab.dot(ca) / ca.dot(ca);
    if (p_v <= eps)
      state = 1;
    else if (p_v >= 1)
      state = 4;
    else
    {
      u = 0; v = p_v;
      sqdist = ((1 - v) * ap + v * cp).squaredNorm();
      return 5;
    }
  }

  if (v <= eps)
  {
    double p_u = u - v * ca.dot(ab) / ab.dot(ab);
    if (p_u <= eps)
      state = 1;
    else if (p_u >= 1)
      state = 2;
    else
    {
      u = p_u; v = 0;
      sqdist = ((1 - u) * ap + u * bp).squaredNorm();
      return 3;
    }
  }

  double uv;
  if ((uv = 1 - u - v) <= eps)
  {
    double s = -ca.dot(bc) / bc.dot(bc);
    double p_u = u + uv * s;
    double p_v = v + uv * (1 - s);
    if (p_v <= eps)
      state = 2;
    else if (p_u <= eps)
      state = 4;
    else
    {
      u = p_u; v = p_v;
      sqdist = (u * bp + v * cp).squaredNorm();
      return 6;
    }
  }

  switch (state)
  {
    case 1:
      u = 0; v = 0;
      sqdist = ap.squaredNorm();
      return 1;

    case 2:
      u = 1; v = 0;
      sqdist = bp.squaredNorm();
      return 2;

    case 4:
      u = 0; v = 1;
      sqdist = cp.squaredNorm();
      return 4;

    default:
      sqdist = ap.dot(n.cast<double>());
      sqdist = sqdist * sqdist * denom;
      return 7;
  }

  return 0;  // should never be reached
}

} // namespace ImlsSurfaceInternal

// //=============================================================================
// // External functions
//
// //: find the zero crossing point by bisection between positive point \a pp and negative point \a pn
// //  Stops searching when $||pp-pn|| < xeps$ or $|f(pm)| < feps$
// Vector3d bisect(ImlsSurface const & f,
//                      Vector3d pp,
//                      Vector3d pn,
//                      double feps, double xeps)
// {
//   assert(f(pp) > 0.0);
//   assert(f(pn) < 0.0);
//   Vector3d pm = centre(pp, pn);
//   const unsigned num_itr =
//       static_cast<unsigned>(std::ceil(std::log((pp-pn).norm()
//                                               / xeps)
//                                      / 0.301029995663981)); // log_2
//   Vector3d dp;
//   double val = f.deriv(pm, dp);
//   val /= dp.norm();
//   for (unsigned int i = 0; i<num_itr; ++i) {
//     if (std::abs(val) < feps)
//       return pm;
//     else if (val > 0.0)
//       pp = pm;
//     else
//       pn = pm;
//     pm = centre(pp, pn);
//     val = f.deriv(pm, dp);
//     val /= dp.norm();
//   }
//   return pm;
// }
//
//
// //: Move the point \a p along the gradient direction until reaching a zero crossing of \a f (within \a eps).
// //  Return true if successful
// bool snap_to_surface(ImlsSurface const & f,
//                      Vector3d & p,
//                      double step, double eps)
// {
//   Vector3d p1(p);
//   Vector3d dp;
//   double val1 = f.deriv(p1, dp);
//   double dl = dp.norm();
//   val1 /= dl;
//   if (std::abs(val1) < eps)
//     return true;
//
//   Vector3d p2 = p1 - (step*val1/dl)*dp;
//   dl = dp.norm();
//   double val2 = f.deriv(p2, dp);
//   val2 /= dl;
//   unsigned int i = 0;
//   for (; i<100 & & val1*val2 > 0.0; ++i) {
//     p1 = p2;
//     val1 = val2;
//     p2 -= (step*val2/dl)*dp;
//     val2 = f.deriv(p2, dp);
//     dl = dp.norm();
//     val2 /= dl;
//     if (std::abs(val2) < eps) {
//       p = p2;
//       return true;
//     }
//   }
//   if (i > = 100)
//     return false;
//
//   if (val1 > 0.0)
//     p = bisect(f, p1, p2, eps);
//   else
//     p = bisect(f, p2, p1, eps);
//
//
//   return true;
// }
//
//
// namespace{
// double func(Vector3d const & n,
//             double v,
//             Vector3d const & dp)
// {
//   v *= v;
//   v /= dp.squaredNorm();
//   double tmp = n.dot(dp)/dp.norm() - 1.0;
//   return v + tmp*tmp;
// }
//
// Vector3d dfunc(Vector3d const & n,
//                     double v,
//                     Vector3d const & dp,
//                     vnl_double_3x3 const & ddp)
// {
//   vnl_double_3 nn(n.x(), n.y(), n.z());
//   vnl_double_3 ndp(dp.x(), dp.y(), dp.z());
//   double sqr_len = dp.squaredNorm();
//   double len = std::sqrt(sqr_len);
//   double n_dot_dp = nn.dot(ndp);
//   vnl_double_3 df = (2*v/sqr_len)*ndp;
//   df += ddp.transpose() * ( ((-2*v*v/(sqr_len*sqr_len))*ndp) +
//                             (2/len*(n_dot_dp/len - 1)*(nn - (n_dot_dp/sqr_len)*ndp)) );
//   return Vector3d(df[0], df[1], df[2]);
// }
// // end of namespace
// }
//
// //: Move the point \a p to minimize $(f^2 + (n*f' - 1)^2)/f'*f'$ a zero crossing of \a f (within \a eps).
// //  Return true if successful
// bool snap_to_surfaceWith_normal(ImlsSurface const & f,
//                                 Vector3d & p,
//                                 Vector3d n,
//                                 double step, double eps)
// {
//   Vector3d p1(p);
//   normalize(n);
//   Vector3d dp;
//   vnl_double_3x3 ddp;
//   double val = f.deriv2(p1, dp, ddp);
//   double f1 = func(n, val, dp);
//   if (f1 < eps)
//     return true;
//
//   Vector3d df = dfunc(n, val, dp, ddp);
//   double dl = df.squaredNorm();
//   unsigned int i;
//   for (i = 0; i<1000; ++i) {
//     Vector3d p2 = p1 - (step*f1/dl)*df;
//     val = f.deriv2(p2, dp, ddp);
//     double f2 = func(n, val, dp);
//     //std::cout << i<<" f: "<<f2<<" step: "<< step<<std::endl;
//     if ( f2 > f1) {
//       step /= 2;
//       continue;
//     }
//     if (f2 < eps || step < eps) {
//       p = p2;
//       return true;
//     }
//     f1 = f2;
//     p1 = p2;
//     df = dfunc(n, val, dp, ddp);
//     dl = df.squaredNorm();
//     step *= 2;
//   }
//   if (i > = 100)
//     return true;
//
//   p = p1;
//   return true;
// }
//
//
// //: Move the point \a p along direction \a dir until reaching a zero crossing of \a f (within \a eps).
// //  Return true if successful
// bool snap_to_surface(ImlsSurface const & f,
//                      Vector3d dir,
//                      Vector3d & p,
//                      double step, double eps)
// {
//   Vector3d p1(p);
//   normalize(dir);
//   Vector3d dp;
//   double val1 = f.deriv(p1, dp);
//   dp = dir * dp.dot(dir);
//   double dl = dp.norm();
//   val1 /= dl;
//   if (std::abs(val1) < eps)
//     return true;
//
//   Vector3d p2 = p1 - (step*val1/dl)*dp;
//   dl = dp.norm();
//   double val2 = f.deriv(p2, dp);
//   val2 /= dl;
//   unsigned int i = 0;
//   for (; i<100 & & val1*val2 > 0.0; ++i) {
//     p1 = p2;
//     val1 = val2;
//     p2 -= (step*val2/dl)*dp;
//     val2 = f.deriv(p2, dp);
//     dp = dir * dp.dot(dir);
//     dl = dp.norm();
//     val2 /= dl;
//     if (std::abs(val2) < eps) {
//       p = p2;
//       return true;
//     }
//   }
//   if (i > = 100)
//     return false;
//
//   if (val1 > 0.0)
//     p = bisect(f, p1, p2, eps);
//   else
//     p = bisect(f, p2, p1, eps);
//
//
//   return true;
// }

} // namespace Algorithms
} // namespace Thea
