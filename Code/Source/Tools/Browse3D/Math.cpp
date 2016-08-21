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

#include "Math.hpp"
#include "../../Array.hpp"
#include "../../Ray3.hpp"
#include "../../Triangle3.hpp"
#include "../../Algorithms/SVD.hpp"

namespace Browse3D {

Matrix3
basisMatrix(Vector3 const & u, Vector3 const & v, Vector3 const & w)
{
  return Matrix3(u.x(), v.x(), w.x(),
                 u.y(), v.y(), w.y(),
                 u.z(), v.z(), w.z());
}

Matrix3
orthonormalBasisMatrix(Vector3 const & u, Vector3 const & v, Vector3 const & w)
{
  Vector3 w2 = w.fastUnit();
  Vector3 u2 = v.cross(w2).fastUnit();
  Vector3 v2 = w2.cross(u2);

  return basisMatrix(u2, v2, w2);
}

Matrix3
orthonormalBasisMatrix(Matrix3 const & m)
{
  return orthonormalBasisMatrix(m.getColumn(0), m.getColumn(1), m.getColumn(2));
}

Real
estimateScalingFactor(Matrix3 const & m)
{
  return (m.getColumn(0).fastLength() + m.getColumn(1).fastLength() + m.getColumn(2).fastLength()) / 3.0f;
}

Vector3
reflectPoint(Vector3 const & p, Plane3 const & mirror_plane)
{
  // Assume Plane3::normal() is unit
  return p - 2 * (p - mirror_plane.getPoint()).dot(mirror_plane.getNormal()) * mirror_plane.getNormal();
}

Vector3
reflectVector(Vector3 const & v, Plane3 const & mirror_plane)
{
  // Assume Plane3::normal() is unit
  return v - 2 * v.dot(mirror_plane.getNormal()) * mirror_plane.getNormal();
}

bool
lineIntersectsTriangle(Line3 const & line, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02, Real * time)
{
  // The code is taken from the ray-triangle intersection test in Dave Eberly's Wild Magic library, v5.3, released under the
  // Boost license: http://www.boost.org/LICENSE_1_0.txt .

  static Real const EPS = 1e-30f;

  Vector3 diff = line.getPoint() - v0;
  Vector3 normal = edge01.cross(edge02);

  // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction, E1 = edge01, E2 = edge02, N = Cross(E1,E2)) by
  //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
  //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
  //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)

  Real DdN = line.getDirection().dot(normal);
  int sign;
  if (DdN > EPS)
    sign = 1;
  else if (DdN < -EPS)
  {
    sign = -1;
    DdN = -DdN;
  }
  else
  {
    // Line and triangle are parallel, call it a "no intersection" even if the line does intersect
    return false;
  }

  Real DdQxE2 = sign * line.getDirection().dot(diff.cross(edge02));
  if (DdQxE2 >= 0)
  {
    Real DdE1xQ = sign * line.getDirection().dot(edge01.cross(diff));
    if (DdE1xQ >= 0)
    {
      if (DdQxE2 + DdE1xQ <= DdN)
      {
        // Line intersects triangle
        Real QdN = -sign * diff.dot(normal);
        if (time) *time = QdN / DdN;
        return true;
      }
      // else: b1 + b2 > 1, no intersection
    }
    // else: b2 < 0, no intersection
  }
  // else: b1 < 0, no intersection

  return false;
}

#define CLOSEST_PT_SEGMENT_SEGMENT(VecType) \
Real \
closestPtSegmentSegment(VecType const & p1, VecType const & q1, VecType const & p2, VecType const & q2, Real & s, Real & t, \
                        VecType & c1, VecType & c2) \
{ \
  static Real const EPSILON = 1e-20f; \
 \
  VecType d1 = q1 - p1; /* Direction vector of segment S1 */ \
  VecType d2 = q2 - p2; /* Direction vector of segment S2 */ \
  VecType r = p1 - p2; \
  Real a = d1.squaredLength(); /* Squared length of segment S1, always nonnegative */ \
  Real e = d2.squaredLength(); /* Squared length of segment S2, always nonnegative */ \
  Real f = d2.dot(r); \
 \
  /* Check if either or both segments degenerate into points */ \
  if (a <= EPSILON && e <= EPSILON) \
  { \
    /* Both segments degenerate into points */ \
    s = t = 0; \
    c1 = p1; \
    c2 = p2; \
    return (c1 - c2).dot(c1 - c2); \
  } \
  if (a <= EPSILON) \
  { \
    /* First segment degenerates into a point */ \
    s = 0; \
    t = f / e; /* s = 0 => t = (b*s + f) / e = f / e */ \
    t = Math::clamp(t, (Real)0, (Real)1); \
  } \
  else \
  { \
    Real c = d1.dot(r); \
    if (e <= EPSILON) \
    { \
        /* Second segment degenerates into a point */ \
        t = 0; \
        s = Math::clamp(-c / a, (Real)0, (Real)1); /* t = 0 => s = (b*t - c) / a = -c / a */ \
    } \
    else \
    { \
      /* The general nondegenerate case starts here */ \
      Real b = d1.dot(d2); \
      Real denom = a * e - b * b; /* Always nonnegative */ \
 \
      /* If segments not parallel, compute closest point on L1 to L2, and */ \
      /* clamp to segment S1. Else pick arbitrary s (here 0) */ \
      if (denom != 0) \
        s = Math::clamp((b * f - c * e) / denom, (Real)0, (Real)1); \
      else \
        s = 0; \
 \
      /* Compute point on L2 closest to S1(s) using */ \
      /* t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e */ \
      t = (b * s + f) / e; \
 \
      /* If t in [0,1] done. Else clamp t, recompute s for the new value */ \
      /* of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a */ \
      /* and clamp s to [0, 1] */ \
      if (t < 0) \
      { \
        t = 0; \
        s = Math::clamp(-c / a, (Real)0, (Real)1); \
      } \
      else if (t > 1) \
      { \
        t = 1; \
        s = Math::clamp((b - c) / a, (Real)0, (Real)1); \
      } \
    } \
  } \
 \
  c1 = p1 + s * d1; \
  c2 = p2 + t * d2; \
  return (c1 - c2).squaredLength(); \
}

CLOSEST_PT_SEGMENT_SEGMENT(Vector2)
CLOSEST_PT_SEGMENT_SEGMENT(Vector3)

#define CLOSEST_PT_SEGMENT_LINE(VecType) \
Real \
closestPtSegmentLine(VecType const & p1, VecType const & q1, VecType const & p2, VecType const & q2, Real & s, Real & t, \
                     VecType & c1, VecType & c2) \
{ \
  static Real const EPSILON = 1e-20f; \
 \
  VecType d1 = q1 - p1; /* Direction vector of segment S1 */ \
  VecType d2 = q2 - p2; /* Direction vector of line L2 */ \
  VecType r = p1 - p2; \
  Real a = d1.squaredLength(); /* Squared length of segment S1, always nonnegative */ \
  Real e = d2.squaredLength(); /* Squared length of segment supporting L2, always nonnegative and assumed to be non-zero */ \
  Real f = d2.dot(r); \
 \
  if (a <= EPSILON) \
  { \
    /* Segment degenerates into a point */ \
    s = 0; \
    t = f / e; /* s = 0 => t = (b*s + f) / e = f / e */ \
  } \
  else \
  { \
    /* The general nondegenerate case starts here */ \
    Real c = d1.dot(r); \
    Real b = d1.dot(d2); \
    Real denom = a * e - b * b; /* Always nonnegative */ \
\
    /* If segments not parallel, compute closest point on L1 to L2, and */ \
    /* clamp to segment S1. Else pick arbitrary s (here 0) */ \
    if (denom != 0) \
      s = Math::clamp((b * f - c * e) / denom, (Real)0, (Real)1); \
    else \
      s = 0; \
\
    /* Compute point on L2 closest to S1(s) using */ \
    /* t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e */ \
    t = (b * s + f) / e; \
  } \
 \
  c1 = p1 + s * d1; \
  c2 = p2 + t * d2; \
  return (c1 - c2).squaredLength(); \
}

CLOSEST_PT_SEGMENT_LINE(Vector2)
CLOSEST_PT_SEGMENT_LINE(Vector3)

Real
closestPtLineTriangle(Line3 const & line, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02, Real & s,
                      Vector3 & c1, Vector3 & c2)
{
  if (lineIntersectsTriangle(line, v0, edge01, edge02, &s))
  {
    c1 = c2 = line.getPoint() + s * line.getDirection();
    return 0;
  }

  // Either (1) the line is not parallel to the triangle and the point of
  // intersection of the line and the plane of the triangle is outside the
  // triangle or (2) the line and triangle are parallel.  Regardless, the
  // closest point on the triangle is on an edge of the triangle.  Compare
  // the line to all three edges of the triangle.
  Real min_sqdist = -1;
  Vector3 v[3] = { v0, v0 + edge01, v0 + edge02 };
  Vector3 tmp_c1, tmp_c2;
  float tmp_s, t;
  for (int i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
  {
    Real sqdist = closestPtSegmentLine(v[i0], v[i1], line.getPoint(), line.getPoint() + line.getDirection(),
                                       t, tmp_s, tmp_c2, tmp_c1);
    if (min_sqdist < 0 || sqdist < min_sqdist)
    {
      min_sqdist = sqdist;
      s = tmp_s;
      c1 = tmp_c1;
      c2 = tmp_c2;
    }
  }

  return min_sqdist;
}

Real
closestPtRayTriangle(Ray3 const & ray, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02, Real & s,
                     Vector3 & c1, Vector3 & c2)
{
  Real sqdist = closestPtLineTriangle(Line3::fromPointAndDirection(ray.getOrigin(), ray.getDirection()), v0, edge01, edge02,
                                      s, c1, c2);
  if (s >= 0)
    return sqdist;

  // Not the most efficient way but the most convenient
  LocalTriangle3 tri(v0, v0 + edge01, v0 + edge02);
  s = 0;
  c1 = ray.getOrigin();
  c2 = tri.closestPoint(c1);

  return (c1 - c2).squaredLength();
}

RigidTransform3
closestRigidTransform(std::vector<CoordinateFrame3> const & src, std::vector<CoordinateFrame3> const & dst)
{
  alwaysAssertM(src.size() == dst.size(), "Source and destination sets must have same number of frames");

  if (src.empty())
    return RigidTransform3::identity();
  else if (src.size() == 1)
    return RigidTransform3(dst[0]) * RigidTransform3(src[0].inverse());
  else if (src.size() == 2)
  {
    // First map line segment to line segment...
    Vector3 src_mean = 0.5f * (src[0].getTranslation() + src[1].getTranslation());
    Vector3 dst_mean = 0.5f * (dst[0].getTranslation() + dst[1].getTranslation());

    static Real const MIN_SQLEN = 1.0e-8f;

    Vector3 src_axis = src[1].getTranslation() - src[0].getTranslation();
    Vector3 dst_axis = dst[1].getTranslation() - dst[0].getTranslation();

    if (src_axis.squaredLength() > MIN_SQLEN && dst_axis.squaredLength() > MIN_SQLEN)
    {
      Matrix3 rot = Matrix3::rotationArc(src_axis, dst_axis);

      // Now rotate around the mapped segment to align the z axes of the frames...
      Vector3 src_dir = rot * (src[0].getRotation().getColumn(2) + src[1].getRotation().getColumn(2));
      Vector3 dst_dir = dst[0].getRotation().getColumn(2) + dst[1].getRotation().getColumn(2);

      // Transform src_dir and dst_dir to be perpendicular to dst_axis
      dst_axis = dst_axis.fastUnit();
      src_dir = src_dir - src_dir.dot(dst_axis) * dst_axis;
      dst_dir = dst_dir - dst_dir.dot(dst_axis) * dst_axis;

      if (src_dir.squaredLength() > MIN_SQLEN && dst_dir.squaredLength() > MIN_SQLEN)
      {
        src_dir = src_dir.fastUnit();
        dst_dir = dst_dir.fastUnit();

        Vector3 src_perp = dst_axis.cross(src_dir);
        Vector3 dst_perp = dst_axis.cross(dst_dir);

        Matrix3 src_basis = basisMatrix(src_dir, src_perp, dst_axis);
        Matrix3 dst_basis = basisMatrix(dst_dir, dst_perp, dst_axis);
        Matrix3 extra_rot = dst_basis * src_basis.transpose();  // inverse == tranpose for rotation matrices

        rot = extra_rot * rot;
      }

      CoordinateFrame3 cf(RigidTransform3::_fromAffine(AffineTransform3(rot, dst_mean - rot * src_mean)));
      // THEA_CONSOLE.nospace() << "R = " << cf.getRotation().toString() << ", T = " << cf.getTranslation().toString();

      return RigidTransform3(cf);
    }
    else
      return RigidTransform3::translation(dst_mean - src_mean);
  }

  // ICP algorithm of Arun et al. 1987.

  size_t num_frames = src.size();
  Vector3 src_mean = Vector3::zero(), dst_mean = Vector3::zero();
  for (size_t i = 0; i < num_frames; ++i)
  {
    CoordinateFrame3 const & src_frame = src[i];
    CoordinateFrame3 const & dst_frame = dst[i];

    // THEA_CONSOLE.nospace() << "src[" << i << "] = R: " << src_frame.getRotation().toString() << ", T: "
    //                    << src_frame.getTranslation().toString();
    // THEA_CONSOLE.nospace() << "dst[" << i << "] = R: " << dst_frame.getRotation().toString() << ", T: "
    //                    << dst_frame.getTranslation().toString();

    src_mean += src_frame.getTranslation();
    dst_mean += dst_frame.getTranslation();
  }

  src_mean /= num_frames;
  dst_mean /= num_frames;

  Matrix3 corr = Matrix3::zero();
  Vector3 src_pt, dst_pt;
  for (size_t i = 0; i < num_frames; ++i)
  {
    CoordinateFrame3 const & src_frame = src[i];
    CoordinateFrame3 const & dst_frame = dst[i];

    src_pt = src_frame.getTranslation() - src_mean;
    dst_pt = dst_frame.getTranslation() - dst_mean;

    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c)
        corr(r, c) += src_pt[r] * dst_pt[c];
  }

  Matrix3 U, V;
  TheaArray<Real> diag;
  Algorithms::SVD::compute(corr, U, diag, V);
  Matrix3 U_T = U.transpose();
  Matrix3 rotation = V * U_T;
  if (rotation.determinant() < 0)
  {
    // FIXME: One of the columns (which?) of V should be negated. See Eggert et al. '97, "Estimating 3D rigid transformations: a
    // comparison of four major algorithms".
    //
    // For now, we'll do the dumb but safe thing and test each option for the minimum error.

    THEA_WARNING << "Estimated rigid transform involves reflection, fixing by negating a column of V";

    Real min_sqerr = -1;
    for (int i = 0; i < 3; ++i)
    {
      // Swap i'th column of V
      Matrix3 V_i = V;
      V_i(0, i) = -V_i(0, i);
      V_i(1, i) = -V_i(1, i);
      V_i(2, i) = -V_i(2, i);

      Matrix3 candidate_rot = V_i * U_T;
      Real sqerr = 0;
      for (size_t j = 0; j < num_frames; ++j)
      {
        CoordinateFrame3 const & src_frame = src[j];
        CoordinateFrame3 const & dst_frame = dst[j];

        src_pt = src_frame.getTranslation() - src_mean;
        dst_pt = dst_frame.getTranslation() - dst_mean;
        sqerr += (candidate_rot * src_pt - dst_pt).squaredLength();
      }

      if (min_sqerr < 0 || sqerr < min_sqerr)
      {
        rotation = candidate_rot;
        min_sqerr = sqerr;
      }
    }
  }

  Vector3 translation = dst_mean - rotation * src_mean;

  CoordinateFrame3 cf(RigidTransform3::_fromAffine(AffineTransform3(rotation, translation)));
  // THEA_CONSOLE.nospace() << "R = " << rotation.toString() << ", T = " << translation.toString();

  return RigidTransform3(cf);
}

int
kdtreeDepth(long num_elems, int max_elems_in_leaf)
{
  alwaysAssertM(num_elems >= 0, "Can't compute kd-tree depth for negative number of elements");
  alwaysAssertM(max_elems_in_leaf > 0, "Can't compute kd-tree depth for non-positive number of elements at leaf");

  int max_depth = (num_elems > 0 ? (int)std::ceil(Math::fastLog2(num_elems / (double)max_elems_in_leaf)) : 0);
  return max_depth;
}

void
getBarycentricCoordinates2(Real qx, Real qy, Real x0, Real y0, Real x1, Real y1, Real x2, Real y2,
                           Real & bc0, Real & bc1, Real & bc2)
{
  Real a = y1 - y2;
  Real b = x0 - x2;
  Real c = y2 - y0;
  Real d = x2 - x1;

  Real det = a * b - c * d;
  if (std::abs(det) < 1.0e-30f)
    bc0 = bc1 = bc2 = 1.0f / 3;
  else
  {
    Real qx_rel_3 = qx - x2;
    Real qy_rel_3 = qy - y2;

    bc0 = (a * qx_rel_3 + d * qy_rel_3) / det;
    bc1 = (c * qx_rel_3 + b * qy_rel_3) / det;
    bc2 = 1 - bc0 - bc1;
  }
}

void
getBarycentricCoordinates2(Vector2 const & q, Vector2 const & v0, Vector2 const & v1, Vector2 const & v2,
                           Real & bc0, Real & bc1, Real & bc2)
{
  getBarycentricCoordinates2(q.x(), q.y(), v0.x(), v0.y(), v1.x(), v1.y(), v2.x(), v2.y(), bc0, bc1, bc2);
}

void
getBarycentricCoordinates3(Vector3 const & q, Vector3 const & v0, Vector3 const & v1, Vector3 const & v2,
                           Real & bc0, Real & bc1, Real & bc2)
{
  switch ((v1 - v0).cross(v2 - v0).maxAbsAxis())
  {
    case 0:           getBarycentricCoordinates2(q.yz(), v0.yz(), v1.yz(), v2.yz(), bc0, bc1, bc2); break;
    case 1:           getBarycentricCoordinates2(q.xz(), v0.xz(), v1.xz(), v2.xz(), bc0, bc1, bc2); break;
    default /* 2 */:  getBarycentricCoordinates2(q.xy(), v0.xy(), v1.xy(), v2.xy(), bc0, bc1, bc2);
  }
}

Vector3
triRandomPoint(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2)
{
  return triRandomPointFromEdges(v0, v1 - v0, v2 - v0);
}

Vector3
triRandomPointFromEdges(Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02)
{
  // From G3D::Triangle

  // Choose a random point in the parallelogram

  float s = Random::common().uniform01();
  float t = Random::common().uniform01();

  if (t > 1.0f - s)
  {
    // Outside the triangle; reflect about the diagonal of the parallelogram
    t = 1.0f - t;
    s = 1.0f - s;
  }

  return s * edge01 + t * edge02 + v0;
}

Vector3
quadRandomPoint(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2, Vector3 const & v3)
{
  Vector3 edge01 = v1 - v0;
  Vector3 edge03 = v3 - v0;
  Vector3 edge21 = v1 - v2;
  Vector3 edge23 = v3 - v2;

  Real tri0_double_area = edge01.cross(edge03).fastLength();
  Real tri1_double_area = edge21.cross(edge23).fastLength();

  float tri = Random::common().uniform(0, tri0_double_area + tri1_double_area);
  if (tri <= tri0_double_area)
    return triRandomPointFromEdges(v0, edge01, edge03);
  else
    return triRandomPointFromEdges(v2, edge21, edge23);
}

Real
raySphereIntersectionTime(Ray3 const & ray, Vector3 const & center, Real radius)
{
  Vector3 pmc = ray.getOrigin() - center;

  double s[3];
  s[2] = ray.getDirection().squaredLength();
  s[1] = 2 * pmc.dot(ray.getDirection());
  s[0] = pmc.squaredLength() - radius * radius;

  double roots[2];
  int num_roots = Math::solveQuadratic(s[0], s[1], s[2], roots);

  double min_root = -1;
  for (int i = 0; i < num_roots; ++i)
    if (roots[i] >= 0 && (min_root < 0 || roots[i] < min_root))
      min_root = roots[i];

#if 0
  if (min_root >= 0)
    THEA_CONSOLE << "min_root =" << min_root;
#endif

  return (Real)min_root;
}

Real
rayTorusIntersectionTime(Ray3 const & ray, Real torus_radius, Real torus_width)
{
  double r2pw2 = torus_radius * torus_radius + torus_width * torus_width;
  double r2mw2 = r2pw2 - 2 * torus_width * torus_width;

  Vector3 p2 = ray.getOrigin() * ray.getOrigin();
  Vector3 pu = ray.getOrigin() * ray.getDirection();
  Vector3 u2 = ray.getDirection() * ray.getDirection();

  double s[5];
  s[4] = u2[0] * (u2[0] + 2 * u2[1]) + u2[1] * (u2[1] + 2 * u2[2]) + u2[2] * (u2[2] + 2 * u2[0]);
  s[3] = 4 * (pu[0] + pu[1] + pu[2]) * (u2[0] + u2[1] + u2[2]);
  s[2] = 2 * (r2mw2 * u2[2] - r2pw2 * (u2[0] + u2[1]))
       + 8 * (pu[0] * pu[1] + pu[1] * pu[2] + pu[2] * pu[0])
       + 6 * (pu[0] * pu[0] + pu[1] * pu[1] + pu[2] * pu[2])
       + 2 * (p2[0] * (u2[1] + u2[2]) + p2[1] * (u2[2] + u2[0]) + p2[2] * (u2[0] + u2[1]));
  s[1] = 4 * (r2mw2 * pu[2] - r2pw2 * (pu[0] + pu[1]) + (p2[0] + p2[1] + p2[2]) * (pu[0] + pu[1] + pu[2]));
  s[0] = 2 * (r2mw2 * p2[2] - r2pw2 * (p2[0] + p2[1]) + p2[0] * p2[1] + p2[1] * p2[2] + p2[2] * p2[0])
       + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2] + r2mw2 * r2mw2;

  double roots[4];
  int num_roots = Math::solveQuartic(s[0], s[1], s[2], s[3], s[4], roots);

  double min_root = -1;
  for (int i = 0; i < num_roots; ++i)
    if (roots[i] >= 0 && (min_root < 0 || roots[i] < min_root))
      min_root = roots[i];

#if 0
  if (min_root >= 0)
    THEA_CONSOLE << "min_root =" << min_root;
#endif

  return (Real)min_root;
}

} // namespace Browse3D
