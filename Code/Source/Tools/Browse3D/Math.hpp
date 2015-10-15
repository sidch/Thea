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

#ifndef __Browse3D_Math_hpp__
#define __Browse3D_Math_hpp__

#include "Common.hpp"
#include "../../AffineTransform3.hpp"
#include "../../CoordinateFrame3.hpp"
#include "../../Line3.hpp"
#include "../../Math.hpp"
#include "../../Plane3.hpp"
#include "../../Ray3.hpp"
#include "../../RigidTransform3.hpp"

namespace Browse3D {

namespace MathInternal {

// exp(-x) approximations from http://www.xnoiz.co.cc/fast-exp-x/
float fastMinuzExp_1(float x);
float fastMinuzExp_2(float x);
float fastMinuzExp_3(float x);
float fastMinuzExpWider_1(float x);
float fastMinuzExpWider_2(float x);
float fastMinuzExpWider_3(float x);

} // namespace MathInternal

/** Computes a fast approximation to exp(-x). */
inline float
fastMinusExp(float x)
{
  return MathInternal::fastMinuzExpWider_1(x);
}

/** Computes a fast approximation to exp(-x) for 0 <= x <= 1. */
inline float
fastMinusExp01(float x)
{
  if (x > 0.69)
  {
    float y = MathInternal::fastMinuzExp_1(0.5 * x);
    return y * y;
  }
  else
    return MathInternal::fastMinuzExp_1(x);
}

/**
 * Compute the normalized Epanechnikov kernel for a precomputed squared distance and bandwidth. See Comaniciu and Meer, Mean
 * Shift: A Robust Approach towards Feature Space Analysis, PAMI 24(5), 2002.
 */
inline float
kernelEpanechnikovSqDist(float squared_dist, float squared_bandwidth)
{
  // These constants assume dim = 6
  static float const VOL_B6 = 5.1677127800499694f;  // volume of 6-D unit ball
  static float const SCALE = 0.5f * (6 + 2) / VOL_B6;

  float nrm_sqdist = squared_dist / squared_bandwidth;
  return nrm_sqdist < 1 ? SCALE * (1 - nrm_sqdist) : 0;
}

/**
 * Compute the normalized Epanechnikov kernel for a given distance and bandwidth. See Comaniciu and Meer, Mean Shift: A Robust
 * Approach towards Feature Space Analysis, PAMI 24(5), 2002.
 */
inline float
kernelEpanechnikov(float dist, float bandwidth)
{
  return kernelEpanechnikovSqDist(dist * dist, bandwidth * bandwidth);
}

/**
 * Compute the unnormalized Epanechnikov kernel for a precomputed squared distance and bandwidth. See Comaniciu and Meer, Mean
 * Shift: A Robust Approach towards Feature Space Analysis, PAMI 24(5), 2002.
 */
inline float
kernelEpanechnikovSqDistUnscaled(float squared_dist, float squared_bandwidth)
{
  float nrm_sqdist = squared_dist / squared_bandwidth;
  return nrm_sqdist < 1 ? (1 - nrm_sqdist) : 0;
}

/**
 * Compute the unnormalized Epanechnikov kernel for a given distance and bandwidth. See Comaniciu and Meer, Mean Shift: A Robust
 * Approach towards Feature Space Analysis, PAMI 24(5), 2002.
 */
inline float
kernelEpanechnikovUnscaled(float dist, float bandwidth)
{
  return kernelEpanechnikovSqDistUnscaled(dist * dist, bandwidth * bandwidth);
}

/**
 * Compute the unnormalized Gaussian kernel for a precomputed squared distance and bandwidth, using a fast approximation of
 * exponentiation.
 */
inline float
kernelFastGaussianSqDistUnscaled(float squared_dist, float squared_bandwidth)
{
  return fastMinusExp(squared_dist / squared_bandwidth);
}

/** Get a matrix whose columns are three given vectors. The resulting matrix transforms *from* basis space *to* world space. */
Matrix3 basisMatrix(Vector3 const & u, Vector3 const & v, Vector3 const & w);

/**
 * Convert three vectors, assumed to be linearly independent, to an orthonormal basis matrix. The direction of the W/Z axis is
 * always preserved.
 */
Matrix3 orthonormalBasisMatrix(Vector3 const & u, Vector3 const & v, Vector3 const & w);

/**
 * Convert three columns of a matrix, assumed to be linearly independent, to an orthonormal basis matrix. The direction of the
 * W/Z axis (third column) is always preserved.
 */
Matrix3 orthonormalBasisMatrix(Matrix3 const & m);

/** Estimate the amount by which a 3x3 matrix scales distances in 3-space. The estimate is exact for a similarity transform. */
Real estimateScalingFactor(Matrix3 const & m);

/** Reflect a point in a plane */
Vector3 reflectPoint(Vector3 const & p, Plane3 const & mirror_plane);

/** Reflect a direction vector (need not be unit) in a plane */
Vector3 reflectVector(Vector3 const & v, Plane3 const & mirror_plane);

/**
 * Check if a line intersects a triangle and optionally return the hit time. Identical to ray-triangle intersection except the
 * hit time in this case is allowed to be negative.
 */
bool lineIntersectsTriangle(Line3 const & line, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02,
                            Real * time = NULL);

/**
 * Closest pair of points betwen two line segments (p1, q1) and (p2, q2) in 2D. The points are returned in c1 and c2.
 * s = ||c1 - p1|| / ||q1 - p1|| and t = ||c2 - p2|| / ||q2 - p2||.
 *
 * From Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.
 */
Real closestPtSegmentSegment(Vector2 const & p1, Vector2 const & q1, Vector2 const & p2, Vector2 const & q2, Real & s, Real & t,
                             Vector2 & c1, Vector2 & c2);

/**
 * Closest pair of points betwen two line segments (p1, q1) and (p2, q2) in 3D. The points are returned in c1 and c2.
 * s = ||c1 - p1|| / ||q1 - p1|| and t = ||c2 - p2|| / ||q2 - p2||.
 *
 * From Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman, 2005.
 */
Real closestPtSegmentSegment(Vector3 const & p1, Vector3 const & q1, Vector3 const & p2, Vector3 const & q2, Real & s, Real & t,
                             Vector3 & c1, Vector3 & c2);

/**
 * Closest pair of points betwen a line segment (p1, q1) and an unbounded line supported by the segment (p2, q2) in 2D. The
 * points are returned in c1 and c2. s = ||c1 - p1|| / ||q1 - p1|| and t = ||c2 - p2|| / ||q2 - p2||.
 *
 * Modified from segment-segment closest points test, from Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman,
 * 2005.
 */
Real closestPtSegmentLine(Vector2 const & p1, Vector2 const & q1, Vector2 const & p2, Vector2 const & q2, Real & s, Real & t,
                          Vector2 & c1, Vector2 & c2);

/**
 * Closest pair of points betwen a line segment (p1, q1) and an unbounded line supported by the segment (p2, q2) in 3D. The
 * points are returned in c1 and c2. s = ||c1 - p1|| / ||q1 - p1|| and t = ||c2 - p2|| / ||q2 - p2||.
 *
 * Modified from segment-segment closest points test, from Christer Ericson, "Real-Time Collision Detection", Morgan-Kaufman,
 * 2005.
 */
Real closestPtSegmentLine(Vector3 const & p1, Vector3 const & q1, Vector3 const & p2, Vector3 const & q2, Real & s, Real & t,
                          Vector3 & c1, Vector3 & c2);

/**
 * Closest pair of points betwen a line and triangle. The points are returned in c1 and c2. s is the parameter of c1 w.r.t. the
 * line.
 *
 * Follows Dave Eberly's Wild Magic library.
 */
Real closestPtLineTriangle(Line3 const & line, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02, Real & s,
                           Vector3 & c1, Vector3 & c2);

/**
 * Closest pair of points betwen a ray and triangle. The points are returned in c1 and c2. s is the hit time of the ray to c1.
 *
 * Follows Dave Eberly's Wild Magic library.
 */
Real closestPtRayTriangle(Ray3 const & ray, Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02, Real & s,
                          Vector3 & c1, Vector3 & c2);

/** Get the rigid transform that best maps one set of coordinate frames to another. */
RigidTransform3 closestRigidTransform(std::vector<CoordinateFrame3> const & src, std::vector<CoordinateFrame3> const & dst);

/** Get the barycentric coordinates of a point (qx, qy) with respect to a triangle [(x0, y0), (x1, y1), (x2, y2)], in 2D */
void getBarycentricCoordinates2(Real qx, Real qy, Real x0, Real y0, Real x1, Real y1, Real x2, Real y2,
                                Real & bc0, Real & bc1, Real & bc2);

/** Get the barycentric coordinates of a point q with respect to a triangle [v0, v1, v2], in 2D */
void getBarycentricCoordinates2(Vector2 const & q, Vector2 const & v0, Vector2 const & v1, Vector2 const & v2,
                                Real & bc0, Real & bc1, Real & bc2);

/**
 * Get the barycentric coordinates of a point q with respect to a triangle [v0, v1, v2], in 3D. q is assumed to lie in the plane
 * of the triangle.
 */
void getBarycentricCoordinates3(Vector3 const & q, Vector3 const & v0, Vector3 const & v1, Vector3 const & v2,
                                Real & bc0, Real & bc1, Real & bc2);

/** Get a uniformly distributed random point in a 3D triangle specified by its three vertices. */
Vector3 triRandomPoint(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2);

/**
 * Get a uniformly distributed random point in a 3D triangle specified by a vertex and the two directed edges emanating from it.
 */
Vector3 triRandomPointFromEdges(Vector3 const & v0, Vector3 const & edge01, Vector3 const & edge02);

/** Get a uniformly distributed random point in a 3D quadrilateral specified by its four vertices. */
Vector3 quadRandomPoint(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2, Vector3 const & v3);

/** Hermite interpolation, given two points and the associated tangents. */
template <typename T>
T
cspline(T const & a, T const & da, T const & b, T const & db, Real t)
{
  Real t2 = t * t;
  Real t3 = t * t2;
  Real c = 2 * t3 - 3 * t2;
  return (c + 1) * a
       + (t3 - 2 * t2 + t) * da
       - c * b
       + (t3 - t2) * db;
}

/** Hit-time of nearest intersection of a ray with a sphere. */
Real raySphereIntersectionTime(Ray3 const & ray, Vector3 const & center, Real radius);

namespace MathInternal {

// exp(-x) approximations from http://www.xnoiz.co.cc/fast-exp-x/

// approx method, returns exp(-x) when 0<=x<=ln(2) {~0.69}
inline float fastMinuzExp_1(float x) {

        // err <= 3e-3
        return 1
                -x*(0.9664
                -x*(0.3536));

}

inline float fastMinuzExp_2(float x) {
        // err <= 3e-5

        return 1
                -x*(0.9998684
                -x*(0.4982926

                -x*(0.1595332
                -x*(0.0293641))));

}

inline float fastMinuzExp_3(float x) {
        // err <= 3e-10

        return 1
        -x*(0.9999999995
        -x*(0.4999999206

        -x*(0.1666653019
        -x*(0.0416573475
        -x*(0.0083013598

        -x*(0.0013298820
        -x*(0.0001413161)))))));

}

// widen up fastMinuzExp
inline float fastMinuzExpWider_1(float x) {

        bool lessZero = false;
        if (x < 0) {

                lessZero = true;
                x = -x;

        }
        int mult = 0;
        while (x > 0.69*2*2*2*2*2*2) {

                mult+=6;
                x /= 64.0f;
        }

        while (x > 0.69*2*2*2) {

                mult+=3;
                x /= 8.0f;
        }

        while (x > 0.69*2*2) {
                mult+=2;

                x /= 4.0f;
        }
        while (x > 0.69) {

                mult++;
                x /= 2.0f;
        }

        x = fastMinuzExp_1(x);
        while (mult) {

                mult--;
                x = x*x;
        }

        if (lessZero) {
                return 1 / x;

        } else {
                return x;
        }
}

// widen up fastMinuzExp
inline float fastMinuzExpWider_2(float x) {

        bool lessZero = false;
        if (x < 0) {

                lessZero = true;
                x = -x;

        }
        int mult = 0;
        while (x > 0.69*2*2*2*2*2*2) {

                mult+=6;
                x /= 64.0f;
        }

        while (x > 0.69*2*2*2) {

                mult+=3;
                x /= 8.0f;
        }

        while (x > 0.69*2*2) {
                mult+=2;

                x /= 4.0f;
        }
        while (x > 0.69) {

                mult++;
                x /= 2.0f;
        }

        x = fastMinuzExp_2(x);
        while (mult) {

                mult--;
                x = x*x;
        }

        if (lessZero) {
                return 1 / x;

        } else {
                return x;
        }
}

// widen up fastMinuzExp
inline float fastMinuzExpWider_3(float x) {

        bool lessZero = false;
        if (x < 0) {

                lessZero = true;
                x = -x;

        }
        int mult = 0;
        while (x > 0.69*2*2*2*2*2*2) {

                mult+=6;
                x /= 64.0f;
        }

        while (x > 0.69*2*2*2) {

                mult+=3;
                x /= 8.0f;
        }

        while (x > 0.69*2*2) {
                mult+=2;

                x /= 4.0f;
        }
        while (x > 0.69) {

                mult++;
                x /= 2.0f;
        }

        x = fastMinuzExp_3(x);
        while (mult) {

                mult--;
                x = x*x;
        }

        if (lessZero) {
                return 1 / x;

        } else {
                return x;
        }
}

} // namespace MathInternal

} // namespace Browse3D

#endif
