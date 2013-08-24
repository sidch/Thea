//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#include "Quat.hpp"

namespace Thea {

Quat::Quat(Matrix3 const & rot)
{
  static int const plus1mod3[] = {1, 2, 0};

  // Find the index of the largest diagonal component. These ? operations hopefully compile to conditional move instructions
  // instead of branches.
  int i = (rot(1, 1) > rot(0, 0)) ? 1 : 0;
  i = (rot(2, 2) > rot(i, i)) ? 2 : i;

  // Find the indices of the other elements
  int j = plus1mod3[i];
  int k = plus1mod3[j];

  // If we attempted to pre-normalize and trusted the matrix to be perfectly orthonormal, the result would be:
  //
  //   double c = sqrt((rot(i, i) - (rot(j, j) + rot(k, k))) + 1.0)
  //   v[i] = -c * 0.5
  //   v[j] = -(rot(i, j) + rot(j, i)) * 0.5 / c
  //   v[k] = -(rot(i, k) + rot(k, i)) * 0.5 / c
  //   s    =  (rot(j, k) - rot(k, j)) * 0.5 / c
  //
  // Since we're going to pay the sqrt anyway, we perform a post normalization, which also fixes any poorly normalized input.
  // Multiply all elements by 2*c in the above, giving nc2 = -c^2
  double nc2 = ((rot(j, j) + rot(k, k)) - rot(i, i)) - 1.0;
  v[i] =  nc2;
  s    =  (rot(j, k) - rot(k, j));
  v[j] = -(rot(i, j) + rot(j, i));
  v[k] = -(rot(i, k) + rot(k, i));

  // We now have the correct result with the wrong magnitude, so normalize it:
  Real len = std::sqrt(v.squaredLength() + s * s);
  if (len > Math::eps<Real>())
  {
    v /= len;
    s /= len;
  }
  else
  {
    // The quaternion is nearly zero. Make it (0 0 0 1)
    *this = identity();
  }
}

Matrix3
Quat::toRotationMatrix() const
{
  Matrix3 out;
  toRotationMatrix(out);
  return out;
}

void
Quat::toRotationMatrix(Matrix3 & rot) const
{
  // Implementation from Watt and Watt, pg 362
  // See also http://www.flipcode.com/documents/matrfaq.html#Q54

  Quat q = unit();

  Real xx = 2 * q.x() * q.x();
  Real xy = 2 * q.x() * q.y();
  Real xz = 2 * q.x() * q.z();
  Real xw = 2 * q.x() * q.w();

  Real yy = 2 * q.y() * q.y();
  Real yz = 2 * q.y() * q.z();
  Real yw = 2 * q.y() * q.w();

  Real zz = 2 * q.z() * q.z();
  Real zw = 2 * q.z() * q.w();

  rot = Matrix3(1 - yy - zz,      xy - zw,      xz + yw,
                    xy + zw,  1 - xx - zz,      yz - xw,
                    xz - yw,      yz + xw,  1 - xx - yy);
}

Quat
Quat::fromAxisAngleRotation(Vector3 const & axis, Real angle)
{
  Quat q;
  q.real() = std::cos(angle / 2.0f);
  q.imag() = axis.unit() * std::sin(angle / 2.0f);

  return q;
}

void
Quat::toAxisAngleRotation(Vector3 & axis, double & angle) const
{
  // Decompose the quaternion into an angle and an axis
  angle = 2 * std::acos(s);

  axis = v;
  Real len = std::sqrt(1 - s * s);
  if (len > Math::eps<Real>())
    axis /= len;

  // Reduce the range of the angle

  if (angle < 0)
  {
    angle = -angle;
    axis = -axis;
  }

  while (angle > Math::twoPi())
  {
    angle -= Math::twoPi();
  }

  if (std::abs(angle) > Math::pi())
  {
    angle -= Math::twoPi();
  }

  // Make the angle positive.
  if (angle < 0)
  {
    angle = -angle;
    axis = -axis;
  }
}

Quat
Quat::slerp(Quat const & target, Real alpha, Real threshold) const
{
  // From: Game Physics -- David Eberly pg 538-540
  // Modified to include lerp for small angles, which is a common practice.
  // See also:
  // http://number-none.com/product/Understanding%20Slerp,%20Then%20Not%20Using%20It/index.html

  Quat const & quat0 = *this;
  Quat quat1 = target;

  debugAssertM(Math::fuzzyEq(quat0.squaredLength(), (Real)1), "Quat: slerp requires unit quaternions");
  debugAssertM(Math::fuzzyEq(quat1.squaredLength(), (Real)1), "Quat: slerp requires unit quaternions");

  // Angle between quaternion rotations
  Real cosphi = quat0.dot(quat1);
  if (cosphi < 0)
  {
    // Change the sign and fix the dot product; we need to
    // loop the other way to get the shortest path
    quat1 = -quat1;
    cosphi = -cosphi;
  }

  // Returns angle between 0 and pi (assumes both quaternions are unit length)
  Real phi = std::acos(cosphi);
  if (phi >= threshold)
  {
    // For large angles, slerp
    Real scale0, scale1;
    scale0 = std::sin((1 - alpha) * phi);
    scale1 = std::sin(alpha * phi);
    return ((quat0 * scale0) + (quat1 * scale1)) / std::sin(phi);
  }
  else
  {
    // For small angles, linear interpolate
    return quat0.nlerp(quat1, alpha);
  }
}

Quat
Quat::nlerp(Quat const & target, Real alpha) const
{
  Quat result = (*this) * (1 - alpha) + target * alpha;
  return result.unit();
}

Quat
Quat::operator*(Quat const & other) const
{
  // Following Watt & Watt, page 360
  Vector3 const & v1 = imag();
  Vector3 const & v2 = other.imag();
  Real            s1 = s;
  Real            s2 = other.s;
  return Quat(s1 * v2 + s2 * v1 + v1.cross(v2), s1 * s2 - v1.dot(v2));
}

Quat
Quat::unitRandom()
{
  // From "Uniform Random Rotations", Ken Shoemake, Graphics Gems III
  Real x0 = Random::common().uniform01();
  Real r1 = std::sqrt(1 - x0),
       r2 = std::sqrt(x0);
  Real t1 = (Real)Math::twoPi() * Random::common().uniform01();
  Real t2 = (Real)Math::twoPi() * Random::common().uniform01();
  Real c1 = std::cos(t1),
       s1 = std::sin(t1);
  Real c2 = std::cos(t2),
       s2 = std::sin(t2);
  return Quat(s1 * r1, c1 * r1, s2 * r2, c2 * r2);
}

} // namespace Thea
