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

/*
  ORIGINAL HEADER

  @file Quat.h

  Quaternion

  @maintainer Morgan McGuire, http://graphics.cs.williams.edu

  @created 2002-01-23
  @edited  2009-05-10
 */

#ifndef __Thea_Quat_hpp__
#define __Thea_Quat_hpp__

#include "Common.hpp"
#include "Math.hpp"
#include "Matrix3.hpp"
#include "Vector3.hpp"

namespace Thea {

/**
 * Arbitrary quaternion (not necessarily unit). Derived from the G3D library: http://g3d.sourceforge.net
 *
 * Unit quaternions are used in computer graphics to represent rotation about an axis. Any 3x3 rotation matrix can be stored
 * as a quaternion.
 *
 * A quaternion represents the sum of a real scalar and an imaginary vector: ix + jy + kz + s. A unit quaternion representing
 * a rotation by A about axis v has the form [sin(A/2)*v, cos(A/2)]. For a unit quaternion, q.conj() == q.inverse() is a
 * rotation by -A about v. -q is the same rotation as q (negate both the axis and angle).
 *
 * A non-unit quaterion q represents the same rotation as q.unitize() (Dam98 pg 28).
 *
 * Although quaternion-vector operations (eg. Quat + Vector3) are well defined, they are not supported by this class because
 * they typically are bugs when they appear in code.
 *
 * @cite Dam98 Erik B. Dam, Martin Koch, Martin Lillholm, Quaternions, Interpolation and Animation. Technical Report
 *   DIKU-TR-98/5, Department of Computer Science, University of Copenhagen, Denmark. 1998.
 */
class THEA_API Quat
{
  private:
    // Components v = (x, y, z), s = w
    //
    // q = [sin(angle / 2) * axis, cos(angle / 2)]
    //
    // In Watt & Watt's notation, v = (x, y, z), s = w
    // In the Real-Time Rendering notation, u = (x, y, z), w = w
    Vector3 v;
    Real s;

  public:
    /** Initialize to a zero degree rotation: (0, 0, 0, 1). */
    Quat() : v(0, 0, 0), s(1) {}

    /** Initialize from a rotation matrix. */
    Quat(Matrix3 const & rot);

    /** Initialize from components. */
    Quat(Real x_, Real y_, Real z_, Real w_) : v(x_, y_, z_), s(w_) {}

    /** Initialize from imaginary and real parts. */
    Quat(Vector3 const & im, Real re = 0) : v(im), s(re) {}

    /** Get the zero quaternion (all components zero). */
    static Quat const & zero() { static Quat const z(0, 0, 0, 0); return z; }

    /** Get the identity quaternion (zero rotation). */
    static Quat const & identity() { static Quat const id(0, 0, 0, 1); return id; }

    /** The real part w of the quaternion. */
    Real const & real() const
    {
      return s;
    }

    /** The real part w of the quaternion. */
    Real & real()
    {
      return s;
    }

    /** The imaginary part (x, y, z) of the quaternion. */
    Vector3 const & imag() const
    {
      return v;
    }

    /** The imaginary part (x, y, z) of the quaternion. */
    Vector3 & imag()
    {
      return v;
    }

    /** Negation operator. */
    Quat operator-() const
    {
      return Quat(-v, -s);
    }

    /** Subtraction operator. */
    Quat operator-(Quat const & other) const
    {
      return Quat(v - other.v, s - other.s);
    }

    /** Subtract-and-assign. */
    Quat & operator-=(Quat const & q)
    {
      v -= q.v;
      s -= q.s;
      return *this;
    }

    /** Addition operator. */
    Quat operator+(Quat const & q) const
    {
      return Quat(v + q.v, s + q.s);
    }

    /** Add-and-assign. */
    Quat & operator+=(Quat const & q)
    {
      v += q.v;
      s += q.s;
      return *this;
    }

    /** Return the conjugate by negating the imaginary part. */
    Quat conj() const
    {
      return Quat(-v, s);
    }

    /** Return the sum of the components. */
    Real sum() const
    {
      return v[0] + v[1] + v[2] + s;
    }

    /** Return the sum of the components. */
    Real average() const
    {
      return sum() / 4.0f;
    }

    /** Multiply by a scalar. */
    Quat operator*(Real t) const
    {
      return Quat(v * t, s * t);
    }

    /** Multiply by a scalar and assign. */
    Quat & operator*=(Real t)
    {
      v *= t;
      s *= t;
      return *this;
    }

    /** Divide by a scalar. */
    Quat operator/(Real t) const
    {
      return Quat(v / t, s / t);
    }

    /** Dot product. */
    Real dot(Quat const & other) const
    {
      return v.dot(other.v) + (s * other.s);
    }

    /**
     * Check if the components of two quaternions are approximately equal. Note: two quats can represent the same rotation and
     * not be equal.
     */
    bool fuzzyEq(Quat const & q) const
    {
      return Math::fuzzyEq(v[0], q.v[0])
          && Math::fuzzyEq(v[1], q.v[1])
          && Math::fuzzyEq(v[2], q.v[2])
          && Math::fuzzyEq(s, q.s);
    }

    /** Check if two quaternions represent the same rotation. Note: every rotation is represented by two values; q and -q. */
    bool sameRotation(Quat const & q) const
    {
      return fuzzyEq(q) || fuzzyEq(-q);
    }

    /** Create the quaternion from a rotation matrix, using the formula q = [sin(angle/2)*axis, cos(angle/2)]. */
    static Quat fromAxisAngleRotation(Vector3 const & axis, Real angle);

    /**
     * Returns the axis and angle of rotation represented by this quaternion, inverting the map q = [sin(angle/2)*axis,
     * cos(angle/2)])
     */
    void toAxisAngleRotation(Vector3 & axis, double & angle) const;

    /**
     * Returns the axis and angle of rotation represented by this quaternion, inverting the map q = [sin(angle/2)*axis,
     * cos(angle/2)])
     */
    void toAxisAngleRotation(Vector3 & axis, Real & angle) const
    {
      double d;
      toAxisAngleRotation(axis, d);
      angle = (Real)d;
    }

    /** Convert the quaternion to a rotation matrix. */
    Matrix3 toRotationMatrix() const;

    /** Convert the quaternion to a rotation matrix. */
    void toRotationMatrix(Matrix3 & rot) const;

    /**
     * Spherical linear interpolation: linear interpolation along the shortest (3D) great-circle route between two quaternions.
     * The quaternions must represent valid rotations, i.e. be of the form (sin(A/2)*axis, cos(A/2)) with unit length.
     *
     * Should correctly rotate between 0 and Pi in the right order.
     *
     * @cite Eberly Based on Game Physics -- David Eberly pg 538-540
     *
     * @param target Interpolation target.
     * @param alpha Interpolation parameter.
     * @param threshold Critical angle between between rotations at which the algorithm switches to normalized lerp, which is
     *   more numerically stable in those situations. 0.0 will always slerp.
     */
    Quat slerp(Quat const & target, Real alpha, Real threshold = 0.05f) const;

    /** Normalized linear interpolation of quaternion components. Returns a unit quaternion. */
    Quat nlerp(Quat const & target, Real alpha) const;

    /**
     * Get the inverse of the quaternion. Note that q.inverse() = q.conj() for a unit quaternion.
     *
     * @cite Dam99 page 13
     */
    Quat inverse() const
    {
      return conj() / dot(*this);
    }

    /** Quaternion multiplication (composition of rotations). Note that this does not commute. */
    Quat operator*(Quat const & other) const;

    /* "Division" for quaternions, defined as (*this) * other.inverse() */
    Quat operator/(Quat const & other) const
    {
      return (*this) * other.inverse();
    }

    /** Get the square of the length of the quaternion. */
    Real squaredLength() const { return dot(*this); }

    /** Get the length of the quaternion. */
    Real length() const { return std::sqrt(squaredLength()); }

    /** Check if the length is nearly 1. */
    bool isUnit(Real tolerance = 1e-5) const
    {
      return std::abs(squaredLength() - 1.0f) < tolerance;
    }

    /** Logarithm of the quaternion. */
    Quat log() const
    {
      if ((v[0] == 0) && (v[1] == 0) && (v[2] == 0))
      {
        if (s > 0)
        {
          return Quat(0, 0, 0, std::log(s));
        }
        else if (s < 0)
        {
          // Log of a negative number.  Multivalued, any number of the form (PI * v, ln(-q.s))
          return Quat((Real)Math::pi(), 0, 0, std::log(-s));
        }
        else
        {
          // Log of zero!
          static Real const n = Math::nan<Real>();
          return Quat(n, n, n, n);
        }
      }
      else
      {
        // Partly imaginary.
        Real imag_len = v.length();
        Real theta = std::atan2(imag_len, s);
        Real t = theta / imag_len;
        return Quat(v * t, std::log(length()));
      }
    }

    /**
     * Exponentiate the quaternion, using the formula exp(q) = [sin(A) * v, cos(A)] where q = [Av, 0]. Only defined for
     * pure-vector quaternions
     */
    Quat exp() const
    {
      alwaysAssertM(s == 0, "Quat: exp only defined for vector quaternions");

      Real imag_len = v.length();
      Vector3 u = v / imag_len;
      return Quat(std::sin(imag_len) * u, std::cos(imag_len));
    }

    /**
     * Raise this quaternion to a power. For a rotation, this is the effect of rotating p times as much as the original
     * quaternion.
     *
     * Note that q.pow(a).pow(b) = q.pow(a + b).
     *
     * @cite Dam98 pg 21
     */
    Quat pow(Real p) const
    {
      return (log() * p).exp();
    }

    /** Scale the quaternion to unit length. */
    void unitize()
    {
      *this = unit();
    }

    /** Get a unit quaternion by dividing by the length. */
    Quat unit() const
    {
      Real len = length();
      if (std::abs(len) < 32 * std::numeric_limits<Real>::min())
        return zero();
      else
        return *this / len;
    }

    /**
     * Generate a uniform random unit quaternion (a random "direction")
     *
     * @cite Shoemake From "Uniform Random Rotations", Ken Shoemake, Graphics Gems III.
     */
    static Quat unitRandom();

    // Element access
    Real x() const { return v[0]; }  ///< Get the first imaginary component.
    Real y() const { return v[1]; }  ///< Get the second imaginary component.
    Real z() const { return v[2]; }  ///< Get the third imaginary component.
    Real w() const { return s;    }  ///< Get the real component.

    //========================================================================================================================
    // 2D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XX. */
    Vector2 xx() const { return Vector2   (x(), x()); }

    /** Swizzle operator, returns YX. */
    Vector2 yx() const { return Vector2   (y(), x()); }

    /** Swizzle operator, returns ZX. */
    Vector2 zx() const { return Vector2   (z(), x()); }

    /** Swizzle operator, returns WX. */
    Vector2 wx() const { return Vector2   (w(), x()); }

    /** Swizzle operator, returns XY. */
    Vector2 xy() const { return Vector2   (x(), y()); }

    /** Swizzle operator, returns YY. */
    Vector2 yy() const { return Vector2   (y(), y()); }

    /** Swizzle operator, returns ZY. */
    Vector2 zy() const { return Vector2   (z(), y()); }

    /** Swizzle operator, returns WY. */
    Vector2 wy() const { return Vector2   (w(), y()); }

    /** Swizzle operator, returns XZ. */
    Vector2 xz() const { return Vector2   (x(), z()); }

    /** Swizzle operator, returns YZ. */
    Vector2 yz() const { return Vector2   (y(), z()); }

    /** Swizzle operator, returns ZZ. */
    Vector2 zz() const { return Vector2   (z(), z()); }

    /** Swizzle operator, returns WZ. */
    Vector2 wz() const { return Vector2   (w(), z()); }

    /** Swizzle operator, returns XW. */
    Vector2 xw() const { return Vector2   (x(), w()); }

    /** Swizzle operator, returns YW. */
    Vector2 yw() const { return Vector2   (y(), w()); }

    /** Swizzle operator, returns ZW. */
    Vector2 zw() const { return Vector2   (z(), w()); }

    /** Swizzle operator, returns WW. */
    Vector2 ww() const { return Vector2   (w(), w()); }

    //========================================================================================================================
    // 3D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XXX. */
    Vector3 xxx() const { return Vector3   (x(), x(), x()); }

    /** Swizzle operator, returns YXX. */
    Vector3 yxx() const { return Vector3   (y(), x(), x()); }

    /** Swizzle operator, returns ZXX. */
    Vector3 zxx() const { return Vector3   (z(), x(), x()); }

    /** Swizzle operator, returns WXX. */
    Vector3 wxx() const { return Vector3   (w(), x(), x()); }

    /** Swizzle operator, returns XYX. */
    Vector3 xyx() const { return Vector3   (x(), y(), x()); }

    /** Swizzle operator, returns YYX. */
    Vector3 yyx() const { return Vector3   (y(), y(), x()); }

    /** Swizzle operator, returns ZYX. */
    Vector3 zyx() const { return Vector3   (z(), y(), x()); }

    /** Swizzle operator, returns WYX. */
    Vector3 wyx() const { return Vector3   (w(), y(), x()); }

    /** Swizzle operator, returns XZX. */
    Vector3 xzx() const { return Vector3   (x(), z(), x()); }

    /** Swizzle operator, returns YZX. */
    Vector3 yzx() const { return Vector3   (y(), z(), x()); }

    /** Swizzle operator, returns ZZX. */
    Vector3 zzx() const { return Vector3   (z(), z(), x()); }

    /** Swizzle operator, returns WZX. */
    Vector3 wzx() const { return Vector3   (w(), z(), x()); }

    /** Swizzle operator, returns XWX. */
    Vector3 xwx() const { return Vector3   (x(), w(), x()); }

    /** Swizzle operator, returns YWX. */
    Vector3 ywx() const { return Vector3   (y(), w(), x()); }

    /** Swizzle operator, returns ZWX. */
    Vector3 zwx() const { return Vector3   (z(), w(), x()); }

    /** Swizzle operator, returns WWX. */
    Vector3 wwx() const { return Vector3   (w(), w(), x()); }

    /** Swizzle operator, returns XXY. */
    Vector3 xxy() const { return Vector3   (x(), x(), y()); }

    /** Swizzle operator, returns YXY. */
    Vector3 yxy() const { return Vector3   (y(), x(), y()); }

    /** Swizzle operator, returns ZXY. */
    Vector3 zxy() const { return Vector3   (z(), x(), y()); }

    /** Swizzle operator, returns WXY. */
    Vector3 wxy() const { return Vector3   (w(), x(), y()); }

    /** Swizzle operator, returns XYY. */
    Vector3 xyy() const { return Vector3   (x(), y(), y()); }

    /** Swizzle operator, returns YYY. */
    Vector3 yyy() const { return Vector3   (y(), y(), y()); }

    /** Swizzle operator, returns ZYY. */
    Vector3 zyy() const { return Vector3   (z(), y(), y()); }

    /** Swizzle operator, returns WYY. */
    Vector3 wyy() const { return Vector3   (w(), y(), y()); }

    /** Swizzle operator, returns XZY. */
    Vector3 xzy() const { return Vector3   (x(), z(), y()); }

    /** Swizzle operator, returns YZY. */
    Vector3 yzy() const { return Vector3   (y(), z(), y()); }

    /** Swizzle operator, returns ZZY. */
    Vector3 zzy() const { return Vector3   (z(), z(), y()); }

    /** Swizzle operator, returns WZY. */
    Vector3 wzy() const { return Vector3   (w(), z(), y()); }

    /** Swizzle operator, returns XWY. */
    Vector3 xwy() const { return Vector3   (x(), w(), y()); }

    /** Swizzle operator, returns YWY. */
    Vector3 ywy() const { return Vector3   (y(), w(), y()); }

    /** Swizzle operator, returns ZWY. */
    Vector3 zwy() const { return Vector3   (z(), w(), y()); }

    /** Swizzle operator, returns WWY. */
    Vector3 wwy() const { return Vector3   (w(), w(), y()); }

    /** Swizzle operator, returns XXZ. */
    Vector3 xxz() const { return Vector3   (x(), x(), z()); }

    /** Swizzle operator, returns YXZ. */
    Vector3 yxz() const { return Vector3   (y(), x(), z()); }

    /** Swizzle operator, returns ZXZ. */
    Vector3 zxz() const { return Vector3   (z(), x(), z()); }

    /** Swizzle operator, returns WXZ. */
    Vector3 wxz() const { return Vector3   (w(), x(), z()); }

    /** Swizzle operator, returns XYZ. */
    Vector3 xyz() const { return Vector3   (x(), y(), z()); }

    /** Swizzle operator, returns YYZ. */
    Vector3 yyz() const { return Vector3   (y(), y(), z()); }

    /** Swizzle operator, returns ZYZ. */
    Vector3 zyz() const { return Vector3   (z(), y(), z()); }

    /** Swizzle operator, returns WYZ. */
    Vector3 wyz() const { return Vector3   (w(), y(), z()); }

    /** Swizzle operator, returns XZZ. */
    Vector3 xzz() const { return Vector3   (x(), z(), z()); }

    /** Swizzle operator, returns YZZ. */
    Vector3 yzz() const { return Vector3   (y(), z(), z()); }

    /** Swizzle operator, returns ZZZ. */
    Vector3 zzz() const { return Vector3   (z(), z(), z()); }

    /** Swizzle operator, returns WZZ. */
    Vector3 wzz() const { return Vector3   (w(), z(), z()); }

    /** Swizzle operator, returns XWZ. */
    Vector3 xwz() const { return Vector3   (x(), w(), z()); }

    /** Swizzle operator, returns YWZ. */
    Vector3 ywz() const { return Vector3   (y(), w(), z()); }

    /** Swizzle operator, returns ZWZ. */
    Vector3 zwz() const { return Vector3   (z(), w(), z()); }

    /** Swizzle operator, returns WWZ. */
    Vector3 wwz() const { return Vector3   (w(), w(), z()); }

    /** Swizzle operator, returns XXW. */
    Vector3 xxw() const { return Vector3   (x(), x(), w()); }

    /** Swizzle operator, returns YXW. */
    Vector3 yxw() const { return Vector3   (y(), x(), w()); }

    /** Swizzle operator, returns ZXW. */
    Vector3 zxw() const { return Vector3   (z(), x(), w()); }

    /** Swizzle operator, returns WXW. */
    Vector3 wxw() const { return Vector3   (w(), x(), w()); }

    /** Swizzle operator, returns XYW. */
    Vector3 xyw() const { return Vector3   (x(), y(), w()); }

    /** Swizzle operator, returns YYW. */
    Vector3 yyw() const { return Vector3   (y(), y(), w()); }

    /** Swizzle operator, returns ZYW. */
    Vector3 zyw() const { return Vector3   (z(), y(), w()); }

    /** Swizzle operator, returns WYW. */
    Vector3 wyw() const { return Vector3   (w(), y(), w()); }

    /** Swizzle operator, returns XZW. */
    Vector3 xzw() const { return Vector3   (x(), z(), w()); }

    /** Swizzle operator, returns YZW. */
    Vector3 yzw() const { return Vector3   (y(), z(), w()); }

    /** Swizzle operator, returns ZZW. */
    Vector3 zzw() const { return Vector3   (z(), z(), w()); }

    /** Swizzle operator, returns WZW. */
    Vector3 wzw() const { return Vector3   (w(), z(), w()); }

    /** Swizzle operator, returns XWW. */
    Vector3 xww() const { return Vector3   (x(), w(), w()); }

    /** Swizzle operator, returns YWW. */
    Vector3 yww() const { return Vector3   (y(), w(), w()); }

    /** Swizzle operator, returns ZWW. */
    Vector3 zww() const { return Vector3   (z(), w(), w()); }

    /** Swizzle operator, returns WWW. */
    Vector3 www() const { return Vector3   (w(), w(), w()); }

    //========================================================================================================================
    // 4D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XXXX. */
    Vector4 xxxx() const { return Vector4   (x(), x(), x(), x()); }

    /** Swizzle operator, returns YXXX. */
    Vector4 yxxx() const { return Vector4   (y(), x(), x(), x()); }

    /** Swizzle operator, returns ZXXX. */
    Vector4 zxxx() const { return Vector4   (z(), x(), x(), x()); }

    /** Swizzle operator, returns WXXX. */
    Vector4 wxxx() const { return Vector4   (w(), x(), x(), x()); }

    /** Swizzle operator, returns XYXX. */
    Vector4 xyxx() const { return Vector4   (x(), y(), x(), x()); }

    /** Swizzle operator, returns YYXX. */
    Vector4 yyxx() const { return Vector4   (y(), y(), x(), x()); }

    /** Swizzle operator, returns ZYXX. */
    Vector4 zyxx() const { return Vector4   (z(), y(), x(), x()); }

    /** Swizzle operator, returns WYXX. */
    Vector4 wyxx() const { return Vector4   (w(), y(), x(), x()); }

    /** Swizzle operator, returns XZXX. */
    Vector4 xzxx() const { return Vector4   (x(), z(), x(), x()); }

    /** Swizzle operator, returns YZXX. */
    Vector4 yzxx() const { return Vector4   (y(), z(), x(), x()); }

    /** Swizzle operator, returns ZZXX. */
    Vector4 zzxx() const { return Vector4   (z(), z(), x(), x()); }

    /** Swizzle operator, returns WZXX. */
    Vector4 wzxx() const { return Vector4   (w(), z(), x(), x()); }

    /** Swizzle operator, returns XWXX. */
    Vector4 xwxx() const { return Vector4   (x(), w(), x(), x()); }

    /** Swizzle operator, returns YWXX. */
    Vector4 ywxx() const { return Vector4   (y(), w(), x(), x()); }

    /** Swizzle operator, returns ZWXX. */
    Vector4 zwxx() const { return Vector4   (z(), w(), x(), x()); }

    /** Swizzle operator, returns WWXX. */
    Vector4 wwxx() const { return Vector4   (w(), w(), x(), x()); }

    /** Swizzle operator, returns XXYX. */
    Vector4 xxyx() const { return Vector4   (x(), x(), y(), x()); }

    /** Swizzle operator, returns YXYX. */
    Vector4 yxyx() const { return Vector4   (y(), x(), y(), x()); }

    /** Swizzle operator, returns ZXYX. */
    Vector4 zxyx() const { return Vector4   (z(), x(), y(), x()); }

    /** Swizzle operator, returns WXYX. */
    Vector4 wxyx() const { return Vector4   (w(), x(), y(), x()); }

    /** Swizzle operator, returns XYYX. */
    Vector4 xyyx() const { return Vector4   (x(), y(), y(), x()); }

    /** Swizzle operator, returns YYYX. */
    Vector4 yyyx() const { return Vector4   (y(), y(), y(), x()); }

    /** Swizzle operator, returns ZYYX. */
    Vector4 zyyx() const { return Vector4   (z(), y(), y(), x()); }

    /** Swizzle operator, returns WYYX. */
    Vector4 wyyx() const { return Vector4   (w(), y(), y(), x()); }

    /** Swizzle operator, returns XZYX. */
    Vector4 xzyx() const { return Vector4   (x(), z(), y(), x()); }

    /** Swizzle operator, returns YZYX. */
    Vector4 yzyx() const { return Vector4   (y(), z(), y(), x()); }

    /** Swizzle operator, returns ZZYX. */
    Vector4 zzyx() const { return Vector4   (z(), z(), y(), x()); }

    /** Swizzle operator, returns WZYX. */
    Vector4 wzyx() const { return Vector4   (w(), z(), y(), x()); }

    /** Swizzle operator, returns XWYX. */
    Vector4 xwyx() const { return Vector4   (x(), w(), y(), x()); }

    /** Swizzle operator, returns YWYX. */
    Vector4 ywyx() const { return Vector4   (y(), w(), y(), x()); }

    /** Swizzle operator, returns ZWYX. */
    Vector4 zwyx() const { return Vector4   (z(), w(), y(), x()); }

    /** Swizzle operator, returns WWYX. */
    Vector4 wwyx() const { return Vector4   (w(), w(), y(), x()); }

    /** Swizzle operator, returns XXZX. */
    Vector4 xxzx() const { return Vector4   (x(), x(), z(), x()); }

    /** Swizzle operator, returns YXZX. */
    Vector4 yxzx() const { return Vector4   (y(), x(), z(), x()); }

    /** Swizzle operator, returns ZXZX. */
    Vector4 zxzx() const { return Vector4   (z(), x(), z(), x()); }

    /** Swizzle operator, returns WXZX. */
    Vector4 wxzx() const { return Vector4   (w(), x(), z(), x()); }

    /** Swizzle operator, returns XYZX. */
    Vector4 xyzx() const { return Vector4   (x(), y(), z(), x()); }

    /** Swizzle operator, returns YYZX. */
    Vector4 yyzx() const { return Vector4   (y(), y(), z(), x()); }

    /** Swizzle operator, returns ZYZX. */
    Vector4 zyzx() const { return Vector4   (z(), y(), z(), x()); }

    /** Swizzle operator, returns WYZX. */
    Vector4 wyzx() const { return Vector4   (w(), y(), z(), x()); }

    /** Swizzle operator, returns XZZX. */
    Vector4 xzzx() const { return Vector4   (x(), z(), z(), x()); }

    /** Swizzle operator, returns YZZX. */
    Vector4 yzzx() const { return Vector4   (y(), z(), z(), x()); }

    /** Swizzle operator, returns ZZZX. */
    Vector4 zzzx() const { return Vector4   (z(), z(), z(), x()); }

    /** Swizzle operator, returns WZZX. */
    Vector4 wzzx() const { return Vector4   (w(), z(), z(), x()); }

    /** Swizzle operator, returns XWZX. */
    Vector4 xwzx() const { return Vector4   (x(), w(), z(), x()); }

    /** Swizzle operator, returns YWZX. */
    Vector4 ywzx() const { return Vector4   (y(), w(), z(), x()); }

    /** Swizzle operator, returns ZWZX. */
    Vector4 zwzx() const { return Vector4   (z(), w(), z(), x()); }

    /** Swizzle operator, returns WWZX. */
    Vector4 wwzx() const { return Vector4   (w(), w(), z(), x()); }

    /** Swizzle operator, returns XXWX. */
    Vector4 xxwx() const { return Vector4   (x(), x(), w(), x()); }

    /** Swizzle operator, returns YXWX. */
    Vector4 yxwx() const { return Vector4   (y(), x(), w(), x()); }

    /** Swizzle operator, returns ZXWX. */
    Vector4 zxwx() const { return Vector4   (z(), x(), w(), x()); }

    /** Swizzle operator, returns WXWX. */
    Vector4 wxwx() const { return Vector4   (w(), x(), w(), x()); }

    /** Swizzle operator, returns XYWX. */
    Vector4 xywx() const { return Vector4   (x(), y(), w(), x()); }

    /** Swizzle operator, returns YYWX. */
    Vector4 yywx() const { return Vector4   (y(), y(), w(), x()); }

    /** Swizzle operator, returns ZYWX. */
    Vector4 zywx() const { return Vector4   (z(), y(), w(), x()); }

    /** Swizzle operator, returns WYWX. */
    Vector4 wywx() const { return Vector4   (w(), y(), w(), x()); }

    /** Swizzle operator, returns XZWX. */
    Vector4 xzwx() const { return Vector4   (x(), z(), w(), x()); }

    /** Swizzle operator, returns YZWX. */
    Vector4 yzwx() const { return Vector4   (y(), z(), w(), x()); }

    /** Swizzle operator, returns ZZWX. */
    Vector4 zzwx() const { return Vector4   (z(), z(), w(), x()); }

    /** Swizzle operator, returns WZWX. */
    Vector4 wzwx() const { return Vector4   (w(), z(), w(), x()); }

    /** Swizzle operator, returns XWWX. */
    Vector4 xwwx() const { return Vector4   (x(), w(), w(), x()); }

    /** Swizzle operator, returns YWWX. */
    Vector4 ywwx() const { return Vector4   (y(), w(), w(), x()); }

    /** Swizzle operator, returns ZWWX. */
    Vector4 zwwx() const { return Vector4   (z(), w(), w(), x()); }

    /** Swizzle operator, returns WWWX. */
    Vector4 wwwx() const { return Vector4   (w(), w(), w(), x()); }

    /** Swizzle operator, returns XXXY. */
    Vector4 xxxy() const { return Vector4   (x(), x(), x(), y()); }

    /** Swizzle operator, returns YXXY. */
    Vector4 yxxy() const { return Vector4   (y(), x(), x(), y()); }

    /** Swizzle operator, returns ZXXY. */
    Vector4 zxxy() const { return Vector4   (z(), x(), x(), y()); }

    /** Swizzle operator, returns WXXY. */
    Vector4 wxxy() const { return Vector4   (w(), x(), x(), y()); }

    /** Swizzle operator, returns XYXY. */
    Vector4 xyxy() const { return Vector4   (x(), y(), x(), y()); }

    /** Swizzle operator, returns YYXY. */
    Vector4 yyxy() const { return Vector4   (y(), y(), x(), y()); }

    /** Swizzle operator, returns ZYXY. */
    Vector4 zyxy() const { return Vector4   (z(), y(), x(), y()); }

    /** Swizzle operator, returns WYXY. */
    Vector4 wyxy() const { return Vector4   (w(), y(), x(), y()); }

    /** Swizzle operator, returns XZXY. */
    Vector4 xzxy() const { return Vector4   (x(), z(), x(), y()); }

    /** Swizzle operator, returns YZXY. */
    Vector4 yzxy() const { return Vector4   (y(), z(), x(), y()); }

    /** Swizzle operator, returns ZZXY. */
    Vector4 zzxy() const { return Vector4   (z(), z(), x(), y()); }

    /** Swizzle operator, returns WZXY. */
    Vector4 wzxy() const { return Vector4   (w(), z(), x(), y()); }

    /** Swizzle operator, returns XWXY. */
    Vector4 xwxy() const { return Vector4   (x(), w(), x(), y()); }

    /** Swizzle operator, returns YWXY. */
    Vector4 ywxy() const { return Vector4   (y(), w(), x(), y()); }

    /** Swizzle operator, returns ZWXY. */
    Vector4 zwxy() const { return Vector4   (z(), w(), x(), y()); }

    /** Swizzle operator, returns WWXY. */
    Vector4 wwxy() const { return Vector4   (w(), w(), x(), y()); }

    /** Swizzle operator, returns XXYY. */
    Vector4 xxyy() const { return Vector4   (x(), x(), y(), y()); }

    /** Swizzle operator, returns YXYY. */
    Vector4 yxyy() const { return Vector4   (y(), x(), y(), y()); }

    /** Swizzle operator, returns ZXYY. */
    Vector4 zxyy() const { return Vector4   (z(), x(), y(), y()); }

    /** Swizzle operator, returns WXYY. */
    Vector4 wxyy() const { return Vector4   (w(), x(), y(), y()); }

    /** Swizzle operator, returns XYYY. */
    Vector4 xyyy() const { return Vector4   (x(), y(), y(), y()); }

    /** Swizzle operator, returns YYYY. */
    Vector4 yyyy() const { return Vector4   (y(), y(), y(), y()); }

    /** Swizzle operator, returns ZYYY. */
    Vector4 zyyy() const { return Vector4   (z(), y(), y(), y()); }

    /** Swizzle operator, returns WYYY. */
    Vector4 wyyy() const { return Vector4   (w(), y(), y(), y()); }

    /** Swizzle operator, returns XZYY. */
    Vector4 xzyy() const { return Vector4   (x(), z(), y(), y()); }

    /** Swizzle operator, returns YZYY. */
    Vector4 yzyy() const { return Vector4   (y(), z(), y(), y()); }

    /** Swizzle operator, returns ZZYY. */
    Vector4 zzyy() const { return Vector4   (z(), z(), y(), y()); }

    /** Swizzle operator, returns WZYY. */
    Vector4 wzyy() const { return Vector4   (w(), z(), y(), y()); }

    /** Swizzle operator, returns XWYY. */
    Vector4 xwyy() const { return Vector4   (x(), w(), y(), y()); }

    /** Swizzle operator, returns YWYY. */
    Vector4 ywyy() const { return Vector4   (y(), w(), y(), y()); }

    /** Swizzle operator, returns ZWYY. */
    Vector4 zwyy() const { return Vector4   (z(), w(), y(), y()); }

    /** Swizzle operator, returns WWYY. */
    Vector4 wwyy() const { return Vector4   (w(), w(), y(), y()); }

    /** Swizzle operator, returns XXZY. */
    Vector4 xxzy() const { return Vector4   (x(), x(), z(), y()); }

    /** Swizzle operator, returns YXZY. */
    Vector4 yxzy() const { return Vector4   (y(), x(), z(), y()); }

    /** Swizzle operator, returns ZXZY. */
    Vector4 zxzy() const { return Vector4   (z(), x(), z(), y()); }

    /** Swizzle operator, returns WXZY. */
    Vector4 wxzy() const { return Vector4   (w(), x(), z(), y()); }

    /** Swizzle operator, returns XYZY. */
    Vector4 xyzy() const { return Vector4   (x(), y(), z(), y()); }

    /** Swizzle operator, returns YYZY. */
    Vector4 yyzy() const { return Vector4   (y(), y(), z(), y()); }

    /** Swizzle operator, returns ZYZY. */
    Vector4 zyzy() const { return Vector4   (z(), y(), z(), y()); }

    /** Swizzle operator, returns WYZY. */
    Vector4 wyzy() const { return Vector4   (w(), y(), z(), y()); }

    /** Swizzle operator, returns XZZY. */
    Vector4 xzzy() const { return Vector4   (x(), z(), z(), y()); }

    /** Swizzle operator, returns YZZY. */
    Vector4 yzzy() const { return Vector4   (y(), z(), z(), y()); }

    /** Swizzle operator, returns ZZZY. */
    Vector4 zzzy() const { return Vector4   (z(), z(), z(), y()); }

    /** Swizzle operator, returns WZZY. */
    Vector4 wzzy() const { return Vector4   (w(), z(), z(), y()); }

    /** Swizzle operator, returns XWZY. */
    Vector4 xwzy() const { return Vector4   (x(), w(), z(), y()); }

    /** Swizzle operator, returns YWZY. */
    Vector4 ywzy() const { return Vector4   (y(), w(), z(), y()); }

    /** Swizzle operator, returns ZWZY. */
    Vector4 zwzy() const { return Vector4   (z(), w(), z(), y()); }

    /** Swizzle operator, returns WWZY. */
    Vector4 wwzy() const { return Vector4   (w(), w(), z(), y()); }

    /** Swizzle operator, returns XXWY. */
    Vector4 xxwy() const { return Vector4   (x(), x(), w(), y()); }

    /** Swizzle operator, returns YXWY. */
    Vector4 yxwy() const { return Vector4   (y(), x(), w(), y()); }

    /** Swizzle operator, returns ZXWY. */
    Vector4 zxwy() const { return Vector4   (z(), x(), w(), y()); }

    /** Swizzle operator, returns WXWY. */
    Vector4 wxwy() const { return Vector4   (w(), x(), w(), y()); }

    /** Swizzle operator, returns XYWY. */
    Vector4 xywy() const { return Vector4   (x(), y(), w(), y()); }

    /** Swizzle operator, returns YYWY. */
    Vector4 yywy() const { return Vector4   (y(), y(), w(), y()); }

    /** Swizzle operator, returns ZYWY. */
    Vector4 zywy() const { return Vector4   (z(), y(), w(), y()); }

    /** Swizzle operator, returns WYWY. */
    Vector4 wywy() const { return Vector4   (w(), y(), w(), y()); }

    /** Swizzle operator, returns XZWY. */
    Vector4 xzwy() const { return Vector4   (x(), z(), w(), y()); }

    /** Swizzle operator, returns YZWY. */
    Vector4 yzwy() const { return Vector4   (y(), z(), w(), y()); }

    /** Swizzle operator, returns ZZWY. */
    Vector4 zzwy() const { return Vector4   (z(), z(), w(), y()); }

    /** Swizzle operator, returns WZWY. */
    Vector4 wzwy() const { return Vector4   (w(), z(), w(), y()); }

    /** Swizzle operator, returns XWWY. */
    Vector4 xwwy() const { return Vector4   (x(), w(), w(), y()); }

    /** Swizzle operator, returns YWWY. */
    Vector4 ywwy() const { return Vector4   (y(), w(), w(), y()); }

    /** Swizzle operator, returns ZWWY. */
    Vector4 zwwy() const { return Vector4   (z(), w(), w(), y()); }

    /** Swizzle operator, returns WWWY. */
    Vector4 wwwy() const { return Vector4   (w(), w(), w(), y()); }

    /** Swizzle operator, returns XXXZ. */
    Vector4 xxxz() const { return Vector4   (x(), x(), x(), z()); }

    /** Swizzle operator, returns YXXZ. */
    Vector4 yxxz() const { return Vector4   (y(), x(), x(), z()); }

    /** Swizzle operator, returns ZXXZ. */
    Vector4 zxxz() const { return Vector4   (z(), x(), x(), z()); }

    /** Swizzle operator, returns WXXZ. */
    Vector4 wxxz() const { return Vector4   (w(), x(), x(), z()); }

    /** Swizzle operator, returns XYXZ. */
    Vector4 xyxz() const { return Vector4   (x(), y(), x(), z()); }

    /** Swizzle operator, returns YYXZ. */
    Vector4 yyxz() const { return Vector4   (y(), y(), x(), z()); }

    /** Swizzle operator, returns ZYXZ. */
    Vector4 zyxz() const { return Vector4   (z(), y(), x(), z()); }

    /** Swizzle operator, returns WYXZ. */
    Vector4 wyxz() const { return Vector4   (w(), y(), x(), z()); }

    /** Swizzle operator, returns XZXZ. */
    Vector4 xzxz() const { return Vector4   (x(), z(), x(), z()); }

    /** Swizzle operator, returns YZXZ. */
    Vector4 yzxz() const { return Vector4   (y(), z(), x(), z()); }

    /** Swizzle operator, returns ZZXZ. */
    Vector4 zzxz() const { return Vector4   (z(), z(), x(), z()); }

    /** Swizzle operator, returns WZXZ. */
    Vector4 wzxz() const { return Vector4   (w(), z(), x(), z()); }

    /** Swizzle operator, returns XWXZ. */
    Vector4 xwxz() const { return Vector4   (x(), w(), x(), z()); }

    /** Swizzle operator, returns YWXZ. */
    Vector4 ywxz() const { return Vector4   (y(), w(), x(), z()); }

    /** Swizzle operator, returns ZWXZ. */
    Vector4 zwxz() const { return Vector4   (z(), w(), x(), z()); }

    /** Swizzle operator, returns WWXZ. */
    Vector4 wwxz() const { return Vector4   (w(), w(), x(), z()); }

    /** Swizzle operator, returns XXYZ. */
    Vector4 xxyz() const { return Vector4   (x(), x(), y(), z()); }

    /** Swizzle operator, returns YXYZ. */
    Vector4 yxyz() const { return Vector4   (y(), x(), y(), z()); }

    /** Swizzle operator, returns ZXYZ. */
    Vector4 zxyz() const { return Vector4   (z(), x(), y(), z()); }

    /** Swizzle operator, returns WXYZ. */
    Vector4 wxyz() const { return Vector4   (w(), x(), y(), z()); }

    /** Swizzle operator, returns XYYZ. */
    Vector4 xyyz() const { return Vector4   (x(), y(), y(), z()); }

    /** Swizzle operator, returns YYYZ. */
    Vector4 yyyz() const { return Vector4   (y(), y(), y(), z()); }

    /** Swizzle operator, returns ZYYZ. */
    Vector4 zyyz() const { return Vector4   (z(), y(), y(), z()); }

    /** Swizzle operator, returns WYYZ. */
    Vector4 wyyz() const { return Vector4   (w(), y(), y(), z()); }

    /** Swizzle operator, returns XZYZ. */
    Vector4 xzyz() const { return Vector4   (x(), z(), y(), z()); }

    /** Swizzle operator, returns YZYZ. */
    Vector4 yzyz() const { return Vector4   (y(), z(), y(), z()); }

    /** Swizzle operator, returns ZZYZ. */
    Vector4 zzyz() const { return Vector4   (z(), z(), y(), z()); }

    /** Swizzle operator, returns WZYZ. */
    Vector4 wzyz() const { return Vector4   (w(), z(), y(), z()); }

    /** Swizzle operator, returns XWYZ. */
    Vector4 xwyz() const { return Vector4   (x(), w(), y(), z()); }

    /** Swizzle operator, returns YWYZ. */
    Vector4 ywyz() const { return Vector4   (y(), w(), y(), z()); }

    /** Swizzle operator, returns ZWYZ. */
    Vector4 zwyz() const { return Vector4   (z(), w(), y(), z()); }

    /** Swizzle operator, returns WWYZ. */
    Vector4 wwyz() const { return Vector4   (w(), w(), y(), z()); }

    /** Swizzle operator, returns XXZZ. */
    Vector4 xxzz() const { return Vector4   (x(), x(), z(), z()); }

    /** Swizzle operator, returns YXZZ. */
    Vector4 yxzz() const { return Vector4   (y(), x(), z(), z()); }

    /** Swizzle operator, returns ZXZZ. */
    Vector4 zxzz() const { return Vector4   (z(), x(), z(), z()); }

    /** Swizzle operator, returns WXZZ. */
    Vector4 wxzz() const { return Vector4   (w(), x(), z(), z()); }

    /** Swizzle operator, returns XYZZ. */
    Vector4 xyzz() const { return Vector4   (x(), y(), z(), z()); }

    /** Swizzle operator, returns YYZZ. */
    Vector4 yyzz() const { return Vector4   (y(), y(), z(), z()); }

    /** Swizzle operator, returns ZYZZ. */
    Vector4 zyzz() const { return Vector4   (z(), y(), z(), z()); }

    /** Swizzle operator, returns WYZZ. */
    Vector4 wyzz() const { return Vector4   (w(), y(), z(), z()); }

    /** Swizzle operator, returns XZZZ. */
    Vector4 xzzz() const { return Vector4   (x(), z(), z(), z()); }

    /** Swizzle operator, returns YZZZ. */
    Vector4 yzzz() const { return Vector4   (y(), z(), z(), z()); }

    /** Swizzle operator, returns ZZZZ. */
    Vector4 zzzz() const { return Vector4   (z(), z(), z(), z()); }

    /** Swizzle operator, returns WZZZ. */
    Vector4 wzzz() const { return Vector4   (w(), z(), z(), z()); }

    /** Swizzle operator, returns XWZZ. */
    Vector4 xwzz() const { return Vector4   (x(), w(), z(), z()); }

    /** Swizzle operator, returns YWZZ. */
    Vector4 ywzz() const { return Vector4   (y(), w(), z(), z()); }

    /** Swizzle operator, returns ZWZZ. */
    Vector4 zwzz() const { return Vector4   (z(), w(), z(), z()); }

    /** Swizzle operator, returns WWZZ. */
    Vector4 wwzz() const { return Vector4   (w(), w(), z(), z()); }

    /** Swizzle operator, returns XXWZ. */
    Vector4 xxwz() const { return Vector4   (x(), x(), w(), z()); }

    /** Swizzle operator, returns YXWZ. */
    Vector4 yxwz() const { return Vector4   (y(), x(), w(), z()); }

    /** Swizzle operator, returns ZXWZ. */
    Vector4 zxwz() const { return Vector4   (z(), x(), w(), z()); }

    /** Swizzle operator, returns WXWZ. */
    Vector4 wxwz() const { return Vector4   (w(), x(), w(), z()); }

    /** Swizzle operator, returns XYWZ. */
    Vector4 xywz() const { return Vector4   (x(), y(), w(), z()); }

    /** Swizzle operator, returns YYWZ. */
    Vector4 yywz() const { return Vector4   (y(), y(), w(), z()); }

    /** Swizzle operator, returns ZYWZ. */
    Vector4 zywz() const { return Vector4   (z(), y(), w(), z()); }

    /** Swizzle operator, returns WYWZ. */
    Vector4 wywz() const { return Vector4   (w(), y(), w(), z()); }

    /** Swizzle operator, returns XZWZ. */
    Vector4 xzwz() const { return Vector4   (x(), z(), w(), z()); }

    /** Swizzle operator, returns YZWZ. */
    Vector4 yzwz() const { return Vector4   (y(), z(), w(), z()); }

    /** Swizzle operator, returns ZZWZ. */
    Vector4 zzwz() const { return Vector4   (z(), z(), w(), z()); }

    /** Swizzle operator, returns WZWZ. */
    Vector4 wzwz() const { return Vector4   (w(), z(), w(), z()); }

    /** Swizzle operator, returns XWWZ. */
    Vector4 xwwz() const { return Vector4   (x(), w(), w(), z()); }

    /** Swizzle operator, returns YWWZ. */
    Vector4 ywwz() const { return Vector4   (y(), w(), w(), z()); }

    /** Swizzle operator, returns ZWWZ. */
    Vector4 zwwz() const { return Vector4   (z(), w(), w(), z()); }

    /** Swizzle operator, returns WWWZ. */
    Vector4 wwwz() const { return Vector4   (w(), w(), w(), z()); }

    /** Swizzle operator, returns XXXW. */
    Vector4 xxxw() const { return Vector4   (x(), x(), x(), w()); }

    /** Swizzle operator, returns YXXW. */
    Vector4 yxxw() const { return Vector4   (y(), x(), x(), w()); }

    /** Swizzle operator, returns ZXXW. */
    Vector4 zxxw() const { return Vector4   (z(), x(), x(), w()); }

    /** Swizzle operator, returns WXXW. */
    Vector4 wxxw() const { return Vector4   (w(), x(), x(), w()); }

    /** Swizzle operator, returns XYXW. */
    Vector4 xyxw() const { return Vector4   (x(), y(), x(), w()); }

    /** Swizzle operator, returns YYXW. */
    Vector4 yyxw() const { return Vector4   (y(), y(), x(), w()); }

    /** Swizzle operator, returns ZYXW. */
    Vector4 zyxw() const { return Vector4   (z(), y(), x(), w()); }

    /** Swizzle operator, returns WYXW. */
    Vector4 wyxw() const { return Vector4   (w(), y(), x(), w()); }

    /** Swizzle operator, returns XZXW. */
    Vector4 xzxw() const { return Vector4   (x(), z(), x(), w()); }

    /** Swizzle operator, returns YZXW. */
    Vector4 yzxw() const { return Vector4   (y(), z(), x(), w()); }

    /** Swizzle operator, returns ZZXW. */
    Vector4 zzxw() const { return Vector4   (z(), z(), x(), w()); }

    /** Swizzle operator, returns WZXW. */
    Vector4 wzxw() const { return Vector4   (w(), z(), x(), w()); }

    /** Swizzle operator, returns XWXW. */
    Vector4 xwxw() const { return Vector4   (x(), w(), x(), w()); }

    /** Swizzle operator, returns YWXW. */
    Vector4 ywxw() const { return Vector4   (y(), w(), x(), w()); }

    /** Swizzle operator, returns ZWXW. */
    Vector4 zwxw() const { return Vector4   (z(), w(), x(), w()); }

    /** Swizzle operator, returns WWXW. */
    Vector4 wwxw() const { return Vector4   (w(), w(), x(), w()); }

    /** Swizzle operator, returns XXYW. */
    Vector4 xxyw() const { return Vector4   (x(), x(), y(), w()); }

    /** Swizzle operator, returns YXYW. */
    Vector4 yxyw() const { return Vector4   (y(), x(), y(), w()); }

    /** Swizzle operator, returns ZXYW. */
    Vector4 zxyw() const { return Vector4   (z(), x(), y(), w()); }

    /** Swizzle operator, returns WXYW. */
    Vector4 wxyw() const { return Vector4   (w(), x(), y(), w()); }

    /** Swizzle operator, returns XYYW. */
    Vector4 xyyw() const { return Vector4   (x(), y(), y(), w()); }

    /** Swizzle operator, returns YYYW. */
    Vector4 yyyw() const { return Vector4   (y(), y(), y(), w()); }

    /** Swizzle operator, returns ZYYW. */
    Vector4 zyyw() const { return Vector4   (z(), y(), y(), w()); }

    /** Swizzle operator, returns WYYW. */
    Vector4 wyyw() const { return Vector4   (w(), y(), y(), w()); }

    /** Swizzle operator, returns XZYW. */
    Vector4 xzyw() const { return Vector4   (x(), z(), y(), w()); }

    /** Swizzle operator, returns YZYW. */
    Vector4 yzyw() const { return Vector4   (y(), z(), y(), w()); }

    /** Swizzle operator, returns ZZYW. */
    Vector4 zzyw() const { return Vector4   (z(), z(), y(), w()); }

    /** Swizzle operator, returns WZYW. */
    Vector4 wzyw() const { return Vector4   (w(), z(), y(), w()); }

    /** Swizzle operator, returns XWYW. */
    Vector4 xwyw() const { return Vector4   (x(), w(), y(), w()); }

    /** Swizzle operator, returns YWYW. */
    Vector4 ywyw() const { return Vector4   (y(), w(), y(), w()); }

    /** Swizzle operator, returns ZWYW. */
    Vector4 zwyw() const { return Vector4   (z(), w(), y(), w()); }

    /** Swizzle operator, returns WWYW. */
    Vector4 wwyw() const { return Vector4   (w(), w(), y(), w()); }

    /** Swizzle operator, returns XXZW. */
    Vector4 xxzw() const { return Vector4   (x(), x(), z(), w()); }

    /** Swizzle operator, returns YXZW. */
    Vector4 yxzw() const { return Vector4   (y(), x(), z(), w()); }

    /** Swizzle operator, returns ZXZW. */
    Vector4 zxzw() const { return Vector4   (z(), x(), z(), w()); }

    /** Swizzle operator, returns WXZW. */
    Vector4 wxzw() const { return Vector4   (w(), x(), z(), w()); }

    /** Swizzle operator, returns XYZW. */
    Vector4 xyzw() const { return Vector4   (x(), y(), z(), w()); }

    /** Swizzle operator, returns YYZW. */
    Vector4 yyzw() const { return Vector4   (y(), y(), z(), w()); }

    /** Swizzle operator, returns ZYZW. */
    Vector4 zyzw() const { return Vector4   (z(), y(), z(), w()); }

    /** Swizzle operator, returns WYZW. */
    Vector4 wyzw() const { return Vector4   (w(), y(), z(), w()); }

    /** Swizzle operator, returns XZZW. */
    Vector4 xzzw() const { return Vector4   (x(), z(), z(), w()); }

    /** Swizzle operator, returns YZZW. */
    Vector4 yzzw() const { return Vector4   (y(), z(), z(), w()); }

    /** Swizzle operator, returns ZZZW. */
    Vector4 zzzw() const { return Vector4   (z(), z(), z(), w()); }

    /** Swizzle operator, returns WZZW. */
    Vector4 wzzw() const { return Vector4   (w(), z(), z(), w()); }

    /** Swizzle operator, returns XWZW. */
    Vector4 xwzw() const { return Vector4   (x(), w(), z(), w()); }

    /** Swizzle operator, returns YWZW. */
    Vector4 ywzw() const { return Vector4   (y(), w(), z(), w()); }

    /** Swizzle operator, returns ZWZW. */
    Vector4 zwzw() const { return Vector4   (z(), w(), z(), w()); }

    /** Swizzle operator, returns WWZW. */
    Vector4 wwzw() const { return Vector4   (w(), w(), z(), w()); }

    /** Swizzle operator, returns XXWW. */
    Vector4 xxww() const { return Vector4   (x(), x(), w(), w()); }

    /** Swizzle operator, returns YXWW. */
    Vector4 yxww() const { return Vector4   (y(), x(), w(), w()); }

    /** Swizzle operator, returns ZXWW. */
    Vector4 zxww() const { return Vector4   (z(), x(), w(), w()); }

    /** Swizzle operator, returns WXWW. */
    Vector4 wxww() const { return Vector4   (w(), x(), w(), w()); }

    /** Swizzle operator, returns XYWW. */
    Vector4 xyww() const { return Vector4   (x(), y(), w(), w()); }

    /** Swizzle operator, returns YYWW. */
    Vector4 yyww() const { return Vector4   (y(), y(), w(), w()); }

    /** Swizzle operator, returns ZYWW. */
    Vector4 zyww() const { return Vector4   (z(), y(), w(), w()); }

    /** Swizzle operator, returns WYWW. */
    Vector4 wyww() const { return Vector4   (w(), y(), w(), w()); }

    /** Swizzle operator, returns XZWW. */
    Vector4 xzww() const { return Vector4   (x(), z(), w(), w()); }

    /** Swizzle operator, returns YZWW. */
    Vector4 yzww() const { return Vector4   (y(), z(), w(), w()); }

    /** Swizzle operator, returns ZZWW. */
    Vector4 zzww() const { return Vector4   (z(), z(), w(), w()); }

    /** Swizzle operator, returns WZWW. */
    Vector4 wzww() const { return Vector4   (w(), z(), w(), w()); }

    /** Swizzle operator, returns XWWW. */
    Vector4 xwww() const { return Vector4   (x(), w(), w(), w()); }

    /** Swizzle operator, returns YWWW. */
    Vector4 ywww() const { return Vector4   (y(), w(), w(), w()); }

    /** Swizzle operator, returns ZWWW. */
    Vector4 zwww() const { return Vector4   (z(), w(), w(), w()); }

    /** Swizzle operator, returns WWWW. */
    Vector4 wwww() const { return Vector4   (w(), w(), w(), w()); }

}; // class Quat

/**
 * Multiply a quaternion by a scalar.
 *
 * @cite Watt Based on Watt & Watt, page 360
 */
inline Quat
operator*(Real s, Quat const & q)
{
  return q * s;
}

} // namespace Thea

#endif
