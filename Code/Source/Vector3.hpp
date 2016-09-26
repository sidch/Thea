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

#ifndef __Thea_Vector3_hpp__
#define __Thea_Vector3_hpp__

#include "Vector2.hpp"
#include "VectorN.hpp"

namespace Thea {

/** 3-dimensional vectors on a field T. */
template <typename T>
class /* THEA_API */ VectorN<3, T> : public Internal::VectorNBase<3, T>
{
  private:
    typedef Internal::VectorNBase<3, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    VectorN() {}

    /** Initialize all components to a single value. */
    explicit VectorN(T const & fill_value) : BaseT(fill_value) {}

    /** Copy constructor. */
    template <typename U> VectorN(VectorN<3, U> const & src) : BaseT(src) {}

    /** Initialize from a column matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<3, 1, U> const & src) : BaseT(src) {}

    /** Initialize from a row matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<1, 3, U> const & src) : BaseT(src) {}

    /** Initialize all components of the vector. */
    VectorN(T const & x_, T const & y_, T const & z_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
      (*this)[2] = z_;
    }

    /** Initialize from a 2D vector (XY) plus a Z component. */
    VectorN(VectorN<2, T> const & xy_, T const & z_)
    {
      (*this)[0] = xy_[0];
      (*this)[1] = xy_[1];
      (*this)[2] = z_;
    }

    /** Set all elements of the vector. */
    void set(T const & x_, T const & y_, T const & z_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
      (*this)[2] = z_;
    }

    /** Cross-product of this vector with another. */
    VectorN cross(VectorN const & rhs) const
    {
      return VectorN((*this)[1] * rhs[2] - (*this)[2] * rhs[1],
                     (*this)[2] * rhs[0] - (*this)[0] * rhs[2],
                     (*this)[0] * rhs[1] - (*this)[1] * rhs[0]);
    }

    /** Get a unit vector perpendicular to this one. */
    VectorN getOrthogonalDirection() const
    {
      return ((this->maxAbsAxis() == 0) ? VectorN(y(), -x(), 0) : VectorN(0, z(), -y())).unit();
    }

    /**
     * Get two unit vectors perpendicular to this one and to each other, forming a right-handed orthonormal basis
     * (\a u, \a v, this->unit()). In other words, if this is the Z axis of the local frame, then the function returns the X and
     * Y axes.
     */
    void createOrthonormalBasis(VectorN & u, VectorN & v) const
    {
      u = getOrthogonalDirection();
      v = this->cross(u).unit();
    }

    /** Get a unit vector along positive X. */
    static VectorN const & unitX() { static VectorN const ux(1, 0, 0); return ux; }

    /** Get a unit vector along positive Y. */
    static VectorN const & unitY() { static VectorN const uy(0, 1, 0); return uy; }

    /** Get a unit vector along positive Z. */
    static VectorN const & unitZ() { static VectorN const uz(0, 0, 1); return uz; }

    /** Get the X coordinate. */
    T const & x() const { return (*this)[0]; }

    /** Get the X coordinate. */
    T & x() { return (*this)[0]; }

    /** Get the Y coordinate. */
    T const & y() const { return (*this)[1]; }

    /** Get the Y coordinate. */
    T & y() { return (*this)[1]; }

    /** Get the Z coordinate. */
    T const & z() const { return (*this)[2]; }

    /** Get the Z coordinate. */
    T & z() { return (*this)[2]; }

    //========================================================================================================================
    // 2D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XX. */
    VectorN<2, T> xx() const { return VectorN<2, T>   (x(), x()); }

    /** Swizzle operator, returns YX. */
    VectorN<2, T> yx() const { return VectorN<2, T>   (y(), x()); }

    /** Swizzle operator, returns ZX. */
    VectorN<2, T> zx() const { return VectorN<2, T>   (z(), x()); }

    /** Swizzle operator, returns XY. */
    VectorN<2, T> xy() const { return VectorN<2, T>   (x(), y()); }

    /** Swizzle operator, returns YY. */
    VectorN<2, T> yy() const { return VectorN<2, T>   (y(), y()); }

    /** Swizzle operator, returns ZY. */
    VectorN<2, T> zy() const { return VectorN<2, T>   (z(), y()); }

    /** Swizzle operator, returns XZ. */
    VectorN<2, T> xz() const { return VectorN<2, T>   (x(), z()); }

    /** Swizzle operator, returns YZ. */
    VectorN<2, T> yz() const { return VectorN<2, T>   (y(), z()); }

    /** Swizzle operator, returns ZZ. */
    VectorN<2, T> zz() const { return VectorN<2, T>   (z(), z()); }

    //========================================================================================================================
    // 3D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XXX. */
    VectorN xxx() const { return VectorN   (x(), x(), x()); }

    /** Swizzle operator, returns YXX. */
    VectorN yxx() const { return VectorN   (y(), x(), x()); }

    /** Swizzle operator, returns ZXX. */
    VectorN zxx() const { return VectorN   (z(), x(), x()); }

    /** Swizzle operator, returns XYX. */
    VectorN xyx() const { return VectorN   (x(), y(), x()); }

    /** Swizzle operator, returns YYX. */
    VectorN yyx() const { return VectorN   (y(), y(), x()); }

    /** Swizzle operator, returns ZYX. */
    VectorN zyx() const { return VectorN   (z(), y(), x()); }

    /** Swizzle operator, returns XZX. */
    VectorN xzx() const { return VectorN   (x(), z(), x()); }

    /** Swizzle operator, returns YZX. */
    VectorN yzx() const { return VectorN   (y(), z(), x()); }

    /** Swizzle operator, returns ZZX. */
    VectorN zzx() const { return VectorN   (z(), z(), x()); }

    /** Swizzle operator, returns XXY. */
    VectorN xxy() const { return VectorN   (x(), x(), y()); }

    /** Swizzle operator, returns YXY. */
    VectorN yxy() const { return VectorN   (y(), x(), y()); }

    /** Swizzle operator, returns ZXY. */
    VectorN zxy() const { return VectorN   (z(), x(), y()); }

    /** Swizzle operator, returns XYY. */
    VectorN xyy() const { return VectorN   (x(), y(), y()); }

    /** Swizzle operator, returns YYY. */
    VectorN yyy() const { return VectorN   (y(), y(), y()); }

    /** Swizzle operator, returns ZYY. */
    VectorN zyy() const { return VectorN   (z(), y(), y()); }

    /** Swizzle operator, returns XZY. */
    VectorN xzy() const { return VectorN   (x(), z(), y()); }

    /** Swizzle operator, returns YZY. */
    VectorN yzy() const { return VectorN   (y(), z(), y()); }

    /** Swizzle operator, returns ZZY. */
    VectorN zzy() const { return VectorN   (z(), z(), y()); }

    /** Swizzle operator, returns XXZ. */
    VectorN xxz() const { return VectorN   (x(), x(), z()); }

    /** Swizzle operator, returns YXZ. */
    VectorN yxz() const { return VectorN   (y(), x(), z()); }

    /** Swizzle operator, returns ZXZ. */
    VectorN zxz() const { return VectorN   (z(), x(), z()); }

    /** Swizzle operator, returns XYZ. */
    VectorN xyz() const { return VectorN   (x(), y(), z()); }

    /** Swizzle operator, returns YYZ. */
    VectorN yyz() const { return VectorN   (y(), y(), z()); }

    /** Swizzle operator, returns ZYZ. */
    VectorN zyz() const { return VectorN   (z(), y(), z()); }

    /** Swizzle operator, returns XZZ. */
    VectorN xzz() const { return VectorN   (x(), z(), z()); }

    /** Swizzle operator, returns YZZ. */
    VectorN yzz() const { return VectorN   (y(), z(), z()); }

    /** Swizzle operator, returns ZZZ. */
    VectorN zzz() const { return VectorN   (z(), z(), z()); }

}; // class VectorN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API VectorN<3, Real>;
#endif

/** The default 3D real vector class. */
typedef VectorN<3, Real> Vector3;

} // namespace Thea

#endif
