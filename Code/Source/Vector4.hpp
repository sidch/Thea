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

#ifndef __Thea_Vector4_hpp__
#define __Thea_Vector4_hpp__

#include "Vector2.hpp"
#include "Vector3.hpp"
#include "VectorN.hpp"

namespace Thea {

/** 4-dimensional vectors on a field T. */
template <typename T>
class /* THEA_API */ VectorN<4, T> : public Internal::VectorNBase<4, T>
{
  private:
    typedef Internal::VectorNBase<4, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    VectorN() {}

    /** Initialize all components to a single value. */
    explicit VectorN(T const & fill_value) : BaseT(fill_value) {}

    /** Copy constructor. */
    template <typename U> VectorN(VectorN<4, U> const & src) : BaseT(src) {}

    /** Initialize from a column matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<4, 1, U> const & src) : BaseT(src) {}

    /** Initialize from a row matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<1, 4, U> const & src) : BaseT(src) {}

    /** Initialize all components of the vector. */
    VectorN(T const & x_, T const & y_, T const & z_, T const & w_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
      (*this)[2] = z_;
      (*this)[3] = w_;
    }

    /** Initialize from a 3D vector (XYZ) plus a W component. */
    VectorN(VectorN<3, T> const & xyz_, T const & w_)
    {
      (*this)[0] = xyz_[0];
      (*this)[1] = xyz_[1];
      (*this)[2] = xyz_[2];
      (*this)[3] = w_;
    }

    /** Set all elements of the vector. */
    void set(T const & x_, T const & y_, T const & z_, T const & w_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
      (*this)[2] = z_;
      (*this)[3] = w_;
    }

    /** Get a unit vector along positive X. */
    static VectorN const & unitX() { static VectorN const ux(1, 0, 0, 0); return ux; }

    /** Get a unit vector along positive Y. */
    static VectorN const & unitY() { static VectorN const uy(0, 1, 0, 0); return uy; }

    /** Get a unit vector along positive Z. */
    static VectorN const & unitZ() { static VectorN const uz(0, 0, 1, 0); return uz; }

    /** Get a unit vector along positive W. */
    static VectorN const & unitW() { static VectorN const uw(0, 0, 0, 1); return uw; }

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

    /** Get the W coordinate. */
    T const & w() const { return (*this)[3]; }

    /** Get the W coordinate. */
    T & w() { return (*this)[3]; }

    //========================================================================================================================
    // 2D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XX. */
    VectorN<2, T> xx() const { return VectorN<2, T>   (x(), x()); }

    /** Swizzle operator, returns YX. */
    VectorN<2, T> yx() const { return VectorN<2, T>   (y(), x()); }

    /** Swizzle operator, returns ZX. */
    VectorN<2, T> zx() const { return VectorN<2, T>   (z(), x()); }

    /** Swizzle operator, returns WX. */
    VectorN<2, T> wx() const { return VectorN<2, T>   (w(), x()); }

    /** Swizzle operator, returns XY. */
    VectorN<2, T> xy() const { return VectorN<2, T>   (x(), y()); }

    /** Swizzle operator, returns YY. */
    VectorN<2, T> yy() const { return VectorN<2, T>   (y(), y()); }

    /** Swizzle operator, returns ZY. */
    VectorN<2, T> zy() const { return VectorN<2, T>   (z(), y()); }

    /** Swizzle operator, returns WY. */
    VectorN<2, T> wy() const { return VectorN<2, T>   (w(), y()); }

    /** Swizzle operator, returns XZ. */
    VectorN<2, T> xz() const { return VectorN<2, T>   (x(), z()); }

    /** Swizzle operator, returns YZ. */
    VectorN<2, T> yz() const { return VectorN<2, T>   (y(), z()); }

    /** Swizzle operator, returns ZZ. */
    VectorN<2, T> zz() const { return VectorN<2, T>   (z(), z()); }

    /** Swizzle operator, returns WZ. */
    VectorN<2, T> wz() const { return VectorN<2, T>   (w(), z()); }

    /** Swizzle operator, returns XW. */
    VectorN<2, T> xw() const { return VectorN<2, T>   (x(), w()); }

    /** Swizzle operator, returns YW. */
    VectorN<2, T> yw() const { return VectorN<2, T>   (y(), w()); }

    /** Swizzle operator, returns ZW. */
    VectorN<2, T> zw() const { return VectorN<2, T>   (z(), w()); }

    /** Swizzle operator, returns WW. */
    VectorN<2, T> ww() const { return VectorN<2, T>   (w(), w()); }

    //========================================================================================================================
    // 3D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XXX. */
    VectorN<3, T> xxx() const { return VectorN<3, T>   (x(), x(), x()); }

    /** Swizzle operator, returns YXX. */
    VectorN<3, T> yxx() const { return VectorN<3, T>   (y(), x(), x()); }

    /** Swizzle operator, returns ZXX. */
    VectorN<3, T> zxx() const { return VectorN<3, T>   (z(), x(), x()); }

    /** Swizzle operator, returns WXX. */
    VectorN<3, T> wxx() const { return VectorN<3, T>   (w(), x(), x()); }

    /** Swizzle operator, returns XYX. */
    VectorN<3, T> xyx() const { return VectorN<3, T>   (x(), y(), x()); }

    /** Swizzle operator, returns YYX. */
    VectorN<3, T> yyx() const { return VectorN<3, T>   (y(), y(), x()); }

    /** Swizzle operator, returns ZYX. */
    VectorN<3, T> zyx() const { return VectorN<3, T>   (z(), y(), x()); }

    /** Swizzle operator, returns WYX. */
    VectorN<3, T> wyx() const { return VectorN<3, T>   (w(), y(), x()); }

    /** Swizzle operator, returns XZX. */
    VectorN<3, T> xzx() const { return VectorN<3, T>   (x(), z(), x()); }

    /** Swizzle operator, returns YZX. */
    VectorN<3, T> yzx() const { return VectorN<3, T>   (y(), z(), x()); }

    /** Swizzle operator, returns ZZX. */
    VectorN<3, T> zzx() const { return VectorN<3, T>   (z(), z(), x()); }

    /** Swizzle operator, returns WZX. */
    VectorN<3, T> wzx() const { return VectorN<3, T>   (w(), z(), x()); }

    /** Swizzle operator, returns XWX. */
    VectorN<3, T> xwx() const { return VectorN<3, T>   (x(), w(), x()); }

    /** Swizzle operator, returns YWX. */
    VectorN<3, T> ywx() const { return VectorN<3, T>   (y(), w(), x()); }

    /** Swizzle operator, returns ZWX. */
    VectorN<3, T> zwx() const { return VectorN<3, T>   (z(), w(), x()); }

    /** Swizzle operator, returns WWX. */
    VectorN<3, T> wwx() const { return VectorN<3, T>   (w(), w(), x()); }

    /** Swizzle operator, returns XXY. */
    VectorN<3, T> xxy() const { return VectorN<3, T>   (x(), x(), y()); }

    /** Swizzle operator, returns YXY. */
    VectorN<3, T> yxy() const { return VectorN<3, T>   (y(), x(), y()); }

    /** Swizzle operator, returns ZXY. */
    VectorN<3, T> zxy() const { return VectorN<3, T>   (z(), x(), y()); }

    /** Swizzle operator, returns WXY. */
    VectorN<3, T> wxy() const { return VectorN<3, T>   (w(), x(), y()); }

    /** Swizzle operator, returns XYY. */
    VectorN<3, T> xyy() const { return VectorN<3, T>   (x(), y(), y()); }

    /** Swizzle operator, returns YYY. */
    VectorN<3, T> yyy() const { return VectorN<3, T>   (y(), y(), y()); }

    /** Swizzle operator, returns ZYY. */
    VectorN<3, T> zyy() const { return VectorN<3, T>   (z(), y(), y()); }

    /** Swizzle operator, returns WYY. */
    VectorN<3, T> wyy() const { return VectorN<3, T>   (w(), y(), y()); }

    /** Swizzle operator, returns XZY. */
    VectorN<3, T> xzy() const { return VectorN<3, T>   (x(), z(), y()); }

    /** Swizzle operator, returns YZY. */
    VectorN<3, T> yzy() const { return VectorN<3, T>   (y(), z(), y()); }

    /** Swizzle operator, returns ZZY. */
    VectorN<3, T> zzy() const { return VectorN<3, T>   (z(), z(), y()); }

    /** Swizzle operator, returns WZY. */
    VectorN<3, T> wzy() const { return VectorN<3, T>   (w(), z(), y()); }

    /** Swizzle operator, returns XWY. */
    VectorN<3, T> xwy() const { return VectorN<3, T>   (x(), w(), y()); }

    /** Swizzle operator, returns YWY. */
    VectorN<3, T> ywy() const { return VectorN<3, T>   (y(), w(), y()); }

    /** Swizzle operator, returns ZWY. */
    VectorN<3, T> zwy() const { return VectorN<3, T>   (z(), w(), y()); }

    /** Swizzle operator, returns WWY. */
    VectorN<3, T> wwy() const { return VectorN<3, T>   (w(), w(), y()); }

    /** Swizzle operator, returns XXZ. */
    VectorN<3, T> xxz() const { return VectorN<3, T>   (x(), x(), z()); }

    /** Swizzle operator, returns YXZ. */
    VectorN<3, T> yxz() const { return VectorN<3, T>   (y(), x(), z()); }

    /** Swizzle operator, returns ZXZ. */
    VectorN<3, T> zxz() const { return VectorN<3, T>   (z(), x(), z()); }

    /** Swizzle operator, returns WXZ. */
    VectorN<3, T> wxz() const { return VectorN<3, T>   (w(), x(), z()); }

    /** Swizzle operator, returns XYZ. */
    VectorN<3, T> xyz() const { return VectorN<3, T>   (x(), y(), z()); }

    /** Swizzle operator, returns YYZ. */
    VectorN<3, T> yyz() const { return VectorN<3, T>   (y(), y(), z()); }

    /** Swizzle operator, returns ZYZ. */
    VectorN<3, T> zyz() const { return VectorN<3, T>   (z(), y(), z()); }

    /** Swizzle operator, returns WYZ. */
    VectorN<3, T> wyz() const { return VectorN<3, T>   (w(), y(), z()); }

    /** Swizzle operator, returns XZZ. */
    VectorN<3, T> xzz() const { return VectorN<3, T>   (x(), z(), z()); }

    /** Swizzle operator, returns YZZ. */
    VectorN<3, T> yzz() const { return VectorN<3, T>   (y(), z(), z()); }

    /** Swizzle operator, returns ZZZ. */
    VectorN<3, T> zzz() const { return VectorN<3, T>   (z(), z(), z()); }

    /** Swizzle operator, returns WZZ. */
    VectorN<3, T> wzz() const { return VectorN<3, T>   (w(), z(), z()); }

    /** Swizzle operator, returns XWZ. */
    VectorN<3, T> xwz() const { return VectorN<3, T>   (x(), w(), z()); }

    /** Swizzle operator, returns YWZ. */
    VectorN<3, T> ywz() const { return VectorN<3, T>   (y(), w(), z()); }

    /** Swizzle operator, returns ZWZ. */
    VectorN<3, T> zwz() const { return VectorN<3, T>   (z(), w(), z()); }

    /** Swizzle operator, returns WWZ. */
    VectorN<3, T> wwz() const { return VectorN<3, T>   (w(), w(), z()); }

    /** Swizzle operator, returns XXW. */
    VectorN<3, T> xxw() const { return VectorN<3, T>   (x(), x(), w()); }

    /** Swizzle operator, returns YXW. */
    VectorN<3, T> yxw() const { return VectorN<3, T>   (y(), x(), w()); }

    /** Swizzle operator, returns ZXW. */
    VectorN<3, T> zxw() const { return VectorN<3, T>   (z(), x(), w()); }

    /** Swizzle operator, returns WXW. */
    VectorN<3, T> wxw() const { return VectorN<3, T>   (w(), x(), w()); }

    /** Swizzle operator, returns XYW. */
    VectorN<3, T> xyw() const { return VectorN<3, T>   (x(), y(), w()); }

    /** Swizzle operator, returns YYW. */
    VectorN<3, T> yyw() const { return VectorN<3, T>   (y(), y(), w()); }

    /** Swizzle operator, returns ZYW. */
    VectorN<3, T> zyw() const { return VectorN<3, T>   (z(), y(), w()); }

    /** Swizzle operator, returns WYW. */
    VectorN<3, T> wyw() const { return VectorN<3, T>   (w(), y(), w()); }

    /** Swizzle operator, returns XZW. */
    VectorN<3, T> xzw() const { return VectorN<3, T>   (x(), z(), w()); }

    /** Swizzle operator, returns YZW. */
    VectorN<3, T> yzw() const { return VectorN<3, T>   (y(), z(), w()); }

    /** Swizzle operator, returns ZZW. */
    VectorN<3, T> zzw() const { return VectorN<3, T>   (z(), z(), w()); }

    /** Swizzle operator, returns WZW. */
    VectorN<3, T> wzw() const { return VectorN<3, T>   (w(), z(), w()); }

    /** Swizzle operator, returns XWW. */
    VectorN<3, T> xww() const { return VectorN<3, T>   (x(), w(), w()); }

    /** Swizzle operator, returns YWW. */
    VectorN<3, T> yww() const { return VectorN<3, T>   (y(), w(), w()); }

    /** Swizzle operator, returns ZWW. */
    VectorN<3, T> zww() const { return VectorN<3, T>   (z(), w(), w()); }

    /** Swizzle operator, returns WWW. */
    VectorN<3, T> www() const { return VectorN<3, T>   (w(), w(), w()); }

    //========================================================================================================================
    // 4D swizzle operators.
    //========================================================================================================================

    /** Swizzle operator, returns XXXX. */
    VectorN xxxx() const { return VectorN   (x(), x(), x(), x()); }

    /** Swizzle operator, returns YXXX. */
    VectorN yxxx() const { return VectorN   (y(), x(), x(), x()); }

    /** Swizzle operator, returns ZXXX. */
    VectorN zxxx() const { return VectorN   (z(), x(), x(), x()); }

    /** Swizzle operator, returns WXXX. */
    VectorN wxxx() const { return VectorN   (w(), x(), x(), x()); }

    /** Swizzle operator, returns XYXX. */
    VectorN xyxx() const { return VectorN   (x(), y(), x(), x()); }

    /** Swizzle operator, returns YYXX. */
    VectorN yyxx() const { return VectorN   (y(), y(), x(), x()); }

    /** Swizzle operator, returns ZYXX. */
    VectorN zyxx() const { return VectorN   (z(), y(), x(), x()); }

    /** Swizzle operator, returns WYXX. */
    VectorN wyxx() const { return VectorN   (w(), y(), x(), x()); }

    /** Swizzle operator, returns XZXX. */
    VectorN xzxx() const { return VectorN   (x(), z(), x(), x()); }

    /** Swizzle operator, returns YZXX. */
    VectorN yzxx() const { return VectorN   (y(), z(), x(), x()); }

    /** Swizzle operator, returns ZZXX. */
    VectorN zzxx() const { return VectorN   (z(), z(), x(), x()); }

    /** Swizzle operator, returns WZXX. */
    VectorN wzxx() const { return VectorN   (w(), z(), x(), x()); }

    /** Swizzle operator, returns XWXX. */
    VectorN xwxx() const { return VectorN   (x(), w(), x(), x()); }

    /** Swizzle operator, returns YWXX. */
    VectorN ywxx() const { return VectorN   (y(), w(), x(), x()); }

    /** Swizzle operator, returns ZWXX. */
    VectorN zwxx() const { return VectorN   (z(), w(), x(), x()); }

    /** Swizzle operator, returns WWXX. */
    VectorN wwxx() const { return VectorN   (w(), w(), x(), x()); }

    /** Swizzle operator, returns XXYX. */
    VectorN xxyx() const { return VectorN   (x(), x(), y(), x()); }

    /** Swizzle operator, returns YXYX. */
    VectorN yxyx() const { return VectorN   (y(), x(), y(), x()); }

    /** Swizzle operator, returns ZXYX. */
    VectorN zxyx() const { return VectorN   (z(), x(), y(), x()); }

    /** Swizzle operator, returns WXYX. */
    VectorN wxyx() const { return VectorN   (w(), x(), y(), x()); }

    /** Swizzle operator, returns XYYX. */
    VectorN xyyx() const { return VectorN   (x(), y(), y(), x()); }

    /** Swizzle operator, returns YYYX. */
    VectorN yyyx() const { return VectorN   (y(), y(), y(), x()); }

    /** Swizzle operator, returns ZYYX. */
    VectorN zyyx() const { return VectorN   (z(), y(), y(), x()); }

    /** Swizzle operator, returns WYYX. */
    VectorN wyyx() const { return VectorN   (w(), y(), y(), x()); }

    /** Swizzle operator, returns XZYX. */
    VectorN xzyx() const { return VectorN   (x(), z(), y(), x()); }

    /** Swizzle operator, returns YZYX. */
    VectorN yzyx() const { return VectorN   (y(), z(), y(), x()); }

    /** Swizzle operator, returns ZZYX. */
    VectorN zzyx() const { return VectorN   (z(), z(), y(), x()); }

    /** Swizzle operator, returns WZYX. */
    VectorN wzyx() const { return VectorN   (w(), z(), y(), x()); }

    /** Swizzle operator, returns XWYX. */
    VectorN xwyx() const { return VectorN   (x(), w(), y(), x()); }

    /** Swizzle operator, returns YWYX. */
    VectorN ywyx() const { return VectorN   (y(), w(), y(), x()); }

    /** Swizzle operator, returns ZWYX. */
    VectorN zwyx() const { return VectorN   (z(), w(), y(), x()); }

    /** Swizzle operator, returns WWYX. */
    VectorN wwyx() const { return VectorN   (w(), w(), y(), x()); }

    /** Swizzle operator, returns XXZX. */
    VectorN xxzx() const { return VectorN   (x(), x(), z(), x()); }

    /** Swizzle operator, returns YXZX. */
    VectorN yxzx() const { return VectorN   (y(), x(), z(), x()); }

    /** Swizzle operator, returns ZXZX. */
    VectorN zxzx() const { return VectorN   (z(), x(), z(), x()); }

    /** Swizzle operator, returns WXZX. */
    VectorN wxzx() const { return VectorN   (w(), x(), z(), x()); }

    /** Swizzle operator, returns XYZX. */
    VectorN xyzx() const { return VectorN   (x(), y(), z(), x()); }

    /** Swizzle operator, returns YYZX. */
    VectorN yyzx() const { return VectorN   (y(), y(), z(), x()); }

    /** Swizzle operator, returns ZYZX. */
    VectorN zyzx() const { return VectorN   (z(), y(), z(), x()); }

    /** Swizzle operator, returns WYZX. */
    VectorN wyzx() const { return VectorN   (w(), y(), z(), x()); }

    /** Swizzle operator, returns XZZX. */
    VectorN xzzx() const { return VectorN   (x(), z(), z(), x()); }

    /** Swizzle operator, returns YZZX. */
    VectorN yzzx() const { return VectorN   (y(), z(), z(), x()); }

    /** Swizzle operator, returns ZZZX. */
    VectorN zzzx() const { return VectorN   (z(), z(), z(), x()); }

    /** Swizzle operator, returns WZZX. */
    VectorN wzzx() const { return VectorN   (w(), z(), z(), x()); }

    /** Swizzle operator, returns XWZX. */
    VectorN xwzx() const { return VectorN   (x(), w(), z(), x()); }

    /** Swizzle operator, returns YWZX. */
    VectorN ywzx() const { return VectorN   (y(), w(), z(), x()); }

    /** Swizzle operator, returns ZWZX. */
    VectorN zwzx() const { return VectorN   (z(), w(), z(), x()); }

    /** Swizzle operator, returns WWZX. */
    VectorN wwzx() const { return VectorN   (w(), w(), z(), x()); }

    /** Swizzle operator, returns XXWX. */
    VectorN xxwx() const { return VectorN   (x(), x(), w(), x()); }

    /** Swizzle operator, returns YXWX. */
    VectorN yxwx() const { return VectorN   (y(), x(), w(), x()); }

    /** Swizzle operator, returns ZXWX. */
    VectorN zxwx() const { return VectorN   (z(), x(), w(), x()); }

    /** Swizzle operator, returns WXWX. */
    VectorN wxwx() const { return VectorN   (w(), x(), w(), x()); }

    /** Swizzle operator, returns XYWX. */
    VectorN xywx() const { return VectorN   (x(), y(), w(), x()); }

    /** Swizzle operator, returns YYWX. */
    VectorN yywx() const { return VectorN   (y(), y(), w(), x()); }

    /** Swizzle operator, returns ZYWX. */
    VectorN zywx() const { return VectorN   (z(), y(), w(), x()); }

    /** Swizzle operator, returns WYWX. */
    VectorN wywx() const { return VectorN   (w(), y(), w(), x()); }

    /** Swizzle operator, returns XZWX. */
    VectorN xzwx() const { return VectorN   (x(), z(), w(), x()); }

    /** Swizzle operator, returns YZWX. */
    VectorN yzwx() const { return VectorN   (y(), z(), w(), x()); }

    /** Swizzle operator, returns ZZWX. */
    VectorN zzwx() const { return VectorN   (z(), z(), w(), x()); }

    /** Swizzle operator, returns WZWX. */
    VectorN wzwx() const { return VectorN   (w(), z(), w(), x()); }

    /** Swizzle operator, returns XWWX. */
    VectorN xwwx() const { return VectorN   (x(), w(), w(), x()); }

    /** Swizzle operator, returns YWWX. */
    VectorN ywwx() const { return VectorN   (y(), w(), w(), x()); }

    /** Swizzle operator, returns ZWWX. */
    VectorN zwwx() const { return VectorN   (z(), w(), w(), x()); }

    /** Swizzle operator, returns WWWX. */
    VectorN wwwx() const { return VectorN   (w(), w(), w(), x()); }

    /** Swizzle operator, returns XXXY. */
    VectorN xxxy() const { return VectorN   (x(), x(), x(), y()); }

    /** Swizzle operator, returns YXXY. */
    VectorN yxxy() const { return VectorN   (y(), x(), x(), y()); }

    /** Swizzle operator, returns ZXXY. */
    VectorN zxxy() const { return VectorN   (z(), x(), x(), y()); }

    /** Swizzle operator, returns WXXY. */
    VectorN wxxy() const { return VectorN   (w(), x(), x(), y()); }

    /** Swizzle operator, returns XYXY. */
    VectorN xyxy() const { return VectorN   (x(), y(), x(), y()); }

    /** Swizzle operator, returns YYXY. */
    VectorN yyxy() const { return VectorN   (y(), y(), x(), y()); }

    /** Swizzle operator, returns ZYXY. */
    VectorN zyxy() const { return VectorN   (z(), y(), x(), y()); }

    /** Swizzle operator, returns WYXY. */
    VectorN wyxy() const { return VectorN   (w(), y(), x(), y()); }

    /** Swizzle operator, returns XZXY. */
    VectorN xzxy() const { return VectorN   (x(), z(), x(), y()); }

    /** Swizzle operator, returns YZXY. */
    VectorN yzxy() const { return VectorN   (y(), z(), x(), y()); }

    /** Swizzle operator, returns ZZXY. */
    VectorN zzxy() const { return VectorN   (z(), z(), x(), y()); }

    /** Swizzle operator, returns WZXY. */
    VectorN wzxy() const { return VectorN   (w(), z(), x(), y()); }

    /** Swizzle operator, returns XWXY. */
    VectorN xwxy() const { return VectorN   (x(), w(), x(), y()); }

    /** Swizzle operator, returns YWXY. */
    VectorN ywxy() const { return VectorN   (y(), w(), x(), y()); }

    /** Swizzle operator, returns ZWXY. */
    VectorN zwxy() const { return VectorN   (z(), w(), x(), y()); }

    /** Swizzle operator, returns WWXY. */
    VectorN wwxy() const { return VectorN   (w(), w(), x(), y()); }

    /** Swizzle operator, returns XXYY. */
    VectorN xxyy() const { return VectorN   (x(), x(), y(), y()); }

    /** Swizzle operator, returns YXYY. */
    VectorN yxyy() const { return VectorN   (y(), x(), y(), y()); }

    /** Swizzle operator, returns ZXYY. */
    VectorN zxyy() const { return VectorN   (z(), x(), y(), y()); }

    /** Swizzle operator, returns WXYY. */
    VectorN wxyy() const { return VectorN   (w(), x(), y(), y()); }

    /** Swizzle operator, returns XYYY. */
    VectorN xyyy() const { return VectorN   (x(), y(), y(), y()); }

    /** Swizzle operator, returns YYYY. */
    VectorN yyyy() const { return VectorN   (y(), y(), y(), y()); }

    /** Swizzle operator, returns ZYYY. */
    VectorN zyyy() const { return VectorN   (z(), y(), y(), y()); }

    /** Swizzle operator, returns WYYY. */
    VectorN wyyy() const { return VectorN   (w(), y(), y(), y()); }

    /** Swizzle operator, returns XZYY. */
    VectorN xzyy() const { return VectorN   (x(), z(), y(), y()); }

    /** Swizzle operator, returns YZYY. */
    VectorN yzyy() const { return VectorN   (y(), z(), y(), y()); }

    /** Swizzle operator, returns ZZYY. */
    VectorN zzyy() const { return VectorN   (z(), z(), y(), y()); }

    /** Swizzle operator, returns WZYY. */
    VectorN wzyy() const { return VectorN   (w(), z(), y(), y()); }

    /** Swizzle operator, returns XWYY. */
    VectorN xwyy() const { return VectorN   (x(), w(), y(), y()); }

    /** Swizzle operator, returns YWYY. */
    VectorN ywyy() const { return VectorN   (y(), w(), y(), y()); }

    /** Swizzle operator, returns ZWYY. */
    VectorN zwyy() const { return VectorN   (z(), w(), y(), y()); }

    /** Swizzle operator, returns WWYY. */
    VectorN wwyy() const { return VectorN   (w(), w(), y(), y()); }

    /** Swizzle operator, returns XXZY. */
    VectorN xxzy() const { return VectorN   (x(), x(), z(), y()); }

    /** Swizzle operator, returns YXZY. */
    VectorN yxzy() const { return VectorN   (y(), x(), z(), y()); }

    /** Swizzle operator, returns ZXZY. */
    VectorN zxzy() const { return VectorN   (z(), x(), z(), y()); }

    /** Swizzle operator, returns WXZY. */
    VectorN wxzy() const { return VectorN   (w(), x(), z(), y()); }

    /** Swizzle operator, returns XYZY. */
    VectorN xyzy() const { return VectorN   (x(), y(), z(), y()); }

    /** Swizzle operator, returns YYZY. */
    VectorN yyzy() const { return VectorN   (y(), y(), z(), y()); }

    /** Swizzle operator, returns ZYZY. */
    VectorN zyzy() const { return VectorN   (z(), y(), z(), y()); }

    /** Swizzle operator, returns WYZY. */
    VectorN wyzy() const { return VectorN   (w(), y(), z(), y()); }

    /** Swizzle operator, returns XZZY. */
    VectorN xzzy() const { return VectorN   (x(), z(), z(), y()); }

    /** Swizzle operator, returns YZZY. */
    VectorN yzzy() const { return VectorN   (y(), z(), z(), y()); }

    /** Swizzle operator, returns ZZZY. */
    VectorN zzzy() const { return VectorN   (z(), z(), z(), y()); }

    /** Swizzle operator, returns WZZY. */
    VectorN wzzy() const { return VectorN   (w(), z(), z(), y()); }

    /** Swizzle operator, returns XWZY. */
    VectorN xwzy() const { return VectorN   (x(), w(), z(), y()); }

    /** Swizzle operator, returns YWZY. */
    VectorN ywzy() const { return VectorN   (y(), w(), z(), y()); }

    /** Swizzle operator, returns ZWZY. */
    VectorN zwzy() const { return VectorN   (z(), w(), z(), y()); }

    /** Swizzle operator, returns WWZY. */
    VectorN wwzy() const { return VectorN   (w(), w(), z(), y()); }

    /** Swizzle operator, returns XXWY. */
    VectorN xxwy() const { return VectorN   (x(), x(), w(), y()); }

    /** Swizzle operator, returns YXWY. */
    VectorN yxwy() const { return VectorN   (y(), x(), w(), y()); }

    /** Swizzle operator, returns ZXWY. */
    VectorN zxwy() const { return VectorN   (z(), x(), w(), y()); }

    /** Swizzle operator, returns WXWY. */
    VectorN wxwy() const { return VectorN   (w(), x(), w(), y()); }

    /** Swizzle operator, returns XYWY. */
    VectorN xywy() const { return VectorN   (x(), y(), w(), y()); }

    /** Swizzle operator, returns YYWY. */
    VectorN yywy() const { return VectorN   (y(), y(), w(), y()); }

    /** Swizzle operator, returns ZYWY. */
    VectorN zywy() const { return VectorN   (z(), y(), w(), y()); }

    /** Swizzle operator, returns WYWY. */
    VectorN wywy() const { return VectorN   (w(), y(), w(), y()); }

    /** Swizzle operator, returns XZWY. */
    VectorN xzwy() const { return VectorN   (x(), z(), w(), y()); }

    /** Swizzle operator, returns YZWY. */
    VectorN yzwy() const { return VectorN   (y(), z(), w(), y()); }

    /** Swizzle operator, returns ZZWY. */
    VectorN zzwy() const { return VectorN   (z(), z(), w(), y()); }

    /** Swizzle operator, returns WZWY. */
    VectorN wzwy() const { return VectorN   (w(), z(), w(), y()); }

    /** Swizzle operator, returns XWWY. */
    VectorN xwwy() const { return VectorN   (x(), w(), w(), y()); }

    /** Swizzle operator, returns YWWY. */
    VectorN ywwy() const { return VectorN   (y(), w(), w(), y()); }

    /** Swizzle operator, returns ZWWY. */
    VectorN zwwy() const { return VectorN   (z(), w(), w(), y()); }

    /** Swizzle operator, returns WWWY. */
    VectorN wwwy() const { return VectorN   (w(), w(), w(), y()); }

    /** Swizzle operator, returns XXXZ. */
    VectorN xxxz() const { return VectorN   (x(), x(), x(), z()); }

    /** Swizzle operator, returns YXXZ. */
    VectorN yxxz() const { return VectorN   (y(), x(), x(), z()); }

    /** Swizzle operator, returns ZXXZ. */
    VectorN zxxz() const { return VectorN   (z(), x(), x(), z()); }

    /** Swizzle operator, returns WXXZ. */
    VectorN wxxz() const { return VectorN   (w(), x(), x(), z()); }

    /** Swizzle operator, returns XYXZ. */
    VectorN xyxz() const { return VectorN   (x(), y(), x(), z()); }

    /** Swizzle operator, returns YYXZ. */
    VectorN yyxz() const { return VectorN   (y(), y(), x(), z()); }

    /** Swizzle operator, returns ZYXZ. */
    VectorN zyxz() const { return VectorN   (z(), y(), x(), z()); }

    /** Swizzle operator, returns WYXZ. */
    VectorN wyxz() const { return VectorN   (w(), y(), x(), z()); }

    /** Swizzle operator, returns XZXZ. */
    VectorN xzxz() const { return VectorN   (x(), z(), x(), z()); }

    /** Swizzle operator, returns YZXZ. */
    VectorN yzxz() const { return VectorN   (y(), z(), x(), z()); }

    /** Swizzle operator, returns ZZXZ. */
    VectorN zzxz() const { return VectorN   (z(), z(), x(), z()); }

    /** Swizzle operator, returns WZXZ. */
    VectorN wzxz() const { return VectorN   (w(), z(), x(), z()); }

    /** Swizzle operator, returns XWXZ. */
    VectorN xwxz() const { return VectorN   (x(), w(), x(), z()); }

    /** Swizzle operator, returns YWXZ. */
    VectorN ywxz() const { return VectorN   (y(), w(), x(), z()); }

    /** Swizzle operator, returns ZWXZ. */
    VectorN zwxz() const { return VectorN   (z(), w(), x(), z()); }

    /** Swizzle operator, returns WWXZ. */
    VectorN wwxz() const { return VectorN   (w(), w(), x(), z()); }

    /** Swizzle operator, returns XXYZ. */
    VectorN xxyz() const { return VectorN   (x(), x(), y(), z()); }

    /** Swizzle operator, returns YXYZ. */
    VectorN yxyz() const { return VectorN   (y(), x(), y(), z()); }

    /** Swizzle operator, returns ZXYZ. */
    VectorN zxyz() const { return VectorN   (z(), x(), y(), z()); }

    /** Swizzle operator, returns WXYZ. */
    VectorN wxyz() const { return VectorN   (w(), x(), y(), z()); }

    /** Swizzle operator, returns XYYZ. */
    VectorN xyyz() const { return VectorN   (x(), y(), y(), z()); }

    /** Swizzle operator, returns YYYZ. */
    VectorN yyyz() const { return VectorN   (y(), y(), y(), z()); }

    /** Swizzle operator, returns ZYYZ. */
    VectorN zyyz() const { return VectorN   (z(), y(), y(), z()); }

    /** Swizzle operator, returns WYYZ. */
    VectorN wyyz() const { return VectorN   (w(), y(), y(), z()); }

    /** Swizzle operator, returns XZYZ. */
    VectorN xzyz() const { return VectorN   (x(), z(), y(), z()); }

    /** Swizzle operator, returns YZYZ. */
    VectorN yzyz() const { return VectorN   (y(), z(), y(), z()); }

    /** Swizzle operator, returns ZZYZ. */
    VectorN zzyz() const { return VectorN   (z(), z(), y(), z()); }

    /** Swizzle operator, returns WZYZ. */
    VectorN wzyz() const { return VectorN   (w(), z(), y(), z()); }

    /** Swizzle operator, returns XWYZ. */
    VectorN xwyz() const { return VectorN   (x(), w(), y(), z()); }

    /** Swizzle operator, returns YWYZ. */
    VectorN ywyz() const { return VectorN   (y(), w(), y(), z()); }

    /** Swizzle operator, returns ZWYZ. */
    VectorN zwyz() const { return VectorN   (z(), w(), y(), z()); }

    /** Swizzle operator, returns WWYZ. */
    VectorN wwyz() const { return VectorN   (w(), w(), y(), z()); }

    /** Swizzle operator, returns XXZZ. */
    VectorN xxzz() const { return VectorN   (x(), x(), z(), z()); }

    /** Swizzle operator, returns YXZZ. */
    VectorN yxzz() const { return VectorN   (y(), x(), z(), z()); }

    /** Swizzle operator, returns ZXZZ. */
    VectorN zxzz() const { return VectorN   (z(), x(), z(), z()); }

    /** Swizzle operator, returns WXZZ. */
    VectorN wxzz() const { return VectorN   (w(), x(), z(), z()); }

    /** Swizzle operator, returns XYZZ. */
    VectorN xyzz() const { return VectorN   (x(), y(), z(), z()); }

    /** Swizzle operator, returns YYZZ. */
    VectorN yyzz() const { return VectorN   (y(), y(), z(), z()); }

    /** Swizzle operator, returns ZYZZ. */
    VectorN zyzz() const { return VectorN   (z(), y(), z(), z()); }

    /** Swizzle operator, returns WYZZ. */
    VectorN wyzz() const { return VectorN   (w(), y(), z(), z()); }

    /** Swizzle operator, returns XZZZ. */
    VectorN xzzz() const { return VectorN   (x(), z(), z(), z()); }

    /** Swizzle operator, returns YZZZ. */
    VectorN yzzz() const { return VectorN   (y(), z(), z(), z()); }

    /** Swizzle operator, returns ZZZZ. */
    VectorN zzzz() const { return VectorN   (z(), z(), z(), z()); }

    /** Swizzle operator, returns WZZZ. */
    VectorN wzzz() const { return VectorN   (w(), z(), z(), z()); }

    /** Swizzle operator, returns XWZZ. */
    VectorN xwzz() const { return VectorN   (x(), w(), z(), z()); }

    /** Swizzle operator, returns YWZZ. */
    VectorN ywzz() const { return VectorN   (y(), w(), z(), z()); }

    /** Swizzle operator, returns ZWZZ. */
    VectorN zwzz() const { return VectorN   (z(), w(), z(), z()); }

    /** Swizzle operator, returns WWZZ. */
    VectorN wwzz() const { return VectorN   (w(), w(), z(), z()); }

    /** Swizzle operator, returns XXWZ. */
    VectorN xxwz() const { return VectorN   (x(), x(), w(), z()); }

    /** Swizzle operator, returns YXWZ. */
    VectorN yxwz() const { return VectorN   (y(), x(), w(), z()); }

    /** Swizzle operator, returns ZXWZ. */
    VectorN zxwz() const { return VectorN   (z(), x(), w(), z()); }

    /** Swizzle operator, returns WXWZ. */
    VectorN wxwz() const { return VectorN   (w(), x(), w(), z()); }

    /** Swizzle operator, returns XYWZ. */
    VectorN xywz() const { return VectorN   (x(), y(), w(), z()); }

    /** Swizzle operator, returns YYWZ. */
    VectorN yywz() const { return VectorN   (y(), y(), w(), z()); }

    /** Swizzle operator, returns ZYWZ. */
    VectorN zywz() const { return VectorN   (z(), y(), w(), z()); }

    /** Swizzle operator, returns WYWZ. */
    VectorN wywz() const { return VectorN   (w(), y(), w(), z()); }

    /** Swizzle operator, returns XZWZ. */
    VectorN xzwz() const { return VectorN   (x(), z(), w(), z()); }

    /** Swizzle operator, returns YZWZ. */
    VectorN yzwz() const { return VectorN   (y(), z(), w(), z()); }

    /** Swizzle operator, returns ZZWZ. */
    VectorN zzwz() const { return VectorN   (z(), z(), w(), z()); }

    /** Swizzle operator, returns WZWZ. */
    VectorN wzwz() const { return VectorN   (w(), z(), w(), z()); }

    /** Swizzle operator, returns XWWZ. */
    VectorN xwwz() const { return VectorN   (x(), w(), w(), z()); }

    /** Swizzle operator, returns YWWZ. */
    VectorN ywwz() const { return VectorN   (y(), w(), w(), z()); }

    /** Swizzle operator, returns ZWWZ. */
    VectorN zwwz() const { return VectorN   (z(), w(), w(), z()); }

    /** Swizzle operator, returns WWWZ. */
    VectorN wwwz() const { return VectorN   (w(), w(), w(), z()); }

    /** Swizzle operator, returns XXXW. */
    VectorN xxxw() const { return VectorN   (x(), x(), x(), w()); }

    /** Swizzle operator, returns YXXW. */
    VectorN yxxw() const { return VectorN   (y(), x(), x(), w()); }

    /** Swizzle operator, returns ZXXW. */
    VectorN zxxw() const { return VectorN   (z(), x(), x(), w()); }

    /** Swizzle operator, returns WXXW. */
    VectorN wxxw() const { return VectorN   (w(), x(), x(), w()); }

    /** Swizzle operator, returns XYXW. */
    VectorN xyxw() const { return VectorN   (x(), y(), x(), w()); }

    /** Swizzle operator, returns YYXW. */
    VectorN yyxw() const { return VectorN   (y(), y(), x(), w()); }

    /** Swizzle operator, returns ZYXW. */
    VectorN zyxw() const { return VectorN   (z(), y(), x(), w()); }

    /** Swizzle operator, returns WYXW. */
    VectorN wyxw() const { return VectorN   (w(), y(), x(), w()); }

    /** Swizzle operator, returns XZXW. */
    VectorN xzxw() const { return VectorN   (x(), z(), x(), w()); }

    /** Swizzle operator, returns YZXW. */
    VectorN yzxw() const { return VectorN   (y(), z(), x(), w()); }

    /** Swizzle operator, returns ZZXW. */
    VectorN zzxw() const { return VectorN   (z(), z(), x(), w()); }

    /** Swizzle operator, returns WZXW. */
    VectorN wzxw() const { return VectorN   (w(), z(), x(), w()); }

    /** Swizzle operator, returns XWXW. */
    VectorN xwxw() const { return VectorN   (x(), w(), x(), w()); }

    /** Swizzle operator, returns YWXW. */
    VectorN ywxw() const { return VectorN   (y(), w(), x(), w()); }

    /** Swizzle operator, returns ZWXW. */
    VectorN zwxw() const { return VectorN   (z(), w(), x(), w()); }

    /** Swizzle operator, returns WWXW. */
    VectorN wwxw() const { return VectorN   (w(), w(), x(), w()); }

    /** Swizzle operator, returns XXYW. */
    VectorN xxyw() const { return VectorN   (x(), x(), y(), w()); }

    /** Swizzle operator, returns YXYW. */
    VectorN yxyw() const { return VectorN   (y(), x(), y(), w()); }

    /** Swizzle operator, returns ZXYW. */
    VectorN zxyw() const { return VectorN   (z(), x(), y(), w()); }

    /** Swizzle operator, returns WXYW. */
    VectorN wxyw() const { return VectorN   (w(), x(), y(), w()); }

    /** Swizzle operator, returns XYYW. */
    VectorN xyyw() const { return VectorN   (x(), y(), y(), w()); }

    /** Swizzle operator, returns YYYW. */
    VectorN yyyw() const { return VectorN   (y(), y(), y(), w()); }

    /** Swizzle operator, returns ZYYW. */
    VectorN zyyw() const { return VectorN   (z(), y(), y(), w()); }

    /** Swizzle operator, returns WYYW. */
    VectorN wyyw() const { return VectorN   (w(), y(), y(), w()); }

    /** Swizzle operator, returns XZYW. */
    VectorN xzyw() const { return VectorN   (x(), z(), y(), w()); }

    /** Swizzle operator, returns YZYW. */
    VectorN yzyw() const { return VectorN   (y(), z(), y(), w()); }

    /** Swizzle operator, returns ZZYW. */
    VectorN zzyw() const { return VectorN   (z(), z(), y(), w()); }

    /** Swizzle operator, returns WZYW. */
    VectorN wzyw() const { return VectorN   (w(), z(), y(), w()); }

    /** Swizzle operator, returns XWYW. */
    VectorN xwyw() const { return VectorN   (x(), w(), y(), w()); }

    /** Swizzle operator, returns YWYW. */
    VectorN ywyw() const { return VectorN   (y(), w(), y(), w()); }

    /** Swizzle operator, returns ZWYW. */
    VectorN zwyw() const { return VectorN   (z(), w(), y(), w()); }

    /** Swizzle operator, returns WWYW. */
    VectorN wwyw() const { return VectorN   (w(), w(), y(), w()); }

    /** Swizzle operator, returns XXZW. */
    VectorN xxzw() const { return VectorN   (x(), x(), z(), w()); }

    /** Swizzle operator, returns YXZW. */
    VectorN yxzw() const { return VectorN   (y(), x(), z(), w()); }

    /** Swizzle operator, returns ZXZW. */
    VectorN zxzw() const { return VectorN   (z(), x(), z(), w()); }

    /** Swizzle operator, returns WXZW. */
    VectorN wxzw() const { return VectorN   (w(), x(), z(), w()); }

    /** Swizzle operator, returns XYZW. */
    VectorN xyzw() const { return VectorN   (x(), y(), z(), w()); }

    /** Swizzle operator, returns YYZW. */
    VectorN yyzw() const { return VectorN   (y(), y(), z(), w()); }

    /** Swizzle operator, returns ZYZW. */
    VectorN zyzw() const { return VectorN   (z(), y(), z(), w()); }

    /** Swizzle operator, returns WYZW. */
    VectorN wyzw() const { return VectorN   (w(), y(), z(), w()); }

    /** Swizzle operator, returns XZZW. */
    VectorN xzzw() const { return VectorN   (x(), z(), z(), w()); }

    /** Swizzle operator, returns YZZW. */
    VectorN yzzw() const { return VectorN   (y(), z(), z(), w()); }

    /** Swizzle operator, returns ZZZW. */
    VectorN zzzw() const { return VectorN   (z(), z(), z(), w()); }

    /** Swizzle operator, returns WZZW. */
    VectorN wzzw() const { return VectorN   (w(), z(), z(), w()); }

    /** Swizzle operator, returns XWZW. */
    VectorN xwzw() const { return VectorN   (x(), w(), z(), w()); }

    /** Swizzle operator, returns YWZW. */
    VectorN ywzw() const { return VectorN   (y(), w(), z(), w()); }

    /** Swizzle operator, returns ZWZW. */
    VectorN zwzw() const { return VectorN   (z(), w(), z(), w()); }

    /** Swizzle operator, returns WWZW. */
    VectorN wwzw() const { return VectorN   (w(), w(), z(), w()); }

    /** Swizzle operator, returns XXWW. */
    VectorN xxww() const { return VectorN   (x(), x(), w(), w()); }

    /** Swizzle operator, returns YXWW. */
    VectorN yxww() const { return VectorN   (y(), x(), w(), w()); }

    /** Swizzle operator, returns ZXWW. */
    VectorN zxww() const { return VectorN   (z(), x(), w(), w()); }

    /** Swizzle operator, returns WXWW. */
    VectorN wxww() const { return VectorN   (w(), x(), w(), w()); }

    /** Swizzle operator, returns XYWW. */
    VectorN xyww() const { return VectorN   (x(), y(), w(), w()); }

    /** Swizzle operator, returns YYWW. */
    VectorN yyww() const { return VectorN   (y(), y(), w(), w()); }

    /** Swizzle operator, returns ZYWW. */
    VectorN zyww() const { return VectorN   (z(), y(), w(), w()); }

    /** Swizzle operator, returns WYWW. */
    VectorN wyww() const { return VectorN   (w(), y(), w(), w()); }

    /** Swizzle operator, returns XZWW. */
    VectorN xzww() const { return VectorN   (x(), z(), w(), w()); }

    /** Swizzle operator, returns YZWW. */
    VectorN yzww() const { return VectorN   (y(), z(), w(), w()); }

    /** Swizzle operator, returns ZZWW. */
    VectorN zzww() const { return VectorN   (z(), z(), w(), w()); }

    /** Swizzle operator, returns WZWW. */
    VectorN wzww() const { return VectorN   (w(), z(), w(), w()); }

    /** Swizzle operator, returns XWWW. */
    VectorN xwww() const { return VectorN   (x(), w(), w(), w()); }

    /** Swizzle operator, returns YWWW. */
    VectorN ywww() const { return VectorN   (y(), w(), w(), w()); }

    /** Swizzle operator, returns ZWWW. */
    VectorN zwww() const { return VectorN   (z(), w(), w(), w()); }

    /** Swizzle operator, returns WWWW. */
    VectorN wwww() const { return VectorN   (w(), w(), w(), w()); }

}; // class VectorN<4, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API VectorN<4, Real>;
#endif

/** The default 4D real vector class. */
typedef VectorN<4, Real> Vector4;

} // namespace Thea

#endif
