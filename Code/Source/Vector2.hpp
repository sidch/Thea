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

#ifndef __Thea_Vector2_hpp__
#define __Thea_Vector2_hpp__

#include "VectorN.hpp"

namespace Thea {

/** 2-dimensional vectors on a field T. */
template <typename T>
class /* THEA_API */ VectorN<2, T> : public Internal::VectorNBase<2, T>
{
  private:
    typedef Internal::VectorNBase<2, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    VectorN() {}

    /** Initialize all components to a single value. */
    explicit VectorN(T const & fill_value) : BaseT(fill_value) {}

    /** Copy constructor. */
    template <typename U> VectorN(VectorN<2, U> const & src) : BaseT(src) {}

    /** Initialize from a column matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<2, 1, U> const & src) : BaseT(src) {}

    /** Initialize from a row matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<1, 2, U> const & src) : BaseT(src) {}

    /** Initialize all components of the vector. */
    VectorN(T const & x_, T const & y_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
    }

    /** Set all elements of the vector. */
    void set(T const & x_, T const & y_)
    {
      (*this)[0] = x_;
      (*this)[1] = y_;
    }

    /**
     * Get a unit vector perpendicular to this one, forming an orthonormal right-handed basis (u, this->unit()). In other words,
     * if this is the Y axis of the local frame, then the function returns the X axis.
     *
     * In 2D, the behavior of this function is identical to createOrthonormalBasis().
     */
    VectorN getOrthogonalDirection() const { return VectorN(y(), -x()); }

    /**
     * Get a unit vector perpendicular to this one, forming an orthonormal right-handed basis (u, this->unit()). In other words,
     * if this is the Y axis of the local frame, then the function returns the X axis.
     *
     * In 2D, the behavior of this function is identical to getOrthogonalDirection().
     */
    void createOrthonormalBasis(VectorN & u) const { u = getOrthogonalDirection(); }

    /** Get a unit vector along positive X. */
    static VectorN const & unitX() { static VectorN const ux(1, 0); return ux; }

    /** Get a unit vector along positive Y. */
    static VectorN const & unitY() { static VectorN const uy(0, 1); return uy; }

    /** Get the X coordinate. */
    T const & x() const { return (*this)[0]; }

    /** Get the X coordinate. */
    T & x() { return (*this)[0]; }

    /** Get the Y coordinate. */
    T const & y() const { return (*this)[1]; }

    /** Get the Y coordinate. */
    T & y() { return (*this)[1]; }

    /** Swizzle operator, returns XX. */
    VectorN xx() const { return VectorN   (x(), x()); }

    /** Swizzle operator, returns YX. */
    VectorN yx() const { return VectorN   (y(), x()); }

    /** Swizzle operator, returns XY. */
    VectorN xy() const { return VectorN   (x(), y()); }

    /** Swizzle operator, returns YY. */
    VectorN yy() const { return VectorN   (y(), y()); }

}; // class VectorN<2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API VectorN<2, Real>;
#endif

/** The default 2D real vector class. */
typedef VectorN<2, Real> Vector2;

} // namespace Thea

#endif
