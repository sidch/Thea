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

#ifndef __Thea_AffineTransform2_hpp__
#define __Thea_AffineTransform2_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"

namespace Thea {

/** An affine transformation in 2-space, consisting of a linear transformation (2x2 matrix) plus a translation. */
template <typename T>
class /* THEA_API */ AffineTransformN<2, T> : public Internal::AffineTransformNBase<2, T>
{
  private:
    typedef Internal::AffineTransformNBase<2, T> BaseT;

  public:
    typedef typename BaseT::VectorT  VectorT;
    typedef typename BaseT::MatrixT  MatrixT;

    /** Default constructor. Does not initialize anything. */
    AffineTransformN() {}

    /** Construct from a linear transform, followed by a translation. */
    AffineTransformN(MatrixT const & linear_, VectorT const & translation_ = VectorT::zero()) : BaseT(linear_, translation_) {}

    /** Copy constructor. */
    AffineTransformN(AffineTransformN const & src) : BaseT(src) {}

    /**
     * Construct from 2 basis vectors, specifying the linear transform, and a translation. The three arguments form the columns
     * of the 2x3 matrix specifying the transform.
     */
    AffineTransformN(VectorT const & x, VectorT const & y, VectorT const & translation_)
    : BaseT(MatrixT(x[0], y[0], x[1], y[1]), translation_) {}

    /** Construct from a 2x3 array. */
    AffineTransformN(T const &  m00, T const &  m01, T const &  m02,
                     T const &  m10, T const &  m11, T const &  m12)
    : BaseT(MatrixT(m00, m01, m10, m11), VectorT(m02, m12)) {}

    using BaseT::scaling;

    /** Construct a scaling transform. */
    static AffineTransformN scaling(T const & sx, T const & sy)
    {
      return BaseT::scaling(VectorT(sx, sy));
    }

    using BaseT::translation;

    /** Construct a translation. */
    static AffineTransformN translation(T const & tx, T const & ty)
    {
      return BaseT::translation(VectorT(tx, ty));
    }

    /** Construct a rotation specified by an angle (in radians) around the origin. */
    static AffineTransformN rotation(Real radians)
    {
      return AffineTransformN(MatrixT::rotation(radians), VectorT::zero());
    }

}; // class AffineTransformN<2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AffineTransformN<2, Real>;
#endif

/** The default affine transform class in real 2-space. */
typedef AffineTransformN<2, Real> AffineTransform2;

} // namespace Thea

#endif
