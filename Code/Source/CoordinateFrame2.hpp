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

#ifndef __Thea_CoordinateFrame2_hpp__
#define __Thea_CoordinateFrame2_hpp__

#include "Common.hpp"
#include "CoordinateFrameN.hpp"

namespace Thea {

/** A coordinate frame in 2-space, defined by two orthonormal vectors. */
template <typename T>
class /* THEA_API */ CoordinateFrameN<2, T> : public Internal::CoordinateFrameNBase<2, T>
{
  private:
    typedef Internal::CoordinateFrameNBase<2, T> BaseT;

  public:
    typedef typename BaseT::RigidTransformT  RigidTransformT;
    typedef typename BaseT::VectorT          VectorT;
    typedef typename BaseT::MatrixT          MatrixT;

    /** Default constructor. Constructs the identity frame. */
    CoordinateFrameN() {}

    /** Construct from a rigid transform. */
    CoordinateFrameN(RigidTransformT const & src) : BaseT(src) {}

    /** Construct from a viewing position (\a eye) and a look-at position (\a center). */
    CoordinateFrameN(VectorT const & eye, VectorT const & center)
    {
      set(eye, center);
    }

    /** Initialize from a viewing position (\a eye) and a look-at position (\a center). */
    void set(VectorT const & eye, VectorT const & center)
    {
      VectorT f = (center - eye).unit();
      this->_setRotation(MatrixT(f[0], -f[1], f[1], f[0]));
      this->setTranslation(eye);
    }

  private:
    // Hide these functions from the default interface
    static RigidTransformT rotation(Real radians) { return RigidTransformT(); }

}; // class CoordinateFrameN<2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API CoordinateFrameN<2, Real>;
#endif

/** The default coordinate frame class in real 2-space. */
typedef CoordinateFrameN<2, Real> CoordinateFrame2;

} // namespace Thea

#endif
