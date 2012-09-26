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

#ifndef __Thea_CoordinateFrameN_hpp__
#define __Thea_CoordinateFrameN_hpp__

#include "Common.hpp"
#include "RigidTransformN.hpp"

namespace Thea {

// Forward declarations
template <long N, typename T> class CoordinateFrameN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for a coordinate frame in N-space, defined by N orthonormal vectors.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ CoordinateFrameNBase : public RigidTransformN<N, T>
{
  public:
    typedef CoordinateFrameN<N, T>             CoordinateFrameT;  ///< N-dimensional coordinate frame type.
    typedef RigidTransformN<N, T>              RigidTransformT;   ///< N-dimensional rigid transform.
    typedef AffineTransformN<N, T>             AffineTransformT;  ///< N-dimensional affine transform.
    typedef typename RigidTransformT::VectorT  VectorT;
    typedef typename RigidTransformT::MatrixT  MatrixT;

    THEA_DEF_POINTER_TYPES(CoordinateFrameT, shared_ptr, weak_ptr)

    /** Default constructor. Constructs the identity frame. */
    CoordinateFrameNBase() {}

    /** Construct from a rigid transform. */
    CoordinateFrameNBase(RigidTransformT const & src) : RigidTransformT(src) {}

    /**
     * Construct from an affine transform, assuming it is rigid (<b>use with caution</b> since it can break the orthonormality
     * guarantee).
     */
    static CoordinateFrameT _fromAffine(AffineTransformT const & aff_)
    {
      return CoordinateFrameT(RigidTransformT::_fromAffine(aff_));
    }

    /** Get an axis of the frame. */
    VectorT getAxis(long i) const { return this->getRotation().getColumn(i); }

    /** Get the inverse transform. */
    CoordinateFrameT inverse() const
    {
      return CoordinateFrameT(RigidTransformT::inverse());
    }

    /** Compose this frame, treated as a rigid transform, with another. The other is applied first. */
    CoordinateFrameT operator*(CoordinateFrameT const & rhs) const
    {
      return CoordinateFrameT(this->RigidTransformT::operator*(rhs));
    }

    /** Transform a point from the local space of the coordinate frame to world space. */
    VectorT operator*(VectorT const & p) const
    {
      return this->RigidTransformT::operator*(p);
    }

    /** Transform a point from the local space of the coordinate frame to world space. */
    VectorT pointToWorldSpace(VectorT const & p) const
    {
      return this->RigidTransformT::operator*(p);
    }

    /** Transform a point from world space to the local space of the coordinate frame. */
    VectorT pointToObjectSpace(VectorT const & p) const
    {
      return (p - this->getTranslation()) * this->getRotation();
    }

    /** Transform a direction vector from the local space of the coordinate frame to world space. */
    VectorT vectorToWorldSpace(VectorT const & v) const
    {
      return this->getRotation() * v;
    }

    /** Transform a direction vector from world space to the local space of the coordinate frame. */
    VectorT vectorToObjectSpace(VectorT const & v) const
    {
      return v * this->getRotation();
    }

    /** Transform a normal from the local space of the coordinate frame to world space. */
    VectorT normalToWorldSpace(VectorT const & n) const
    {
      return this->getRotation() * n;
    }

    /** Transform a normal from world space to the local space of the coordinate frame. */
    VectorT normalToObjectSpace(VectorT const & n) const
    {
      return n * this->getRotation();
    }

    /** Get the identity frame (same as the world frame). */
    static CoordinateFrameT const & identity()
    {
      static CoordinateFrameT const idty(RigidTransformT::identity());
      return idty;
    }

  private:
    // Hide these functions from the default interface
    RigidTransformT operator*(RigidTransformT const & rhs) const      { return RigidTransformT(); }
    static RigidTransformT translation(VectorT const & translation_)  { return RigidTransformT(); }

}; // class CoordinateFrameNBase

} // namespace Internal

/** A coordinate frame in N-space, defined by N orthonormal vectors. */
template <long N, typename T>
class /* THEA_API */ CoordinateFrameN : public Internal::CoordinateFrameNBase<N, T>
{
  private:
    typedef Internal::CoordinateFrameNBase<N, T> BaseT;

  public:
    typedef typename BaseT::RigidTransformT  RigidTransformT;
    typedef typename BaseT::VectorT          VectorT;
    typedef typename BaseT::MatrixT          MatrixT;

    /** Default constructor. Constructs the identity frame. */
    CoordinateFrameN() {}

    /** Construct from a rigid transform. */
    CoordinateFrameN(RigidTransformT const & src) : BaseT(src) {}

}; // class CoordinateFrameN

/** Pipe a textual representation of a coordinate frame to a <code>std::ostream</code>. */
template <long N, typename T>
std::ostream &
operator<<(std::ostream & os, CoordinateFrameN<N, T> const & cf)
{
  return os << cf.toString();
}

} // namespace Thea

#include "CoordinateFrame2.hpp"
#include "CoordinateFrame3.hpp"

#endif
