//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
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
      VectorT f = (center - eye).normalized();
      this->_setRotation((MatrixT() << f[0], -f[1], f[1], f[0]).finished());
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
