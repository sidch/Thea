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

#ifndef __Thea_CoordinateFrame3_hpp__
#define __Thea_CoordinateFrame3_hpp__

#include "Common.hpp"
#include "CoordinateFrameN.hpp"

namespace Thea {

/** A coordinate frame in 3-space, defined by three orthonormal vectors. */
template <typename T>
class /* THEA_API */ CoordinateFrameN<3, T> : public Internal::CoordinateFrameNBase<3, T>
{
  private:
    typedef Internal::CoordinateFrameNBase<3, T> BaseT;

  public:
    typedef typename BaseT::RigidTransformT  RigidTransformT;
    typedef typename BaseT::VectorT          VectorT;
    typedef typename BaseT::MatrixT          MatrixT;

    /** Default constructor. Constructs the identity frame. */
    CoordinateFrameN() {}

    /** Construct from a rigid transform. */
    CoordinateFrameN(RigidTransformT const & src) : BaseT(src) {}

    /**
     * Construct from a viewing position (\a eye), look-at position (\a look_at) and up direction (\a up). The up direction need
     * not be exactly orthogonal to \a look_at - \a eye or be a unit vector.
     */
    static CoordinateFrameN fromViewFrame(VectorT const & eye, VectorT const & look_at, VectorT const & up)
    {
      CoordinateFrameN tr;
      tr.setViewFrame(eye, look_at, up);
      return tr;
    }

    /**
     * Initialize from a viewing position (\a eye), look-at position (\a look_at) and up direction (\a up). The up direction
     * need not be exactly orthogonal to \a look_at - \a eye or be a unit vector.
     */
    void setViewFrame(VectorT const & eye, VectorT const & look_at, VectorT const & up)
    {
      // See documentation of gluLookAt
      VectorT f = (look_at - eye).normalized();
      VectorT s = f.cross(up).normalized();
      VectorT u = s.cross(f);

      this->_setRotation((MatrixT() << s[0], u[0], -f[0],
                                       s[1], u[1], -f[1],
                                       s[2], u[2], -f[2]).finished());

      this->setTranslation(eye);
    }

    /** Get the viewing direction (the negative Z axis of the frame). */
    VectorT lookVector() const
    {
      return -this->getAxis(2);
    }

    /** Get the up direction (the Y axis of the frame). */
    VectorT upVector() const
    {
      return this->getAxis(1);
    }

    /** Get the right-hand direction (the X axis of the frame). */
    VectorT rightVector() const
    {
      return this->getAxis(0);
    }

  private:
    // Hide these functions from the default interface
    static RigidTransformT rotationAxisAngle(VectorT const & axis, Real radians)    { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesXYZ(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesXZY(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesYXZ(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesYZX(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesZXY(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }
    static RigidTransformT rotationEulerAnglesZYX(Real yaw, Real pitch, Real roll)  { return RigidTransformT(); }

    static RigidTransformT rotationArc(VectorT const & start_dir, VectorT const & end_dir, bool normalize_dirs = true)
    { return RigidTransformT(); }

}; // class CoordinateFrameN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API CoordinateFrameN<3, Real>;
#endif

/** The default coordinate frame class in real 3-space. */
typedef CoordinateFrameN<3, Real> CoordinateFrame3;

} // namespace Thea

#endif
