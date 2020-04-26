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

#ifndef __Thea_RigidTransform3_hpp__
#define __Thea_RigidTransform3_hpp__

#include "Common.hpp"
#include "AffineTransform3.hpp"
#include "RigidTransformN.hpp"

namespace Thea {

/**
 * A rigid transformation in 3-space, consisting of a rotation followed by a translation.
 *
 * @note While this is technically an affine transform, it restricts enough functionality to make a separate implementation
 * preferable. It can be trivially converted to an AffineTransform3 using the implicit conversion operator or toAffine().
 */
template <typename T>
class /* THEA_API */ RigidTransformN<3, T> : public Internal::RigidTransformNBase<3, T>
{
  private:
    typedef Internal::RigidTransformNBase<3, T> BaseT;

  public:
    typedef typename BaseT::AffineTransformT  AffineTransformT;
    typedef typename BaseT::VectorT           VectorT;
    typedef typename BaseT::MatrixT           MatrixT;

    /** Default constructor. Constructs the identity transform. */
    RigidTransformN() {}

    using BaseT::translation;

    /** Construct a translation. */
    static RigidTransformN translation(T const & tx, T const & ty, T const & tz)
    {
      return BaseT::translation(VectorT(tx, ty, tz));
    }

    /** Construct a rotation specified by an angle (in radians) around an axis. */
    static RigidTransformN rotationAxisAngle(VectorT const & axis, Real radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationAxisAngle(axis, radians));
    }

    /**
     * Rotate about the Z axis by the roll angle, then the Y axis by the pitch angle, and finally the X axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesXYZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesXYZ(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Rotate about the Y axis by the roll angle, then the Z axis by the pitch angle, and finally the X axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesXZY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesXZY(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Rotate about the Z axis by the roll angle, then the X axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesYXZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesYXZ(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Rotate about the X axis by the roll angle, then the Z axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesYZX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesYZX(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Rotate about the Y axis by the roll angle, then the X axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesZXY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesZXY(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Rotate about the X axis by the roll angle, then the Y axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static RigidTransformN rotationEulerAnglesZYX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationEulerAnglesZYX(yaw_radians, pitch_radians, roll_radians));
    }

    /**
     * Return the transformation corresponding to the rotation from one direction vector to another.
     *
     * @param start_dir The vector to rotate from.
     * @param end_dir The vector to rotate to.
     * @param normalize_dirs If false, the directions will be assumed to have been pre-normalized to unit length before being
     *   passed to this function.
     */
    static RigidTransformN rotationArc(VectorT const & start_dir, VectorT const & end_dir, bool normalize_dirs = true)
    {
      return BaseT::_fromAffine(AffineTransformT::rotationArc(start_dir, end_dir, normalize_dirs));
    }

}; // class RigidTransformN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API RigidTransformN<3, Real>;
#endif

/** The default affine transform class in real 3-space. */
typedef RigidTransformN<3, Real> RigidTransform3;

} // namespace Thea

#endif
