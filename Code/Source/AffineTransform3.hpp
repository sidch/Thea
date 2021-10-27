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

#ifndef __Thea_AffineTransform3_hpp__
#define __Thea_AffineTransform3_hpp__

#include "Common.hpp"
#include "AffineTransformN.hpp"

namespace Thea {

/** An affine transformation in 3-space, consisting of a linear transformation (3x3 matrix) plus a translation. */
template <typename T>
class /* THEA_API */ AffineTransformN<3, T> : public Internal::AffineTransformNBase<3, T>
{
  private:
    typedef Internal::AffineTransformNBase<3, T> BaseT;

  public:
    typedef typename BaseT::VectorT  VectorT;
    typedef typename BaseT::MatrixT  MatrixT;

    /** Default constructor. Does not initialize anything. */
    AffineTransformN() {}

    /** Construct from a linear transform, followed by a translation. */
    AffineTransformN(MatrixT const & linear_, VectorT const & translation_ = VectorT::Zero()) : BaseT(linear_, translation_) {}

    /** Copy constructor. */
    AffineTransformN(AffineTransformN const & src) : BaseT(src) {}

    /**
     * Construct from 3 basis vectors, specifying the linear transform, and a translation. The four arguments form the columns
     * of the 3x4 matrix specifying the transform.
     */
    AffineTransformN(VectorT const & x, VectorT const & y, VectorT const & z, VectorT const & translation_)
    : BaseT((MatrixT() << x, y, z).finished(), translation_) {}

    /** Construct from a 3x4 array. */
    AffineTransformN(T const &  m00, T const &  m01, T const &  m02, T const & m03,
                     T const &  m10, T const &  m11, T const &  m12, T const & m13,
                     T const &  m20, T const &  m21, T const &  m22, T const & m23)
    : BaseT((MatrixT() << m00, m01, m02, m10, m11, m12, m20, m21, m22).finished(), VectorT(m03, m13, m23)) {}

    using BaseT::scaling;

    /** Construct a scaling transform. */
    static AffineTransformN scaling(T const & sx, T const & sy, T const & sz)
    {
      return BaseT::scaling(VectorT(sx, sy, sz));
    }

    using BaseT::translation;

    /** Construct a translation. */
    static AffineTransformN translation(T const & tx, T const & ty, T const & tz)
    {
      return BaseT::translation(VectorT(tx, ty, tz));
    }

    /** Construct a rotation specified by an angle (in radians) around an axis. */
    static AffineTransformN rotationAxisAngle(VectorT const & axis, Real radians)
    {
      return AffineTransformN(Math::rotationAxisAngle<T>(axis, radians), VectorT::Zero());
    }

    /**
     * Rotate about the Z axis by the roll angle, then the Y axis by the pitch angle, and finally the X axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesXYZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesXYZ<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Rotate about the Y axis by the roll angle, then the Z axis by the pitch angle, and finally the X axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesXZY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesXZY<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Rotate about the Z axis by the roll angle, then the X axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesYXZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesYXZ<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Rotate about the X axis by the roll angle, then the Z axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesYZX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesYZX<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Rotate about the Y axis by the roll angle, then the X axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesZXY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesZXY<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Rotate about the X axis by the roll angle, then the Y axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static AffineTransformN rotationEulerAnglesZYX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      return AffineTransformN(Math::rotationEulerAnglesZYX<T>(yaw_radians, pitch_radians, roll_radians), VectorT::Zero());
    }

    /**
     * Return the transformation corresponding to the rotation from one direction vector to another.
     *
     * @param start_dir The vector to rotate from.
     * @param end_dir The vector to rotate to.
     * @param normalize_dirs If false, the directions will be assumed to have been pre-normalized to unit length before being
     *   passed to this function.
     */
    static AffineTransformN rotationArc(VectorT const & start_dir, VectorT const & end_dir, bool normalize_dirs = true)
    {
      return AffineTransformN(Math::rotationArc<T>(start_dir, end_dir, normalize_dirs), VectorT::Zero());
    }

}; // class AffineTransformN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AffineTransformN<3, Real>;
#endif

/** The default affine transform class in real 3-space. */
typedef AffineTransformN<3, Real> AffineTransform3;

} // namespace Thea

#endif
