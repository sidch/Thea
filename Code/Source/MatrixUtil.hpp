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
// First version: 2011
//
//============================================================================

#ifndef __Thea_MatrixUtil_hpp__
#define __Thea_MatrixUtil_hpp__

#include "Common.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include <algorithm>
#include <type_traits>

namespace Thea {
namespace Math {

namespace Internal {

template <typename MatrixT, typename Enable = void>
struct IsSquareImpl
{
  // Default implementation
  static bool isSquare(MatrixT const & m)
  {
    return m.rows() == m.cols();
  }
};

} // namespace Internal

/** Check if a matrix is square, that is, has the same number of rows and columns. */
template <typename MatrixT>
bool
isSquare(MatrixT const & m)
{
  return Internal::IsSquareImpl<MatrixT>::isSquare(m);
}

/** Get the elements of a matrix in row-major order. */
template <typename Derived, typename OutT>
void
getElementsRowMajor(Eigen::MatrixBase<Derived> const & m, OutT * buf)
{
  alwaysAssertM(buf, "Math::getElementsRowMajor: Output buffer should be non-null");

  Eigen::Map< MatrixX<OutT, MatrixLayout::ROW_MAJOR> > bm(buf, m.rows(), m.cols());
  bm = m.template cast<OutT>();
}

/** Get the elements of a matrix in column-major order. */
template <typename Derived, typename OutT>
void
getElementsColumnMajor(Eigen::MatrixBase<Derived> const & m, OutT * buf)
{
  alwaysAssertM(buf, "Math::getElementsColumnMajor: Output buffer should be non-null");

  Eigen::Map< MatrixX<OutT, MatrixLayout::COLUMN_MAJOR> > bm(buf, m.rows(), m.cols());
  bm = m.template cast<OutT>();
}

/** Get the coordinate of a vector with the least value. The behavior is undefined if \a v is not a (row or column) vector. */
template <typename Derived>
intx
minAxis(Eigen::MatrixBase<Derived> const & v)
{
  intx min_axis;
  v.minCoeff(&min_axis);
  return min_axis;
}

/** Get the coordinate of a vector with the largest value. The behavior is undefined if \a v is not a (row or column) vector. */
template <typename Derived>
intx
maxAxis(Eigen::MatrixBase<Derived> const & v)
{
  intx max_axis;
  v.maxCoeff(&max_axis);
  return max_axis;
}

/**
 * Get the coordinate of a vector with the least absolute value. The behavior is undefined if \a v is not a (row or column)
 * vector.
 */
template <typename Derived>
intx
minAbsAxis(Eigen::MatrixBase<Derived> const & v)
{
  intx min_axis;
  v.cwiseAbs().minCoeff(&min_axis);
  return min_axis;
}

/**
 * Get the coordinate of a vector with the largest absolute value. The behavior is undefined if \a v is not a (row or column)
 * vector.
 */
template <typename Derived>
intx
maxAbsAxis(Eigen::MatrixBase<Derived> const & v)
{
  intx max_axis;
  v.cwiseAbs().maxCoeff(&max_axis);
  return max_axis;
}

/**
 * Traits class that checks whether a matrix has compile-time dimensions corresponding to a transform in homogeneous
 * coordinates.
 */
template <typename MatrixT, int N, typename Enable = void>
struct IsHomMatrix
{
  static bool const value = false;
};

template <typename MatrixT, int N>
struct IsHomMatrix< MatrixT, N, typename std::enable_if< MatrixT::RowsAtCompileTime == N + 1
                                                      && MatrixT::ColsAtCompileTime == N + 1 >::type >
{
  static bool const value = true;
};

/**
 * Convenience function to multiply an (N + 1) x (N + 1) matrix by an N-vector, by converting to homogeneous coordinates,
 * multiplying, and converting back.
 */
template < typename HomMatrixT, typename VectorT,
           typename std::enable_if< VectorT::ColsAtCompileTime == 1
                                 && HomMatrixT::RowsAtCompileTime == HomMatrixT::ColsAtCompileTime
                                 && HomMatrixT::ColsAtCompileTime == VectorT::RowsAtCompileTime + 1 >::type * = nullptr >
Matrix< VectorT::RowsAtCompileTime, 1, typename VectorT::value_type,
        HomMatrixT::IsRowMajor ? MatrixLayout::ROW_MAJOR : MatrixLayout::COLUMN_MAJOR,
        VectorT::MaxRowsAtCompileTime, VectorT::MaxColsAtCompileTime >
hmul(HomMatrixT const & m, VectorT const & v)
{
  return (m * v.homogeneous()).hnormalized();
}

/**
 * Get the one-hot vector with all entries 0 except a single entry which is 1. Values of the CoordinateAxis enum may be used as
 * arguments. This version returns a fixed-length vector, and is equivalent to <tt>Vector<N, T>::Unit(coord)</tt>.
 */
template <int N, typename T = Real> Vector<N, T> oneHot(int coord) { return Vector<N, T>::Unit(coord); }

/**
 * Get the one-hot vector with all entries 0 except a single entry which is 1. Values of the CoordinateAxis enum may be used as
 * arguments. This version initializes a provided vector to be one-hot.
 */
template <typename Derived>
void
oneHot(int coord, Eigen::MatrixBase<Derived> & v)
{
  static_assert(Derived::IsVectorAtCompileTime != 0, "Math::oneHot: Output must be compile-time vector");
  theaAssertM(coord >= 0 && coord < v.size(), "Math::oneHot: Coordinate index out of bounds");

  v.fill(0); v[coord] = 1;
}

/**
 * Get the one-cold vector with all entries 1 except a single entry which is 0. Values of the CoordinateAxis enum may be used as
 * arguments. This version creates and returns a fixed-length vector.
 */
template <int N, typename T = Real>
Vector<N, T>
oneCold(int coord)
{
  theaAssertM(coord >= 0 && coord < N, "Math::oneCold: Coordinate index out of bounds");

  Vector<N, T> v;
  v.fill(1); v[coord] = 0;
  return v;
}

/**
 * Get the one-cold vector with all entries 1 except a single entry which is 0. Values of the CoordinateAxis enum may be used as
 * arguments. This version initializes a provided vector to be one-cold.
 */
template <typename Derived>
void
oneCold(int coord, Eigen::MatrixBase<Derived> & v)
{
  static_assert(Derived::IsVectorAtCompileTime != 0, "Math::oneCold: Output must be compile-time vector");
  theaAssertM(coord >= 0 && coord < v.size(), "Math::oneCold: Coordinate index out of bounds");

  v.fill(1); v[coord] = 0;
}

/**
 * Given a 2D vector \a v, get the vector <tt>u</tt> perpendicular to it and of the same length, forming a right-handed basis
 * <tt>(u, v)</tt>.
 */
template <typename T>
Vector<2, T>
orthogonalVector(Vector<2, T> const & v)
{
  return Vector<2, T>(v[1], -v[0]);
}

/**
 * Given a 2D vector \a v, get the unit vector <tt>u</tt> perpendicular to it, forming an orthonormal right-handed basis
 * <tt>(u, v.normalized())</tt>. In other words, if \a v is the Y axis of the local frame, then the function returns the unit X
 * axis.
 */
template <typename T>
Vector<2, T>
orthogonalDirection(Vector<2, T> const & v)
{
  return orthogonalVector(v).normalized();
}

/** Given a 3D vector, get an arbitrary unit vector perpendicular to it. */
template <typename T>
Vector<3, T>
orthogonalDirection(Vector<3, T> const & v)
{
  if (maxAbsAxis(v) == 0)
    return Vector<3, T>(v[1], -v[0], 0).normalized();
  else
    return Vector<3, T>(0, v[2], -v[1]).normalized();
}

/**
 * Given a 3D vector \a w (need not be unit length), construct a 3x3 rotation matrix whose last column is the normalized
 * direction of \a w. The first two columns are arbitrarily chosen to be unit vectors perpendicular to each other and to \a w.
 * The columns define a local right-handed coordinate frame with Z axis along \a w, and X and Y axes along the first and second
 * columns respectively (X x Y == Z).
 */
template <typename T>
Matrix<3, 3, T>
orthonormalBasis(Vector<3, T> const & w)
{
  Vector<3, T>  wnrm = w.normalized();
  Vector<3, T>  u = orthogonalDirection(w);
  Vector<3, T>  v = wnrm.cross(u);

  Matrix<3, 3, T> m; m << u, v, wnrm;
  return m;
}

/**
 * Given three 3D vectors \a u, \a v and \a w (need not be unit length), construct a 3x3 rotation matrix whose last column is
 * the normalized direction of \a w, and whose first and second columns are "close to" the normalized directions of \a u and
 * \a v respectively. The columns of the returned matrix will form a right-handed orthonormal basis. If the input approximates a
 * left-handed basis, one of the first two columns of the returned matrix will oppose the corresponding input vector.
 */
template <typename T>
Matrix<3, 3, T>
orthonormalBasis(Vector<3, T> const & u, Vector<3, T> const & v, Vector<3, T> const & w)
{
  Vector<3, T> w2 = w.normalized();
  Vector<3, T> u2 = v.cross(w2).normalized();
  Vector<3, T> v2 = w2.cross(u2);

  Matrix<3, 3, T> m; m << u2, v2, w2;
  return m;
}

/**
 * Given an arbitrary input 3x3 matrix \a m, construct another 3x3 rotation matrix whose last column is
 * the normalized direction of the last column of \a m, and whose first and second columns are "close to" the normalized
 * directions of the first and second columns respectively of \a m. The columns of the returned matrix will form a right-handed
 * orthonormal basis. If the input approximates a left-handed basis, one of the first two columns of the returned matrix will
 * oppose the corresponding input column.
 */
template <typename T>
Matrix<3, 3, T>
orthonormalBasis(Matrix<3, 3, T> const & m)
{
  return orthonormalBasis(Vector<3, T>(m.col(0)), Vector<3, T>(m.col(1)), Vector<3, T>(m.col(2)));
}

/**
 * Matrix to scale a point by the corresponding scaling parameter along each dimension.
 *
 * @param s The vector of scaling factors per dimension.
 */
template <int N, typename T>
Matrix<N, N, T>
scaling(Vector<N, T> const & s)
{
  return Matrix<N, N, T>(s.asDiagonal());
}

/** Matrix to uniformaly scale a point by a scaling factor \a s. */
template < int N, typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<N, N, T>
scaling(T const & s)
{
  Vector<N, T> v; v.fill(s);
  return scaling(v);
}

/** Matrix to rotate a 2D vector about the origin by an angle (in radians). */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<2, 2, T>
rotation(T const & radians)
{
  T s = std::sin(radians);
  T c = std::cos(radians);

  Matrix<2, 2, T> m; m << c, -s,
                          s,  c;
  return m;
}

/** Rotate around the given 3D axis (need not be a unit vector) by a given angle. */
template <typename T>
Matrix<3, 3, T>
rotationAxisAngle(Vector<3, T> const & axis, T const & radians)
{
  Vector<3, T> uaxis = axis.normalized();

  T cos_val = std::cos(radians);
  T sin_val = std::sin(radians);
  T one_minus_cos = 1 - cos_val;
  T x2   = uaxis[0] * uaxis[0];
  T y2   = uaxis[1] * uaxis[1];
  T z2   = uaxis[2] * uaxis[2];
  T xym  = uaxis[0] * uaxis[1] * one_minus_cos;
  T xzm  = uaxis[0] * uaxis[2] * one_minus_cos;
  T yzm  = uaxis[1] * uaxis[2] * one_minus_cos;
  T xsin = uaxis[0] * sin_val;
  T ysin = uaxis[1] * sin_val;
  T zsin = uaxis[2] * sin_val;

  Matrix<3, 3, T> m;
  m(0, 0) = x2 * one_minus_cos + cos_val;
  m(0, 1) = xym - zsin;
  m(0, 2) = xzm + ysin;

  m(1, 0) = xym + zsin;
  m(1, 1) = y2 * one_minus_cos + cos_val;
  m(1, 2) = yzm - xsin;

  m(2, 0) = xzm - ysin;
  m(2, 1) = yzm + xsin;
  m(2, 2) = z2 * one_minus_cos + cos_val;

  return m;
}

/**
 * Rotate about the Z axis by the roll angle, then the Y axis by the pitch angle, and finally the X axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesXYZ(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0.0, sin_val, cos_val;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  return xmat * (ymat * zmat);
}

/**
 * Rotate about the Y axis by the roll angle, then the Z axis by the pitch angle, and finally the X axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesXZY(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  return xmat * (zmat * ymat);
}

/**
 * Rotate about the Z axis by the roll angle, then the X axis by the pitch angle, and finally the Y axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesYXZ(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  return ymat * (xmat * zmat);
}

/**
 * Rotate about the X axis by the roll angle, then the Z axis by the pitch angle, and finally the Y axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesYZX(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  return ymat * (zmat * xmat);
}

/**
 * Rotate about the Y axis by the roll angle, then the X axis by the pitch angle, and finally the Z axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesZXY(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  return zmat * (xmat * ymat);
}

/**
 * Rotate about the X axis by the roll angle, then the Y axis by the pitch angle, and finally the Z axis by the yaw angle.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<3, 3, T>
rotationEulerAnglesZYX(T const & yaw_radians, T const & pitch_radians, T const & roll_radians)
{
  T cos_val, sin_val;

  cos_val = std::cos(yaw_radians);
  sin_val = std::sin(yaw_radians);
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = std::cos(pitch_radians);
  sin_val = std::sin(pitch_radians);
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = std::cos(roll_radians);
  sin_val = std::sin(roll_radians);
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  return zmat * (ymat * xmat);
}

/**
 * Return the quaternion corresponding to the rotation from one direction vector to another.
 *
 * Use this instead of Eigen's Quaternion::FromTwoVectors() when the input vectors are pre-normalized, since the latter function
 * will always normalize the input vectors, wasting two sqrt's. However, for most other practical purposes the two functions
 * should have the same behavior.
 *
 * @param start_dir The vector to rotate from.
 * @param end_dir The vector to rotate to.
 * @param normalize_dirs If false, the directions will be assumed to have been pre-normalized to unit length before being
 *   passed to this function.
 */
template <typename T>
Quaternion<T>
rotationArcQuat(Vector<3, T> const & start_dir, Vector<3, T> const & end_dir, bool normalize_dirs = true)
{
  // From John Ratcliff's Code Suppository.
  //
  // Reference, from Stan Melax in Game Gems I
  //  Quaternion q;
  //  vector3 c = CrossProduct(v0,v1);
  //  REAL d = DotProduct(v0,v1);
  //  REAL s = (REAL)sqrt((1+d)*2);
  //  q.x = c.x / s;
  //  q.y = c.y / s;
  //  q.z = c.z / s;
  //  q.w = s / 2.0f;
  //  return q;

  Vector<3, T> u = start_dir, v = end_dir;
  if (normalize_dirs)
  {
    u.normalize();
    v.normalize();
  }

  T d = u.dot(v);
  if (d < -0.9999)
  {
    Vector<3, T> perp;
    switch (maxAbsAxis(u))
    {
      case 0:           perp = Vector<3, T>(u.y(),  -u.x(),       0).normalized(); break;
      case 1:           perp = Vector<3, T>(u.y(),  -u.x(),       0).normalized(); break;
      default /* 2 */:  perp = Vector<3, T>(u.z(),       0,  -u.x()).normalized();
    }

    return Quaternion<T>(Eigen::AngleAxis<T>(static_cast<T>(Math::pi()), perp));
  }

  T s = std::sqrt((1 + d) * 2);
  T recip = 1 / s;
  Vector<3, T> rcc = recip * u.cross(v);

  return Quaternion<T>(static_cast<T>(0.5 * s), rcc.x(), rcc.y(), rcc.z());
}

/**
 * Return the matrix corresponding to the rotation from one direction vector to another.
 *
 * @param start_dir The vector to rotate from.
 * @param end_dir The vector to rotate to.
 * @param normalize_dirs If false, the directions will be assumed to have been pre-normalized to unit length before being
 *   passed to this function.
 */
template <typename T>
Matrix<3, 3, T>
rotationArc(Vector<3, T> const & start_dir, Vector<3, T> const & end_dir, bool normalize_dirs = true)
{
  return rotationArcQuat(start_dir, end_dir, normalize_dirs).toRotationMatrix();
}

/**
 * Constructs an orthogonal 3D projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
 * the  near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases downwards
 * instead, e.g. for screen pixel space.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<4, 4, T>
orthogonalProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval, T const & farval,
                     bool y_increases_upwards = true)
{
  // Adapted from Mesa.
  // Note that Microsoft (http://msdn.microsoft.com/library/default.asp?url=/library/en-us/opengl/glfunc03_8qnj.asp) and
  // Linux (http://www.xfree86.org/current/glOrtho.3.html) have different matrices shown in their documentation.

  alwaysAssertM(left < right,     format("Math::orthogonalProjection: Left (%lf) must be less than right (%lf)",
                                         (double)left, (double)right));
  alwaysAssertM(bottom < top,     format("Math::orthogonalProjection: Bottom (%lf) must be less than top (%lf)",
                                         (double)bottom, (double)top));
  alwaysAssertM(nearval < farval, format("Math::orthogonalProjection: Near (%lf) must be less than far (%lf)",
                                         (double)nearval, (double)farval));

  T x =  2 / (right - left);
  T y =  2 / (top - bottom);
  T z = -2 / (farval - nearval);
  T tx = -(right + left) / (right - left);
  T ty = -(top + bottom) / (top - bottom);
  T tz = -(farval + nearval) / (farval - nearval);

  if (!y_increases_upwards)
  {
    y  = -y;
    ty = -ty;
  }

  return (Matrix<4, 4, T>() << x, 0, 0, tx,
                               0, y, 0, ty,
                               0, 0, z, tz,
                               0, 0, 0, 1).finished();
}

/**
 * Infer the parameters of a 3D orthographic projection, given the final projection matrix which must be exactly consistent with
 * the conventions of orthographicProjection().
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
bool
inferOrthogonalProjectionParams(Matrix<4, 4, T> const & m,
                                T & left, T & right, T & bottom, T & top, T & nearval, T & farval, bool & y_increases_upwards)
{
  // Given the matrix:
  //
  // x 0 0 tx
  // 0 y 0 ty
  // 0 0 z tz
  // 0 0 0 1
  //
  // ... solve the system of linear equations for (L, R, B, T, N, F):
  //
  //   -Lx + Rx = 2
  //   -By + Ty = 2
  //   -Nz + Fz = -2
  //   L(1 - tx) + R(1 + tx) = 0
  //   B(1 - ty) + T(1 + ty) = 0
  //   N(1 - tz) + F(1 + tz) = 0

  if (Math::fuzzyNe(m(0, 1), static_cast<T>(0)) || Math::fuzzyNe(m(0, 2), static_cast<T>(0))
   || Math::fuzzyNe(m(1, 0), static_cast<T>(0)) || Math::fuzzyNe(m(1, 2), static_cast<T>(0))
   || Math::fuzzyNe(m(2, 0), static_cast<T>(0)) || Math::fuzzyNe(m(2, 1), static_cast<T>(0))
   || Math::fuzzyNe(m(3, 0), static_cast<T>(0)) || Math::fuzzyNe(m(3, 1), static_cast<T>(0))
   || Math::fuzzyNe(m(3, 2), static_cast<T>(0)) || Math::fuzzyNe(m(3, 3), static_cast<T>(1)))
  {
    THEA_ERROR << "Math::inferOrthogonalProjectionParams: Matrix is not an orthogonal projection in the expected form";
    return false;
  }

  Matrix<6, 6, T> coeffs; coeffs.setZero();
  Vector<6, T> constants;

  coeffs(0, 0) = -m(0, 0);    coeffs(0, 1) = m(0, 0);     constants[0] = 2;
  coeffs(1, 2) = -m(1, 1);    coeffs(1, 3) = m(1, 1);     constants[1] = 2;
  coeffs(2, 4) = -m(2, 2);    coeffs(2, 5) = m(2, 2);     constants[2] = -2;
  coeffs(3, 0) = 1 - m(0, 3); coeffs(3, 1) = 1 + m(0, 3); constants[3] = 0;
  coeffs(4, 2) = 1 - m(1, 3); coeffs(4, 3) = 1 + m(1, 3); constants[4] = 0;
  coeffs(5, 4) = 1 - m(2, 3); coeffs(5, 5) = 1 + m(2, 3); constants[5] = 0;

  Vector<6, T> sol = coeffs.colPivHouseholderQr().solve(constants);
  if (!(coeffs * sol).isApprox(constants))
  {
    THEA_ERROR << "Math::inferOrthogonalProjectionParams: Could not solve linear system for projection parameters";
    return false;
  }

  left = sol[0]; right = sol[1]; bottom = sol[2]; top = sol[3]; nearval = sol[4]; farval = sol[5];

  y_increases_upwards = true;
  if (bottom > top) { std::swap(bottom, top); y_increases_upwards = false; }

  if (left >= right)
  { THEA_ERROR << "Math::inferOrthogonalProjectionParams: Left not less than right -- check the matrix"; return false; }

  if (bottom >= top)
  { THEA_ERROR << "Math::inferOrthogonalProjectionParams: Bottom not less than top -- check the matrix"; return false; }

  if (nearval >= farval)
  { THEA_ERROR << "Math::inferOrthogonalProjectionParams: Near not less than far -- check the matrix"; return false; }

  return true;
}

/**
 * Constructs a 3D perspective projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
 * the near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases downwards
 * instead, e.g. for screen pixel space.
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
Matrix<4, 4, T>
perspectiveProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval, T const & farval,
                      bool y_increases_upwards = true)
{
  alwaysAssertM(left < right,     format("Math::perspectiveProjection: Left (%lf) must be less than right (%lf)",
                                         (double)left, (double)right));
  alwaysAssertM(bottom < top,     format("Math::perspectiveProjection: Bottom (%lf) must be less than top (%lf)",
                                         (double)bottom, (double)top));
  alwaysAssertM(nearval < farval, format("Math::perspectiveProjection: Near (%lf) must be less than far (%lf)",
                                         (double)nearval, (double)farval));
  alwaysAssertM(nearval > 0,      format("Math::perspectiveProjection: Near (%lf) must be positive", (double)nearval));

  T x = (2 * nearval) / (right - left);
  T y = (2 * nearval) / (top - bottom);
  T a = (right + left) / (right - left);
  T b = (top + bottom) / (top - bottom);

  T c, d;
  if (Math::isInfinite(farval))
  {
    // Infinite view frustum
    c = -1;
    d = -2 * nearval;
  }
  else
  {
    c = -(farval + nearval) / (farval - nearval);
    d = (-2 * farval * nearval) / (farval - nearval);
  }

  if (!y_increases_upwards)
  {
    y = -y;
    b = -b;
  }

  return (Matrix<4, 4, T>() << x,  0,  a,  0,
                               0,  y,  b,  0,
                               0,  0,  c,  d,
                               0,  0, -1,  0).finished();
}

/**
 * Infer the parameters of a 3D perspective projection, given the final projection matrix which must be exactly consistent with
 * the conventions of perspectiveProjection().
 */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
bool
inferPerspectiveProjectionParams(Matrix<4, 4, T> const & m,
                                 T & left, T & right, T & bottom, T & top, T & nearval, T & farval, bool & y_increases_upwards)
{
  // Given the matrix:
  //
  // x  0  a  0
  // 0  y  b  0
  // 0  0  c  d
  // 0  0 -1  0
  //
  // ... solve the system of equations for (L, R, B, T, N, F):
  //
  //   -Lx + Rx - 2N = 0
  //   -By + Ty - 2N = 0
  //   -L(a + 1) + R(a - 1) = 0
  //   -B(b + 1) + T(b - 1) = 0
  //   N(1 - c) + F(1 + c) = 0
  //   -Nd + Fd + 2FN = 0   <-- quadratic!
  //
  // ... by first solving this quadratic equation for its non-zero solution:
  //
  // Nd (1 + (1 - c) / (1 + c)) + N^2 2(1 - c) / (1 + c) = 0

  if (Math::fuzzyNe(m(0, 1), static_cast<T>(0))  || Math::fuzzyNe(m(0, 3), static_cast<T>(0))
   || Math::fuzzyNe(m(1, 0), static_cast<T>(0))  || Math::fuzzyNe(m(1, 3), static_cast<T>(0))
   || Math::fuzzyNe(m(2, 0), static_cast<T>(0))  || Math::fuzzyNe(m(2, 1), static_cast<T>(0))
   || Math::fuzzyNe(m(3, 0), static_cast<T>(0))  || Math::fuzzyNe(m(3, 1), static_cast<T>(0))
   || Math::fuzzyNe(m(3, 2), static_cast<T>(-1)) || Math::fuzzyNe(m(3, 3), static_cast<T>(0)))
  {
    THEA_ERROR << "Math::inferPerspectiveProjectionParams: Matrix is not a perspective projection in the expected form";
    return false;
  }

  if (Math::fuzzyEq(m(2, 2), static_cast<T>(-1)))  // infinite view frustum
  {
    nearval = m(2, 3) / (-2);
    farval = Math::inf<T>();
  }
  else
  {
    auto x = (1 - m(2, 2)) / (1 + m(2, 2));
    nearval = -m(2, 3) * (1 + x) / (2 * x);
    farval = -nearval * x;
  }

  if (nearval <= 0)
  { THEA_ERROR << "Math::inferPerspectiveProjectionParams: Near not positive -- check the matrix"; return false; }

  if (nearval >= farval)
  { THEA_ERROR << "Math::inferPerspectiveProjectionParams: Near not less than far -- check the matrix"; return false; }

  Matrix<4, 4, T> coeffs;
  coeffs << -m(0, 0),        m(0, 0),        0,              0,
             0,              0,             -m(1, 1),        m(1, 1),
            -(m(0, 2) + 1),  (m(0, 2) - 1),  0,              0,
             0,              0,             -(m(1, 2) + 1),  (m(1, 2) - 1);
  Vector<4, T> constants;
  constants << 2 * nearval, 2 * nearval, 0, 0;

  Vector<4, T> sol = coeffs.colPivHouseholderQr().solve(constants);
  if (!(coeffs * sol).isApprox(constants))
  {
    THEA_ERROR << "Math::inferPerspectiveProjectionParams: Could not solve linear system for projection parameters";
    return false;
  }

  left = sol[0]; right = sol[1]; bottom = sol[2]; top = sol[3];;

  y_increases_upwards = true;
  if (bottom > top) { std::swap(bottom, top); y_increases_upwards = false; }

  if (left >= right)
  { THEA_ERROR << "Math::inferPerspectiveProjectionParams: Left must be less than right -- check the matrix"; return false; }

  if (bottom >= top)
  { THEA_ERROR << "Math::inferPerspectiveProjectionParams: Bottom must be less than top -- check the matrix"; return false; }

  return true;
}

/**
 * Solve for the real eigenvalues and eigenvectors of a 2x2 matrix.
 *
 * @param m The matrix whose eigenvalues/vectors will be found.
 * @param eigenvalues Used to return the eigenvalues of the matrix. Must be preallocated to at least 2 elements.
 * @param eigenvectors Used to return the eigenvectors of the matrix. Must be preallocated to at least 2 elements.
 * @param tol Numerical tolerance, negative for default.
 *
 * @return The number of real eigenvalues found.
 */
template <typename T>
int
eigenSolve(Matrix<2, 2, T> const & m, T * eigenvalues, Vector<2, T> * eigenvectors, T const & tol = -1)
{
  T a = m(0, 0), b = m(0, 1);
  T c = m(1, 0), d = m(1, 1);

  T trace  =  a + d;
  T det    =  a * d - b * c;

  T disc = trace * trace / 4 - det;
  if (disc < 0)
    return 0;

  T s = std::sqrt(disc);
  eigenvalues[0] = trace / 2 - s;
  eigenvalues[1] = trace / 2 + s;

  if (!Math::fuzzyEq(c, static_cast<T>(0), (tol >= 0 ? tol : Math::eps(c, static_cast<T>(0)))))
  {
    eigenvectors[0][0] = eigenvalues[0] - d;
    eigenvectors[0][1] = c;

    eigenvectors[1][0] = eigenvalues[1] - d;
    eigenvectors[1][1] = c;
  }
  else if (!Math::fuzzyEq(b, static_cast<T>(0), (tol >= 0 ? tol : Math::eps(b, static_cast<T>(0)))))
  {
    eigenvectors[0][0] = b;
    eigenvectors[0][1] = eigenvalues[0] - a;

    eigenvectors[1][0] = b;
    eigenvectors[1][1] = eigenvalues[1] - a;
  }
  else
  {
    eigenvectors[0][0] = 1;
    eigenvectors[0][1] = 0;

    eigenvectors[1][0] = 0;
    eigenvectors[1][1] = 1;
  }

  return 2;
}

} // namespace Math

/**
 * Get a single-line string representation of a matrix or vector. If the number of rows (or columns) of the matrix is larger
 * than \a max_rows (resp. \a max_cols), the middle elements will be elided via ellipsis.
 *
 * \code
 * Matrix3 m; m << 1, 2, 3,
 *                 4, 5, 6,
 *                 7, 8, 9;
 *
 * cout << toString(m, 3, 3) << endl;
 *    // Output: [1 2 3; 4 5 6; 7 8 9]
 * cout << toString(m, 2, 2) << endl;
 *    // Output: [1 ... 3; ... ; 7 ... 9]
 * cout << toString(m, 2, 1) << endl;
 *    // Output: [1 ... ; ... ; 7 ... ]
 * cout << toString(m, 1, 2) << endl;
 *    // Output: [1 ... 3; ... ]
 * cout << toString(m, 1, 1) << endl;
 *    // Output: [1 ... ; ... ]
 * \endcode
 */
template <typename Derived>
std::string
toString(Eigen::DenseBase<Derived> const & m, intx max_rows = 4, intx max_cols = 4)
{
  intx first_rows = m.rows(), last_rows = 0;
  intx first_cols = m.cols(), last_cols = 0;
  if (m.rows() > max_rows)
  {
    first_rows = max_rows / 2 + max_rows % 2;
    last_rows  = max_rows / 2;
  }
  if (m.cols() > max_cols)
  {
    first_cols = max_cols / 2 + max_cols % 2;
    last_cols  = max_cols / 2;
  }

  std::ostringstream oss;
  oss << '[';

  for (int i = 0; i < 2; ++i)
  {
    intx row_begin = (i == 0 ? 0 : m.rows() - last_rows);
    intx row_end   = (i == 0 ? first_rows : m.rows());

    for (intx r = row_begin; r < row_end; ++r)
    {
      for (int j = 0; j < 2; ++j)
      {
        intx col_begin = (j == 0 ? 0 : m.cols() - last_cols);
        intx col_end   = (j == 0 ? first_cols : m.cols());

        for (intx c = col_begin; c < col_end; ++c)
        {
          if (c > col_begin) oss << ' ';
          oss << m(r, c);
        }

        if (j == 0 && (first_cols + last_cols) < m.cols())
          oss << " ... ";
      }

      if (r != m.rows() - 1) oss << "; ";
    }

    if (i == 0)
    {
      if ((first_rows + last_rows) < m.rows()) oss << " ... ";
      if (last_rows > 0) oss << "; ";
    }
  }

  oss << ']';
  return oss.str();
}

} // namespace Thea

#endif
