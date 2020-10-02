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
Eigen::Index
minAxis(Eigen::MatrixBase<Derived> const & v)
{
  Eigen::Index min_axis;
  v.minCoeff(&min_axis);
  return min_axis;
}

/** Get the coordinate of a vector with the largest value. The behavior is undefined if \a v is not a (row or column) vector. */
template <typename Derived>
Eigen::Index
maxAxis(Eigen::MatrixBase<Derived> const & v)
{
  Eigen::Index max_axis;
  v.maxCoeff(&max_axis);
  return max_axis;
}

/**
 * Get the coordinate of a vector with the least absolute value. The behavior is undefined if \a v is not a (row or column)
 * vector.
 */
template <typename Derived>
Eigen::Index
minAbsAxis(Eigen::MatrixBase<Derived> const & v)
{
  Eigen::Index min_axis;
  v.cwiseAbs().minCoeff(&min_axis);
  return min_axis;
}

/**
 * Get the coordinate of a vector with the largest absolute value. The behavior is undefined if \a v is not a (row or column)
 * vector.
 */
template <typename Derived>
Eigen::Index
maxAbsAxis(Eigen::MatrixBase<Derived> const & v)
{
  Eigen::Index max_axis;
  v.cwiseAbs().maxCoeff(&max_axis);
  return max_axis;
}

/**
 * Convenience function to multiply a 4x4 matrix by a 3-vector, by converting to homogeneous coordinates, multiplying, and
 * converting back.
 */
template <typename T, int N, int O1, int R1, int C1, int O2, int R2, int C2>
Eigen::Matrix<T, N - 1, 1, O2, R2, C2>
hmul(Eigen::Matrix<T, N,     N, O1, R1, C1> const & m,
     Eigen::Matrix<T, N - 1, 1, O2, R2, C2> const & v)
{
  return (m * v.homogeneous()).hnormalized();
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
 * Given a 3D vector \a w, construct a rotation matrix whose last column is the normalized direction of \a w. The first two
 * columns are arbitrarily chosen to be unit vectors perpendicular to each other and to \a w. Thus, the columns define a local
 * coordinate frame with Z axis along \a w.
 */
template <typename T>
Matrix<3, 3, T>
orthonormalBasis(Vector<3, T> const & w)
{
  Vector<3, T>  wnrm = w.normalized();
  Vector<3, T>  u = orthogonalDirection(w);
  Vector<3, T>  v = wnrm.cross(u);

  Matrix<3, 3, T> m;
  m << u, v, wnrm;
  return m;
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
template <int N, typename T>
Matrix<N, N, T>
scaling(T const & s)
{
  Vector<N, T> v; v.fill(s);
  return scaling(v);
}

/** Matrix to rotate a 2D vector about the origin by an angle (in radians). */
template <typename T>
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
rotationAxisAngle(Vector<3, T> const & axis, Real radians)
{
  Vector<3, T> uaxis = axis.normalized();

  T cos_val = static_cast<T>(std::cos(radians));
  T sin_val = static_cast<T>(std::sin(radians));
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
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesXYZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0.0, sin_val, cos_val;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  return xmat * (ymat * zmat);
}

/**
 * Rotate about the Y axis by the roll angle, then the Z axis by the pitch angle, and finally the X axis by the yaw angle.
 */
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesXZY(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  return xmat * (zmat * ymat);
}

/**
 * Rotate about the Z axis by the roll angle, then the X axis by the pitch angle, and finally the Y axis by the yaw angle.
 */
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesYXZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  return ymat * (xmat * zmat);
}

/**
 * Rotate about the X axis by the roll angle, then the Z axis by the pitch angle, and finally the Y axis by the yaw angle.
 */
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesYZX(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  return ymat * (zmat * xmat);
}

/**
 * Rotate about the Y axis by the roll angle, then the X axis by the pitch angle, and finally the Z axis by the yaw angle.
 */
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesZXY(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  return zmat * (xmat * ymat);
}

/**
 * Rotate about the X axis by the roll angle, then the Y axis by the pitch angle, and finally the Z axis by the yaw angle.
 */
template <typename T>
Matrix<3, 3, T>
rotationEulerAnglesZYX(Real yaw_radians, Real pitch_radians, Real roll_radians)
{
  T cos_val, sin_val;

  cos_val = static_cast<T>(std::cos(yaw_radians));
  sin_val = static_cast<T>(std::sin(yaw_radians));
  Matrix<3, 3, T> zmat; zmat << cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1;

  cos_val = static_cast<T>(std::cos(pitch_radians));
  sin_val = static_cast<T>(std::sin(pitch_radians));
  Matrix<3, 3, T> ymat; ymat << cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val;

  cos_val = static_cast<T>(std::cos(roll_radians));
  sin_val = static_cast<T>(std::sin(roll_radians));
  Matrix<3, 3, T> xmat; xmat << 1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val;

  return zmat * (ymat * xmat);
}

/** Construct a rotation matrix from a quaternion. */
template <typename T>
Matrix<3, 3, T>
rotationQuat(T const & x, T const & y, T const & z, T const & w)
{
  // From G3D::Matrix3:
  // Implementation from Watt and Watt, pg 362
  // See also http://www.flipcode.com/documents/matrfaq.html#Q54

  Vector<4, T> q(x, y, z, w);
  q.normalize();

  T xx = 2 * q.x() * q.x();
  T xy = 2 * q.x() * q.y();
  T xz = 2 * q.x() * q.z();
  T xw = 2 * q.x() * q.w();

  T yy = 2 * q.y() * q.y();
  T yz = 2 * q.y() * q.z();
  T yw = 2 * q.y() * q.w();

  T zz = 2 * q.z() * q.z();
  T zw = 2 * q.z() * q.w();

  Matrix<3, 3, T> ret;
  ret << 1 - yy - zz,      xy - zw,      xz + yw,
             xy + zw,  1 - xx - zz,      yz - xw,
             xz - yw,      yz + xw,  1 - xx - yy;

  return ret;
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
      case 0:           perp = Vector<3, T>(u.y(),  -u.x(),       0); break;
      case 1:           perp = Vector<3, T>(u.y(),  -u.x(),       0); break;
      default /* 2 */:  perp = Vector<3, T>(u.z(),       0,  -u.x());
    }

    return rotationAxisAngle(perp, (Real)Math::pi());
  }

  T s = std::sqrt((1 + d) * 2);
  T recip = 1 / s;
  Vector<3, T> rcc = recip * u.cross(v);

  return rotationQuat(rcc.x(), rcc.y(), rcc.z(), static_cast<T>(0.5 * s));
}

/**
 * Constructs an orthogonal 3D projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
 * the  near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases downwards
 * instead, e.g. for screen pixel space.
 */
template <typename T>
Matrix<4, 4, T>
orthogonalProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval, T const & farval,
                     bool y_increases_upwards = true)
{
  // Adapted from Mesa.
  // Note that Microsoft (http://msdn.microsoft.com/library/default.asp?url=/library/en-us/opengl/glfunc03_8qnj.asp) and
  // Linux (http://www.xfree86.org/current/glOrtho.3.html) have different matrices shown in their documentation.

  T x =  2 / (right - left);
  T y =  2 / (top - bottom);
  T z = -2 / (farval - nearval);
  T tx = -(right + left) / (right - left);
  T ty = -(top + bottom) / (top - bottom);
  T tz = -(farval + nearval) / (farval - nearval);

  if (!y_increases_upwards)
  {
    y  *= -1;
    ty *= -1;
  }

  return (Matrix<4, 4, T>() << x, 0, 0, tx,
                               0, y, 0, ty,
                               0, 0, z, tz,
                               0, 0, 0, 1).finished();
}

/**
 * Constructs a 3D perspective projection matrix from the given parameters. \a nearval and \a farval are the <i>negative</i> of
 * the near and far plane Z values (to follow OpenGL conventions). Set \a y_increases_upwards to false if Y increases downwards
 * instead, e.g. for screen pixel space.
 */
template <typename T>
Matrix<4, 4, T>
perspectiveProjection(T const & left, T const & right, T const & bottom, T const & top, T const & nearval, T const & farval,
                      bool y_increases_upwards = true)
{
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
    d = -(2 * farval * nearval) / (farval - nearval);
  }

  if (!y_increases_upwards)
  {
    y *= -1;
    b *= -1;
  }

  return (Matrix<4, 4, T>() << x,  0,  a,  0,
                               0,  y,  b,  0,
                               0,  0,  c,  d,
                               0,  0, -1,  0).finished();
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
