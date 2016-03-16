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

#ifndef __Thea_Matrix3_hpp__
#define __Thea_Matrix3_hpp__

#include "Math.hpp"
#include "Matrix2.hpp"
#include "MatrixMN.hpp"
#include "VectorN.hpp"
#include <limits>

namespace Thea {

/**
 * Square 3 x 3 matrices on a field T. The matrices are stored internally in row-major form, so row-major access is recommended.
 *
 * @note Several of the functions in this class are taken from the G3D Matrix3 class.
 */
template <typename T>
class /* THEA_API */ MatrixMN<3, 3, T> : public Internal::SquareMatrixN<3, T>
{
  private:
    typedef Internal::SquareMatrixN<3, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /** Initialize all 9 components of the matrix. */
    MatrixMN(T const & m00, T const & m01, T const & m02,
             T const & m10, T const & m11, T const & m12,
             T const & m20, T const & m21, T const & m22)
    {
      (*this)(0, 0) = m00; (*this)(0, 1) = m01; (*this)(0, 2) = m02;
      (*this)(1, 0) = m10; (*this)(1, 1) = m11; (*this)(1, 2) = m12;
      (*this)(2, 0) = m20; (*this)(2, 1) = m21; (*this)(2, 2) = m22;
    }

    /**
     * Initialize from the upper 2x2 submatrix and the first 2 entries of the last column ("translation" terms). The last row is
     * set to [0, 0, 1].
     */
    MatrixMN(MatrixMN<2, 2, T> const & u2x2, VectorN<2, T> const & trans)
    {
      (*this)(0, 0) = u2x2(0, 0); (*this)(0, 1) = u2x2(0, 1); (*this)(0, 2) = trans[0];
      (*this)(1, 0) = u2x2(1, 0); (*this)(1, 1) = u2x2(1, 1); (*this)(1, 2) = trans[1];
      (*this)(2, 0) = 0; (*this)(2, 1) = 0; (*this)(2, 2) = 1;
    }

    /** Get the matrix of cofactors. */
    MatrixMN cofactor() const
    {
      return adjoint().transpose();
    }

    /** Get the adjoint of the matrix. */
    MatrixMN adjoint() const
    {
      MatrixMN adj;
      adj(0, 0) = (*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1);
      adj(0, 1) = (*this)(0, 2) * (*this)(2, 1) - (*this)(0, 1) * (*this)(2, 2);
      adj(0, 2) = (*this)(0, 1) * (*this)(1, 2) - (*this)(0, 2) * (*this)(1, 1);
      adj(1, 0) = (*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2);
      adj(1, 1) = (*this)(0, 0) * (*this)(2, 2) - (*this)(0, 2) * (*this)(2, 0);
      adj(1, 2) = (*this)(0, 2) * (*this)(1, 0) - (*this)(0, 0) * (*this)(1, 2);
      adj(2, 0) = (*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0);
      adj(2, 1) = (*this)(0, 1) * (*this)(2, 0) - (*this)(0, 0) * (*this)(2, 1);
      adj(2, 2) = (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);

      return adj;
    }

    /** Get the determinant of the matrix. */
    T determinant() const
    {
      // Determinant is the dot product of the first row and the first row of cofactors (the first column of the adjoint matrix)

      T cofactor00 = (*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1);
      T cofactor10 = (*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2);
      T cofactor20 = (*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0);

      return (*this)(0, 0) * cofactor00 + (*this)(0, 1) * cofactor10 + (*this)(0, 2) * cofactor20;
    }

    /**
     * Invert the matrix in-place.
     *
     * @param tolerance The numerical tolerance of the computation. If the determinant has smaller absolute value than the
     *   tolerance, the computation is aborted and false is returned.
     *
     * @return True if the computation succeeded with the given tolerance, else false.
     */
    bool invert(double tolerance = 1.0e-30)
    {
      // Invert a 3x3 matrix using cofactors. This is about 8 times faster than the Numerical Recipes code which uses Gaussian
      // elimination.

      MatrixMN adj = adjoint();
      Real det = (*this)(0, 0) * adj(0, 0) + (*this)(0, 1) * adj(1, 0) + (*this)(0, 2) * adj(2, 0);
      if (std::abs(det) <= tolerance)
        return false;

      Real inv_det = 1 / det;
      for (long r = 0; r < 3; ++r)
        for (long c = 0; c < 3; ++c)
          (*this)(r, c) = adj(r, c) * inv_det;

      return true;
    }

    /**
     * Get the inverse of the matrix.
     *
     * @param tolerance The numerical tolerance of the computation. If the determinant has smaller absolute value than the
     *   tolerance, the computation is aborted and an error is thrown.
     */
    MatrixMN inverse(double tolerance = 1.0e-30) const
    {
      MatrixMN result = *this;
      if (!result.invert(tolerance))
        throw Error("MatrixMN<3, 3, T>: Could not invert matrix " + this->toString() + " with given tolerance");

      return result;
    }

    /** Get the upper 2x2 submatrix. */
    MatrixMN<2, 2, T> upper2x2() const
    {
      return MatrixMN<2, 2, T>((*this)(0, 0), (*this)(0, 1),
                               (*this)(1, 0), (*this)(1, 1));
    }

    /**
     * Solve for the real eigenvalues and eigenvectors of the matrix, assuming it is symmetric.
     *
     * @param eigenvalues Used to return the eigenvalues of the matrix. Must be preallocated to at least 3 elements.
     * @param eigenvectors Used to return the eigenvectors of the matrix. Must be preallocated to at least 3 elements.
     *
     * @return The number of real eigenvalues found.
     */
    int eigenSolveSymmetric(T * eigenvalues, VectorN<3, T> * eigenvectors) const
    {
      // From G3D::Matrix3

      MatrixMN scratch = *this;
      T subdiag[3];
      scratch.tridiagonal(eigenvalues, subdiag);
      if (!scratch.qlAlgorithm(eigenvalues, subdiag))
      {
        THEA_ERROR << "MatrixMN<3, 3, T>: Could not apply QL algorithm to matrix";
        return 0;
      }

      for (int i = 0; i < 3; i++)
        eigenvectors[i] = scratch.getColumn(i);

      // make eigenvectors form a right--handed system
      VectorN<3, T> cross = eigenvectors[1].cross(eigenvectors[2]);
      T det = eigenvectors[0].dot(cross);
      if (det < 0)
        eigenvectors[2] = -eigenvectors[2];

      return 3;
    }

    /**
     * Create a matrix from its columns.
     *
     * @param cv0 First column of the matrix.
     * @param cv1 Second column of the matrix.
     * @param cv2 Third column of the matrix.
     */
    static MatrixMN fromColumns(VectorN<3, T> const & cv0, VectorN<3, T> const & cv1, VectorN<3, T> const & cv2)
    {
      return MatrixMN(cv0[0], cv1[0], cv2[0],
                      cv0[1], cv1[1], cv2[1],
                      cv0[2], cv1[2], cv2[2]);
    }

    /**
     * Create a matrix from its rows.
     *
     * @param rv0 First row of the matrix.
     * @param rv1 Second row of the matrix.
     * @param rv2 Third row of the matrix.
     */
    static MatrixMN fromRows(VectorN<3, T> const & rv0, VectorN<3, T> const & rv1, VectorN<3, T> const & rv2)
    {
      return MatrixMN(rv0[0], rv0[1], rv0[2],
                      rv1[0], rv1[1], rv1[2],
                      rv2[0], rv2[1], rv2[2]);
    }

    /** 2D scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(T const & sx, T const & sy)
    {
      return MatrixMN(sx,  0,  0,
                       0, sy,  0,
                       0,  0,  1);
    }

    /** 2D scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(VectorN<2, T> const & v)
    {
      return homScaling(v[0], v[1]);
    }

    /** 2D uniform scaling matrix in homogeneous coordinates. */
    static MatrixMN homScaling(T const & s)
    {
      return homScaling(s, s);
    }

    /** 2D translation matrix in homogeneous coordinates. */
    static MatrixMN homTranslation(T const & tx, T const & ty)
    {
      return MatrixMN(1, 0, tx,
                      0, 1, ty,
                      0, 0, 1);
    }

    /** 2D translation matrix in homogeneous coordinates. */
    static MatrixMN homTranslation(VectorN<2, T> const & v)
    {
      return homTranslation(v[0], v[1]);
    }

    /** Rotate around the given 3D axis (need not be a unit vector) by a given angle. */
    static MatrixMN rotationAxisAngle(VectorN<3, T> const & axis, Real radians)
    {
        VectorN<3, T> uaxis = axis.unit();

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

        MatrixMN m;
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
    static MatrixMN rotationEulerAnglesXYZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0.0, sin_val, cos_val);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      return xmat * (ymat * zmat);
    }

    /**
     * Rotate about the Y axis by the roll angle, then the Z axis by the pitch angle, and finally the X axis by the yaw angle.
     */
    static MatrixMN rotationEulerAnglesXZY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      return xmat * (zmat * ymat);
    }

    /**
     * Rotate about the Z axis by the roll angle, then the X axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static MatrixMN rotationEulerAnglesYXZ(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      return ymat * (xmat * zmat);
    }

    /**
     * Rotate about the X axis by the roll angle, then the Z axis by the pitch angle, and finally the Y axis by the yaw angle.
     */
    static MatrixMN rotationEulerAnglesYZX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val);

      return ymat * (zmat * xmat);
    }

    /**
     * Rotate about the Y axis by the roll angle, then the X axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static MatrixMN rotationEulerAnglesZXY(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      return zmat * (xmat * ymat);
    }

    /**
     * Rotate about the X axis by the roll angle, then the Y axis by the pitch angle, and finally the Z axis by the yaw angle.
     */
    static MatrixMN rotationEulerAnglesZYX(Real yaw_radians, Real pitch_radians, Real roll_radians)
    {
      T cos_val, sin_val;

      cos_val = static_cast<T>(std::cos(yaw_radians));
      sin_val = static_cast<T>(std::sin(yaw_radians));
      MatrixMN zmat(cos_val, -sin_val, 0, sin_val, cos_val, 0, 0, 0, 1);

      cos_val = static_cast<T>(std::cos(pitch_radians));
      sin_val = static_cast<T>(std::sin(pitch_radians));
      MatrixMN ymat(cos_val, 0, sin_val, 0, 1, 0, -sin_val, 0, cos_val);

      cos_val = static_cast<T>(std::cos(roll_radians));
      sin_val = static_cast<T>(std::sin(roll_radians));
      MatrixMN xmat(1, 0, 0, 0, cos_val, -sin_val, 0, sin_val, cos_val);

      return zmat * (ymat * xmat);
    }

    /**
     * Return the matrix corresponding to the rotation from one direction vector to another.
     *
     * @param start_dir The vector to rotate from.
     * @param end_dir The vector to rotate to.
     * @param unitize_dirs If false, the directions will be assumed to have been pre-normalized to unit length before being
     *   passed to this function.
     */
    static MatrixMN rotationArc(VectorN<3, T> const & start_dir, VectorN<3, T> const & end_dir, bool unitize_dirs = true)
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

      VectorN<3, T> u = start_dir, v = end_dir;
      if (unitize_dirs)
      {
        u.unitize();
        v.unitize();
      }

      T d = u.dot(v);
      if (d < -0.9999)
      {
        VectorN<3, T> perp;
        switch (u.maxAbsAxis())
        {
          case 0:           perp = VectorN<3, T>(u.y(),  -u.x(),       0); break;
          case 1:           perp = VectorN<3, T>(u.y(),  -u.x(),       0); break;
          default /* 2 */:  perp = VectorN<3, T>(u.z(),       0,  -u.x());
        }

        return rotationAxisAngle(perp, (Real)Math::pi());
      }

      T s = std::sqrt((1 + d) * 2);
      T recip = 1 / s;
      VectorN<3, T> rcc = recip * u.cross(v);

      return rotationQuat(rcc.x(), rcc.y(), rcc.z(), static_cast<T>(0.5 * s));
    }

  private:
    /**
     * Householder reduction T = Q^t M Q.
     *   Input:
     *     mat, symmetric 3x3 matrix M
     *   Output:
     *     mat, orthogonal matrix Q
     *     diag, diagonal entries of T (preallocated with at least 3 elements)
     *     subd, subdiagonal entries of T (T is symmetric) (preallocated with at least 3 elements)
     */
    void tridiagonal(T * afDiag, T * afSubDiag)
    {
      // From G3D::Matrix3

      T fA = (*this)(0, 0);
      T fB = (*this)(0, 1);
      T fC = (*this)(0, 2);
      T fD = (*this)(1, 1);
      T fE = (*this)(1, 2);
      T fF = (*this)(2, 2);

      afDiag[0] = fA;
      afSubDiag[2] = 0;

      if (Math::fuzzyNe(fC, static_cast<T>(0)))
      {
        T fLength = std::sqrt(fB * fB + fC * fC);
        T fInvLength = 1 / fLength;
        fB *= fInvLength;
        fC *= fInvLength;
        T fQ = 2 * fB * fE + fC * (fF - fD);
        afDiag[1] = fD + fC * fQ;
        afDiag[2] = fF - fC * fQ;
        afSubDiag[0] = fLength;
        afSubDiag[1] = fE - fB * fQ;
        (*this)(0, 0) =   1;
        (*this)(0, 1) =   0;
        (*this)(0, 2) =   0;
        (*this)(1, 0) =   0;
        (*this)(1, 1) =  fB;
        (*this)(1, 2) =  fC;
        (*this)(2, 0) =   0;
        (*this)(2, 1) =  fC;
        (*this)(2, 2) = -fB;
      }
      else
      {
        afDiag[1] = fD;
        afDiag[2] = fF;
        afSubDiag[0] = fB;
        afSubDiag[1] = fE;
        this->makeIdentity();
      }
    }

    /** QL iteration with implicit shifting to reduce matrix from tridiagonal to diagonal. */
    bool qlAlgorithm(T * afDiag, T * afSubDiag)
    {
      // From G3D::Matrix3

      for (int i0 = 0; i0 < 3; i0++)
      {
        int const iMaxIter = 32;
        int iIter;

        for (iIter = 0; iIter < iMaxIter; iIter++)
        {
          int i1;

          for (i1 = i0; i1 <= 1; i1++)
          {
            T fSum = std::abs(afDiag[i1]) + std::abs(afDiag[i1 + 1]);
            if (std::abs(afSubDiag[i1]) + fSum == fSum)
              break;
          }

          if ( i1 == i0 )
            break;

          T fTmp0 = (afDiag[i0 + 1] - afDiag[i0]) / (2 * afSubDiag[i0]);
          T fTmp1 = std::sqrt(fTmp0 * fTmp0 + 1);

          if (fTmp0 < 0)
            fTmp0 = afDiag[i1] - afDiag[i0] + afSubDiag[i0] / (fTmp0 - fTmp1);
          else
            fTmp0 = afDiag[i1] - afDiag[i0] + afSubDiag[i0] / (fTmp0 + fTmp1);

          T fSin = 1;
          T fCos = 1;
          T fTmp2 = 0;

          for (int i2 = i1 - 1; i2 >= i0; i2--)
          {
            T fTmp3 = fSin * afSubDiag[i2];
            T fTmp4 = fCos * afSubDiag[i2];

            if (std::abs(fTmp3) >= std::abs(fTmp0))
            {
              fCos = fTmp0 / fTmp3;
              fTmp1 = std::sqrt(fCos * fCos + 1);
              afSubDiag[i2 + 1] = fTmp3 * fTmp1;
              fSin = 1 / fTmp1;
              fCos *= fSin;
            }
            else
            {
              fSin = fTmp3 / fTmp0;
              fTmp1 = std::sqrt(fSin * fSin + 1);
              afSubDiag[i2 + 1] = fTmp0 * fTmp1;
              fCos = 1 / fTmp1;
              fSin *= fCos;
            }

            fTmp0 = afDiag[i2 + 1] - fTmp2;
            fTmp1 = (afDiag[i2] - fTmp0) * fSin + 2 * fTmp4 * fCos;
            fTmp2 = fSin * fTmp1;
            afDiag[i2 + 1] = fTmp0 + fTmp2;
            fTmp0 = fCos * fTmp1 - fTmp4;

            for (int iRow = 0; iRow < 3; iRow++)
            {
              fTmp3 = (*this)(iRow, i2 + 1);
              (*this)(iRow, i2 + 1) = fSin * (*this)(iRow, i2) + fCos * fTmp3;
              (*this)(iRow, i2) = fCos * (*this)(iRow, i2) - fSin * fTmp3;
            }
          }

          afDiag[i0] -= fTmp2;
          afSubDiag[i0] = fTmp0;
          afSubDiag[i1] = 0;
        }

        if (iIter == iMaxIter)
        {
          // should not get here under normal circumstances
          return false;
        }
      }

      return true;
    }

    /** Construct a rotation matrix from a quaternion. */
    static MatrixMN rotationQuat(T const & x, T const & y, T const & z, T const & w)
    {
      // From G3D::Matrix3:
      // Implementation from Watt and Watt, pg 362
      // See also http://www.flipcode.com/documents/matrfaq.html#Q54

      VectorN<4, T> q(x, y, z, w);
      q.unitize();

      T xx = 2 * q.x() * q.x();
      T xy = 2 * q.x() * q.y();
      T xz = 2 * q.x() * q.z();
      T xw = 2 * q.x() * q.w();

      T yy = 2 * q.y() * q.y();
      T yz = 2 * q.y() * q.z();
      T yw = 2 * q.y() * q.w();

      T zz = 2 * q.z() * q.z();
      T zw = 2 * q.z() * q.w();

      return MatrixMN(1 - yy - zz,      xy - zw,      xz + yw,
                          xy + zw,  1 - xx - zz,      yz - xw,
                          xz - yw,      yz + xw,  1 - xx - yy);
    }

}; // class MatrixMN<3, 3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API MatrixMN<3, 3, Real>;
#endif

/** The default 3x3 real matrix class. */
typedef MatrixMN<3, 3, Real> Matrix3;

/**
 * Multiply a 2-vector by a homogenous 3x3 matrix (assuming the vector has z = 1) and scale back to non-homogeneous coordinates.
 */
template <typename T>
VectorN<2, T>
operator*(MatrixMN<3, 3, T> const & m, VectorN<2, T> const & v)
{
  VectorN<3, T> result = m * VectorN<3, T>(v, 1);

  if (Math::fuzzyEq(result.z(), 0))
    return VectorN<2, T>(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
  else
    return result.xy() / result.z();
}

} // namespace Thea

#endif
