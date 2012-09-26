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

#ifndef __Thea_Matrix2_hpp__
#define __Thea_Matrix2_hpp__

#include "MatrixMN.hpp"
#include <cmath>

namespace Thea {

/**
 * Square 2 x 2 matrices on a field T. The matrices are stored internally in row-major form, so row-major access is recommended.
 */
template <typename T>
class /* THEA_API */ MatrixMN<2, 2, T> : public Internal::SquareMatrixN<2, T>
{
  private:
    typedef Internal::SquareMatrixN<2, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    MatrixMN() {}

    /** Initialize all components to a single value. */
    explicit MatrixMN(T const & fill_value) : BaseT(fill_value) {}

    /** Initialize all 4 components of the matrix. */
    MatrixMN(T const & m00, T const & m01,
             T const & m10, T const & m11)
    {
      (*this)(0, 0) = m00; (*this)(0, 1) = m01;
      (*this)(1, 0) = m10; (*this)(1, 1) = m11;
    }

    /** Matrix to rotate about the origin by an angle (in radians). */
    static MatrixMN rotation(Real radians)
    {
      T s = std::sin(radians);
      T c = std::cos(radians);
      return MatrixMN(c, -s, s, c);
    }

    /** Get the matrix of cofactors. */
    MatrixMN cofactor() const
    {
      return adjoint().transpose();
    }

    /** Get the adjoint of the matrix. */
    MatrixMN adjoint() const
    {
      return MatrixMN((*this)(1, 1), -(*this)(0, 1),
                     -(*this)(1, 0),  (*this)(0, 0));
    }

    /** Get the determinant of the matrix. */
    T determinant() const
    {
      return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
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
      // Invert a 2x2 matrix using cofactors.

      T det = determinant();
      if (std::abs(det) <= tolerance)
        return false;

      T inv_det = 1 / det;
      *this = MatrixMN(inv_det * (*this)(1, 1), -inv_det * (*this)(0, 1),
                      -inv_det * (*this)(1, 0),  inv_det * (*this)(0, 0));

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
      if (!result.invert())
        throw Error("MatrixMN<2, 2, T>: Could not invert matrix " + this->toString() + " with given tolerance");

      return result;
    }

}; // class MatrixMN<2, 2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API MatrixMN<2, 2, Real>;
#endif

/** The default 2x2 real matrix class. */
typedef MatrixMN<2, 2, Real> Matrix2;

} // namespace Thea

#endif
