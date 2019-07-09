//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2019, Siddhartha Chaudhuri
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

#ifndef __Thea_AbstractDenseMatrix_hpp__
#define __Thea_AbstractDenseMatrix_hpp__

#include "AbstractAddressableMatrix.hpp"
#include "MatVec.hpp"

namespace Thea {

/**
 * Abstract base interface for a 2D dense matrix, assumed to be packed with no gaps or padding in a contiguous memory block.
 * Useful for passing matrices across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ AbstractDenseMatrix : public virtual AbstractAddressableMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(AbstractDenseMatrix, std::shared_ptr, std::weak_ptr)

    /** Is the matrix stored in row-major format? */
    virtual bool isRowMajor() const = 0;

    /** Is the matrix stored in column-major format? */
    virtual bool isColumnMajor() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T const * data() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T * data() = 0;

    /** Set all elements of the matrix to a given value. */
    virtual void fill(T const & value) = 0;

}; // class AbstractDenseMatrix

//=============================================================================================================================
// Define conversion functions for interpreting AbstractDenseMatrix objects as Eigen::Map objects.
//
// E.g.:
//     AbstractDenseMatrix<Real> const * d = <... get a matrix e.g. from across a DLL boundary ...>
//     MatrixXConstMap m = Math::toMatrixXConstMap(*d);
//     ... treat m as a normal Eigen dynamic-size matrix ...
//=============================================================================================================================

namespace Math {

#define THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER1(scalar, matrix_type, const_flag)          \
  inline matrix_type ## Map to ## matrix_type(AbstractDenseMatrix<scalar> const_flag & m)  \
  { return matrix_type ## Map(m.data(), m.rows(), m.cols()); }  // Eigen should automatically do dimension checks

#define THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, matrix_type)            \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER1(scalar, matrix_type, )                \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER1(scalar, matrix_type ## Const, const)

#define THEA_DEFINE_MATRIX_MAP_FUNCTIONS(scalar, suffix)                  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix2    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix3    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix4    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Vector2    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Vector3    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Vector4    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, RowVector2 ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, RowVector3 ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, RowVector4 ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix2X   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix3X   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, Matrix4X   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, MatrixX2   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, MatrixX3   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, MatrixX4   ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, MatrixX    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, VectorX    ## suffix)  \
  THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2(scalar, RowVectorX ## suffix)

THEA_DEFINE_MATRIX_MAP_FUNCTIONS(Real, )
THEA_DEFINE_MATRIX_MAP_FUNCTIONS(float, f)
THEA_DEFINE_MATRIX_MAP_FUNCTIONS(double, d)
THEA_DEFINE_MATRIX_MAP_FUNCTIONS(std::complex<float>, cf)
THEA_DEFINE_MATRIX_MAP_FUNCTIONS(std::complex<double>, cd)
THEA_DEFINE_MATRIX_MAP_FUNCTIONS(int, i)

#undef THEA_DEFINE_MATRIX_MAP_FUNCTIONS
#undef THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER2
#undef THEA_DEFINE_MATRIX_MAP_FUNCTIONS_HELPER1

} // namespace Math

} // namespace Thea

#endif
