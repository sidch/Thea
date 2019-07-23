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
#include <type_traits>

namespace Thea {

/**
 * Abstract base interface for a 2D dense matrix, assumed to be packed with no gaps or padding in a contiguous memory block.
 * Useful for passing matrices across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ AbstractDenseMatrix : public virtual AbstractAddressableMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(AbstractDenseMatrix)

    /** Is the matrix stored in row-major format? */
    virtual int8 isRowMajor() const = 0;

    /** Is the matrix stored in column-major format? */
    virtual int8 isColumnMajor() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T const * data() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T * data() = 0;

    /** Set all elements of the matrix to a given value. */
    virtual void fill(T const & value) = 0;

}; // class AbstractDenseMatrix

namespace Math {

/**
 * Convenience function for interpreting raw buffers wrapped in AbstractDenseMatrix objects as Eigen::Map objects. It avoids
 * having to pass the matrix dimensions separately as required by the Eigen::Map constructors for non-fixed-size matrices. Note
 * that the function signature is equivalent to:
 *
 * \code
 *   Eigen::Map<MatrixT const> mapTo(AbstractDenseMatrix<ScalarT> const & m);
 * \endcode
 *
 * for const-qualified matrix types, and
 *
 * \code
 *   Eigen::Map<MatrixT> mapTo(AbstractDenseMatrix<ScalarT> & m);
 * \endcode
 *
 * for non-const types.
 *
 * E.g.:
 * \code
 *   AbstractDenseMatrix<Real> const * d = <... get a matrix e.g. from across a DLL boundary ...>
 *   MatrixXConstMap<> m = Math::mapTo<MatrixX<> const>(*d);
 *   ... treat m as a normal Eigen dynamic-size matrix ...
 * \endcode
 */
template <typename MatrixT, typename AbstractDenseMatrixT>
Eigen::Map<MatrixT>  // this should ensure only valid Eigen types will match this function
mapTo(AbstractDenseMatrixT & m,
      typename std::enable_if< std::is_base_of< AbstractDenseMatrix<typename AbstractDenseMatrixT::value_type>,
                                                AbstractDenseMatrixT >::value
                            && std::is_same<typename AbstractDenseMatrixT::value_type,
                                            typename MatrixT::value_type>::value
                            && std::is_const<AbstractDenseMatrixT>::value == std::is_const<MatrixT>::value
                             >::type * dummy = 0)
{
  alwaysAssertM(MatrixT::IsVectorAtCompileTime || MatrixT::IsRowMajor == m.isRowMajor(),
                "mapTo: AbstractDenseMatrix layout does not match target Eigen::Map matrix type");
  return Eigen::Map<MatrixT>(m.data(), m.rows(), m.cols());   // Eigen should automatically do dimension checks
};

} // namespace Math

} // namespace Thea

#endif
