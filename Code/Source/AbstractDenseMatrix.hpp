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
// First version: 2019
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
template < typename MatrixT, typename AbstractDenseMatrixT,
           typename std::enable_if< std::is_base_of< AbstractDenseMatrix<typename AbstractDenseMatrixT::value_type>,
                                                 AbstractDenseMatrixT >::value
                                 && std::is_same<typename AbstractDenseMatrixT::value_type,
                                                 typename MatrixT::value_type>::value
                                 && std::is_const<AbstractDenseMatrixT>::value == std::is_const<MatrixT>::value,
                                    int >::type = 0 >
Eigen::Map<MatrixT>  // this should ensure only valid Eigen types will match this function
mapTo(AbstractDenseMatrixT & m)
{
  alwaysAssertM(MatrixT::IsVectorAtCompileTime || MatrixT::IsRowMajor == m.isRowMajor(),
                "mapTo: AbstractDenseMatrix layout does not match target Eigen::Map matrix type");
  return Eigen::Map<MatrixT>(m.data(), m.rows(), m.cols());   // Eigen should automatically do dimension checks
};

} // namespace Math

} // namespace Thea

#endif
