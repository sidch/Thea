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

#ifndef __Thea_IDenseMatrix_hpp__
#define __Thea_IDenseMatrix_hpp__

#include "IAddressableMatrix.hpp"
#include "MatVec.hpp"
#include "IRowOrColumnMajorMatrix.hpp"
#include <type_traits>

namespace Thea {

/**
 * Interface for a 2D dense matrix, assumed to be packed with no gaps or padding in a contiguous memory block. Useful for
 * passing matrices across shared library boundaries.
 */
template <typename T>
class /* THEA_API */ IDenseMatrix : public virtual IAddressableMatrix<T>, public virtual IRowOrColumnMajorMatrix<T>
{
  public:
    THEA_DECL_SMART_POINTERS(IDenseMatrix)

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T const * THEA_ICALL data() const = 0;

    /** Get a pointer to the beginning of the matrix's data block. */
    virtual T * THEA_ICALL data() = 0;

    /**
     * Set all elements of the matrix to a given value.
     *
     * @warning This is safe to call across shared library boundaries only if T is a POD type.
     */
    virtual void THEA_ICALL fill(T value) = 0;

}; // class IDenseMatrix

namespace Math {

/**
 * Convenience function for interpreting raw buffers wrapped in IDenseMatrix objects as Eigen::Map objects. It avoids
 * having to pass the matrix dimensions separately as required by the Eigen::Map constructors for non-fixed-size matrices. Note
 * that the function signature is equivalent to:
 *
 * \code
 *   Eigen::Map<MatrixT const> mapTo(IDenseMatrix<ScalarT> const & m);
 * \endcode
 *
 * for const-qualified matrix types, and
 *
 * \code
 *   Eigen::Map<MatrixT> mapTo(IDenseMatrix<ScalarT> & m);
 * \endcode
 *
 * for non-const types.
 *
 * E.g.:
 * \code
 *   IDenseMatrix<Real> const * d = <... get a matrix e.g. from across a DLL boundary ...>
 *   auto m = Math::mapTo<MatrixX<> const>(*d);
 *   ... treat m as a normal Eigen dynamic-size matrix ...
 * \endcode
 */
template < typename MatrixT, typename IDenseMatrixT,
           typename std::enable_if< std::is_base_of< IDenseMatrix<typename IDenseMatrixT::value_type>,
                                                 IDenseMatrixT >::value
                                 && std::is_same<typename IDenseMatrixT::value_type,
                                                 typename MatrixT::value_type>::value
                                 && std::is_const<IDenseMatrixT>::value == std::is_const<MatrixT>::value,
                                    int >::type = 0 >
Eigen::Map<MatrixT>  // this should ensure only valid Eigen types will match this function
mapTo(IDenseMatrixT & m)
{
  alwaysAssertM(MatrixT::IsVectorAtCompileTime || MatrixT::IsRowMajor == m.isRowMajor(),
                "mapTo: IDenseMatrix layout does not match target Eigen::Map matrix type");
  return Eigen::Map<MatrixT>(m.data(), m.rows(), m.cols());   // Eigen should automatically do dimension checks
};

} // namespace Math

} // namespace Thea

#endif
