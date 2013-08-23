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

#ifndef __Thea_MatrixUtil_hpp__
#define __Thea_MatrixUtil_hpp__

#include "Common.hpp"
#include "CompressedSparseMatrix.hpp"
#include "Matrix.hpp"
#include "MatrixFormat.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace Thea {

/** Utility functions for matrices. */
namespace MatrixUtil {

namespace MatrixUtilInternal {

// isSquare helper classes

// Default
template <typename MatrixT, typename Enable = void>
struct IsSquareImpl { static bool isSquare(MatrixT const & m) { return m.size1() == m.size2(); } };

// Check if a standard matrix is square
template <typename MatrixT>
struct IsSquareImpl< MatrixT, typename boost::enable_if< boost::is_base_of< BasicMatrix<typename MatrixT::Value>,
                                                                            MatrixT > >::type >
{ static bool isSquare(MatrixT const & m) { return m.isSquare(); } };

// getFormat helper classes

// Default
template <typename MatrixT, typename Enable = void>
struct GetFormatImpl { static MatrixFormat getFormat(MatrixT const & m) { return MatrixFormat::UNKNOWN; } };

// Get the format of a dense matrix. */
template <typename MatrixT>
struct GetFormatImpl< MatrixT, typename boost::enable_if< boost::is_base_of< Matrix<typename MatrixT::Value, MatrixT::Layout>,
                                                                             MatrixT > >::type >
{
  static MatrixFormat getFormat(MatrixT const & m)
  { return m.getLayout() == MatrixLayout::ROW_MAJOR ? MatrixFormat::DENSE_ROW_MAJOR : MatrixFormat::DENSE_COLUMN_MAJOR; }
};

// Get the format of a sparse matrix. */
template <typename MatrixT>
struct GetFormatImpl< MatrixT, typename boost::enable_if< boost::is_base_of< CompressedSparseMatrix<typename MatrixT::Value,
                                                                                                    MatrixT::Layout,
                                                                                                    typename MatrixT::Index2D,
                                                                                                    typename MatrixT::Index1D>,
                                                                             MatrixT > >::type >
{
  static MatrixFormat getFormat(MatrixT const & m)
  { return m.getLayout() == MatrixLayout::ROW_MAJOR ? MatrixFormat::SPARSE_ROW_MAJOR : MatrixFormat::SPARSE_COLUMN_MAJOR; }
};

} // namespace MatrixUtilInternal

/** Test if a matrix is square. */
template <typename MatrixT>
bool
isSquare(MatrixT const & m)
{
  return MatrixUtilInternal::IsSquareImpl<MatrixT>::isSquare(m);
}

/** Get the format of a matrix. */
template <typename MatrixT>
MatrixFormat
getFormat(MatrixT const & m)
{
  return MatrixUtilInternal::GetFormatImpl<MatrixT>::getFormat(m);
}

} // namespace MatrixUtil

} // namespace Thea

#endif
