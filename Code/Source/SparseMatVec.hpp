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

#ifndef __Thea_SparseMatVec_hpp__
#define __Thea_SparseMatVec_hpp__

#include "Common.hpp"
#include "MatrixFormat.hpp"
#include <Eigen/SparseCore>

namespace Thea {

#ifdef THEA_SPARSE_ROW_MAJOR
  int const DEFAULT_SPARSE_MATRIX_LAYOUT = MatrixLayout::ROW_MAJOR;
#else  // nothing provided or THEA_SPARSE_COLUMN_MAJOR defined
  int const DEFAULT_SPARSE_MATRIX_LAYOUT = MatrixLayout::COLUMN_MAJOR;
#endif

/**
 * General 2D sparse matrix template, alias for <code>Eigen::SparseMatrix</code> with a custom default layout (row or column
 * major).
 */
template <typename T = Real,
          int Options = DEFAULT_SPARSE_MATRIX_LAYOUT,
          typename StorageIndex = int>
using SparseMatrix = Eigen::SparseMatrix<T, Options, StorageIndex>;

/** 2D sparse column matrix template, alias for <code>Eigen::SparseMatrix</code> with <code>Eigen::ColMajor</code>. */
template <typename T = Real,
          typename StorageIndex = int>
using SparseColumnMatrix = Eigen::SparseMatrix<T, MatrixLayout::COLUMN_MAJOR, StorageIndex>;

/** 2D sparse row matrix template, alias for <code>Eigen::SparseMatrix</code> with <code>Eigen::RowMajor</code>. */
template <typename T = Real,
          typename StorageIndex = int>
using SparseRowMatrix = Eigen::SparseMatrix<T, MatrixLayout::ROW_MAJOR, StorageIndex>;

/**
 * General 1D sparse vector template, alias for <code>Eigen::SparseVector</code> with a custom default layout (row or column
 * major).
 */
template <typename T = Real,
          int Options = DEFAULT_SPARSE_MATRIX_LAYOUT,
          typename StorageIndex = int>
using SparseVector = Eigen::SparseVector<T, Options, StorageIndex>;

/** 1D sparse column vector template, alias for <code>Eigen::SparseVector</code> with <code>Eigen::ColMajor</code>. */
template <typename T = Real,
          typename StorageIndex = int>
using SparseColumnVector = Eigen::SparseVector<T, MatrixLayout::COLUMN_MAJOR, StorageIndex>;

/** 1D sparse row vector template, alias for <code>Eigen::SparseVector</code> with <code>Eigen::RowMajor</code>. */
template <typename T = Real,
          typename StorageIndex = int>
using SparseRowVector = Eigen::SparseVector<T, MatrixLayout::ROW_MAJOR, StorageIndex>;

} // namespace Thea

#include "MatrixUtil.hpp"

#endif
