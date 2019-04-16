//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, * except for portions
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
