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

#ifndef __Thea_MatVec_hpp__
#define __Thea_MatVec_hpp__

#include "Common.hpp"
#include "MatrixFormat.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>  // for some weird reason Eigen defines cross() here
#include <complex>

namespace Thea {

#ifdef THEA_ROW_MAJOR
  int const DEFAULT_MATRIX_LAYOUT = MatrixLayout::ROW_MAJOR;
#else  // nothing provided or THEA_COLUMN_MAJOR defined
  int const DEFAULT_MATRIX_LAYOUT = MatrixLayout::COLUMN_MAJOR;
#endif

// Typedef common instantiations. These are fully instantiated and have no further template parameters.
#define THEA_DECL_MATRIX_TYPEDEFS(scalar, suffix) \
    typedef Eigen::Matrix< scalar, 2,              2             , DEFAULT_MATRIX_LAYOUT      >   Matrix2    ## suffix; \
    typedef Eigen::Matrix< scalar, 2,              Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT      >   Matrix2X   ## suffix; \
    typedef Eigen::Matrix< scalar, 3,              3             , DEFAULT_MATRIX_LAYOUT      >   Matrix3    ## suffix; \
    typedef Eigen::Matrix< scalar, 3,              Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT      >   Matrix3X   ## suffix; \
    typedef Eigen::Matrix< scalar, 4,              4             , DEFAULT_MATRIX_LAYOUT      >   Matrix4    ## suffix; \
    typedef Eigen::Matrix< scalar, 4,              Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT      >   Matrix4X   ## suffix; \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 2             , DEFAULT_MATRIX_LAYOUT      >   MatrixX2   ## suffix; \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 3             , DEFAULT_MATRIX_LAYOUT      >   MatrixX3   ## suffix; \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 4             , DEFAULT_MATRIX_LAYOUT      >   MatrixX4   ## suffix; \
    typedef Eigen::Matrix< scalar, 2,              1             , MatrixLayout::COLUMN_MAJOR >   Vector2    ## suffix; \
    typedef Eigen::Matrix< scalar, 3,              1             , MatrixLayout::COLUMN_MAJOR >   Vector3    ## suffix; \
    typedef Eigen::Matrix< scalar, 4,              1             , MatrixLayout::COLUMN_MAJOR >   Vector4    ## suffix; \
    typedef Eigen::Matrix< scalar, 1,              2             , MatrixLayout::ROW_MAJOR    >   RowVector2 ## suffix; \
    typedef Eigen::Matrix< scalar, 1,              3             , MatrixLayout::ROW_MAJOR    >   RowVector3 ## suffix; \
    typedef Eigen::Matrix< scalar, 1,              4             , MatrixLayout::ROW_MAJOR    >   RowVector4 ## suffix;

THEA_DECL_MATRIX_TYPEDEFS(Real, )
THEA_DECL_MATRIX_TYPEDEFS(float, f)
THEA_DECL_MATRIX_TYPEDEFS(double, d)
THEA_DECL_MATRIX_TYPEDEFS(std::complex<float>, cf)
THEA_DECL_MATRIX_TYPEDEFS(std::complex<double>, cd)
THEA_DECL_MATRIX_TYPEDEFS(int, i)

#undef THEA_DECL_MATRIX_TYPEDEFS

// Typedef additional instantiations of fully resizable matrices, NOT including the "plain X" versions which will be declared
// below with (fully optional) template arguments to avoid having to repeatedly type
// <code>Matrix<Eigen::Dynamic, Eigen::Dynamic, T></code> in templated code.
#define THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(scalar, suffix) \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT      >   MatrixX    ## suffix; \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 1             , MatrixLayout::COLUMN_MAJOR >   VectorX    ## suffix; \
    typedef Eigen::Matrix< scalar, 1,              Eigen::Dynamic, MatrixLayout::ROW_MAJOR    >   RowVectorX ## suffix;

THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(float, f)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(double, d)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(std::complex<float>, cf)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(std::complex<double>, cd)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(int, i)

#undef THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS

/**
 * General 2D dense matrix template, alias for <code>Eigen::Matrix</code> with a custom default layout (row or column major).
 */
template <int Rows, int Cols, typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Rows,
          int MaxColsAtCompileTime = Cols>
using Matrix = Eigen::Matrix<T, Rows, Cols, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/**
 * General 2D dense matrix template with dynamic resizing, alias for <code>Eigen::Matrix</code> with Eigen::Dynamic and a custom
 * default layout (row or column major).
 */
template <typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Eigen::Dynamic,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/** General 1D dense column vector template, alias for <code>Eigen::Matrix<T, Size, 1,...></code>. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Size,
          int MaxColsAtCompileTime = 1>
using Vector = Eigen::Matrix<T, Size, 1, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/**
 * General 1D dense column vector template with dynamic resizing, alias for <code>Eigen::Matrix<T, Eigen::Dynamic, 1,...></code>
 * with a custom default layout (row or column major).
 */
template <typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Eigen::Dynamic,
          int MaxColsAtCompileTime = 1>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/** General 1D dense row vector template, alias for <code>Eigen::Matrix<T, 1, Size,...></code>. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxRowsAtCompileTime = 1,
          int MaxColsAtCompileTime = Size>
using RowVector = Eigen::Matrix<T, 1, Size, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/**
 * General 1D dense row vector template with dynamic resizing, alias for <code>Eigen::Matrix<T, 1, Eigen::Dynamic,...></code>
 * with a custom default layout (row or column major).
 */
template <typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxRowsAtCompileTime = 1,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using RowVectorX = Eigen::Matrix<T, 1, Eigen::Dynamic, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

} // namespace Thea

#include "MatrixUtil.hpp"

#endif
