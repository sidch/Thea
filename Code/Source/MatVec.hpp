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

#ifndef __Thea_MatVec_hpp__
#define __Thea_MatVec_hpp__

#include "Common.hpp"
#include "MatrixFormat.hpp"

// Resolve a conflict between X.h (which #defines Success) and Eigen (which has an enum called Success). Fortunately both
// headers currently assign the same value to Success (0) so the #undef should be all that's needed.
// https://stackoverflow.com/questions/22400905/eigen-and-cimg-compatibility-issues
#ifdef Success
#  undef Success
#endif
#include <Eigen/Core>
#include <Eigen/Geometry>  // for Quaternion, and because for some weird reason Eigen defines cross() here

#include <complex>

namespace Thea {

#ifdef THEA_ROW_MAJOR
  int const DEFAULT_MATRIX_LAYOUT = MatrixLayout::ROW_MAJOR;
#else  // nothing provided or THEA_COLUMN_MAJOR defined
  int const DEFAULT_MATRIX_LAYOUT = MatrixLayout::COLUMN_MAJOR;
#endif

// Typedef common instantiations. These are fully instantiated and have no further template parameters.
// Alignment is DISABLED for fixed-size vectorizable types (Vector4f, Matrix4f, Vector2d etc). These require special handling
// for every class which directly or indirectly has such a member, and could lead to unexpected bugs and considerably less ease
// of use.
//
// See:
// - http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
// - https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
// - http://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
// - https://stackoverflow.com/q/41087043
//
#define THEA_DECL_MATRIX_TYPEDEFS(scalar, suffix)                                                                             \
    typedef Eigen::Matrix< scalar, 2, 2,              DEFAULT_MATRIX_LAYOUT      | Eigen::DontAlign >  Matrix2    ## suffix;  \
    typedef Eigen::Matrix< scalar, 3, 3,              DEFAULT_MATRIX_LAYOUT      | Eigen::DontAlign >  Matrix3    ## suffix;  \
    typedef Eigen::Matrix< scalar, 4, 4,              DEFAULT_MATRIX_LAYOUT      | Eigen::DontAlign >  Matrix4    ## suffix;  \
    typedef Eigen::Matrix< scalar, 2, 1,              MatrixLayout::COLUMN_MAJOR | Eigen::DontAlign >  Vector2    ## suffix;  \
    typedef Eigen::Matrix< scalar, 3, 1,              MatrixLayout::COLUMN_MAJOR | Eigen::DontAlign >  Vector3    ## suffix;  \
    typedef Eigen::Matrix< scalar, 4, 1,              MatrixLayout::COLUMN_MAJOR | Eigen::DontAlign >  Vector4    ## suffix;  \
    typedef Eigen::Matrix< scalar, 1, 2,              MatrixLayout::ROW_MAJOR    | Eigen::DontAlign >  RowVector2 ## suffix;  \
    typedef Eigen::Matrix< scalar, 1, 3,              MatrixLayout::ROW_MAJOR    | Eigen::DontAlign >  RowVector3 ## suffix;  \
    typedef Eigen::Matrix< scalar, 1, 4,              MatrixLayout::ROW_MAJOR    | Eigen::DontAlign >  RowVector4 ## suffix;  \
    typedef Eigen::Matrix< scalar, 2, Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT                         >  Matrix2X   ## suffix;  \
    typedef Eigen::Matrix< scalar, 3, Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT                         >  Matrix3X   ## suffix;  \
    typedef Eigen::Matrix< scalar, 4, Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT                         >  Matrix4X   ## suffix;  \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 2, DEFAULT_MATRIX_LAYOUT                         >  MatrixX2   ## suffix;  \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 3, DEFAULT_MATRIX_LAYOUT                         >  MatrixX3   ## suffix;  \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 4, DEFAULT_MATRIX_LAYOUT                         >  MatrixX4   ## suffix;

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
#define THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(scalar, suffix)                                                                   \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, Eigen::Dynamic, DEFAULT_MATRIX_LAYOUT      >  MatrixX    ## suffix;        \
    typedef Eigen::Matrix< scalar, Eigen::Dynamic, 1             , MatrixLayout::COLUMN_MAJOR >  VectorX    ## suffix;        \
    typedef Eigen::Matrix< scalar, 1,              Eigen::Dynamic, MatrixLayout::ROW_MAJOR    >  RowVectorX ## suffix;

THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(float, f)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(double, d)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(std::complex<float>, cf)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(std::complex<double>, cd)
THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS(int, i)

#undef THEA_DECL_RESIZABLE_MATRIX_TYPEDEFS

/**
 * General 2D dense matrix template, alias for <code>Eigen::Matrix</code> with a custom default layout (row or column major) and
 * alignment preference.
 *
 * This alias currently <b>disables alignment</b> for all fixed-size types: fixed-size vectorizable types (Vector4f, Matrix4f,
 * Vector2d etc) require special handling for every class which directly or indirectly has such a member, and could lead to
 * unexpected bugs and considerably less ease of use.
 *
 * @see
 *   - http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
 *   - https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
 *   - http://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
 *   - https://stackoverflow.com/q/41087043
 */
template <int Rows, int Cols, typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Rows,
          int MaxColsAtCompileTime = Cols>
using Matrix = Eigen::Matrix<T, Rows, Cols,
                             Options | ((Options & Eigen::DontAlign) == 0 && (Rows == Eigen::Dynamic || Cols == Eigen::Dynamic)
                                      ? Eigen::AutoAlign : Eigen::DontAlign),
                             MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/**
 * General 2D dense matrix template with dynamic resizing, alias for <code>Eigen::Matrix</code> with Eigen::Dynamic and a custom
 * default layout (row or column major).
 */
template <typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Eigen::Dynamic,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;

/**
 * General 1D dense column vector template, alias for <code>Eigen::Matrix<T, Size, 1,...></code>, with a custom alignment
 * preference.
 *
 * This alias currently <b>disables alignment</b> for all fixed-size types: fixed-size vectorizable types (Vector4f, Vector2d
 * etc) require special handling for every class which directly or indirectly has such a member, and could lead to unexpected
 * bugs and considerably less ease of use.
 *
 * @see
 *   - http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
 *   - https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
 *   - http://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
 *   - https://stackoverflow.com/q/41087043
 */
template <int Size, typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Size>
using Vector = Eigen::Matrix<T, Size, 1,
                             Options | ((Options & Eigen::DontAlign) == 0 && Size == Eigen::Dynamic
                                      ? Eigen::AutoAlign : Eigen::DontAlign),
                             MaxRowsAtCompileTime, 1>;

/**
 * General 1D dense column vector template with dynamic resizing, alias for <code>Eigen::Matrix<T, Eigen::Dynamic, 1,...></code>
 * with a custom default layout (row or column major).
 */
template <typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Eigen::Dynamic>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1, Options, MaxRowsAtCompileTime, 1>;

/**
 * General 1D dense row vector template, alias for <code>Eigen::Matrix<T, 1, Size,...></code>, with a custom alignment
 * preference.
 *
 * This alias currently <b>disables alignment</b> for all fixed-size types: fixed-size vectorizable types (RowVector4f,
 * RowVector2d etc) require special handling for every class which directly or indirectly has such a member, and could lead to
 * unexpected bugs and considerably less ease of use.
 *
 * @see
 *   - http://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
 *   - https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
 *   - http://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
 *   - https://stackoverflow.com/q/41087043
 */
template <int Size, typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Size>
using RowVector = Eigen::Matrix<T, 1, Size,
                                Options | ((Options & Eigen::DontAlign) == 0 && Size == Eigen::Dynamic
                                         ? Eigen::AutoAlign : Eigen::DontAlign),
                                1, MaxColsAtCompileTime>;

/**
 * General 1D dense row vector template with dynamic resizing, alias for <code>Eigen::Matrix<T, 1, Eigen::Dynamic,...></code>
 * with a custom default layout (row or column major).
 */
template <typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using RowVectorX = Eigen::Matrix<T, 1, Eigen::Dynamic, Options, 1, MaxColsAtCompileTime>;

//=============================================================================================================================
// Typedef Eigen::Map wrappers for interpreting raw data as common Eigen types.
//=============================================================================================================================

#define THEA_DECL_MATRIX_MAP_TYPEDEFS(suffix)                                 \
    typedef Eigen::Map< Matrix2    ## suffix >  Matrix2    ## suffix ## Map;  \
    typedef Eigen::Map< Matrix3    ## suffix >  Matrix3    ## suffix ## Map;  \
    typedef Eigen::Map< Matrix4    ## suffix >  Matrix4    ## suffix ## Map;  \
    typedef Eigen::Map< Vector2    ## suffix >  Vector2    ## suffix ## Map;  \
    typedef Eigen::Map< Vector3    ## suffix >  Vector3    ## suffix ## Map;  \
    typedef Eigen::Map< Vector4    ## suffix >  Vector4    ## suffix ## Map;  \
    typedef Eigen::Map< RowVector2 ## suffix >  RowVector2 ## suffix ## Map;  \
    typedef Eigen::Map< RowVector3 ## suffix >  RowVector3 ## suffix ## Map;  \
    typedef Eigen::Map< RowVector4 ## suffix >  RowVector4 ## suffix ## Map;  \
    typedef Eigen::Map< Matrix2X   ## suffix >  Matrix2X   ## suffix ## Map;  \
    typedef Eigen::Map< Matrix3X   ## suffix >  Matrix3X   ## suffix ## Map;  \
    typedef Eigen::Map< Matrix4X   ## suffix >  Matrix4X   ## suffix ## Map;  \
    typedef Eigen::Map< MatrixX2   ## suffix >  MatrixX2   ## suffix ## Map;  \
    typedef Eigen::Map< MatrixX3   ## suffix >  MatrixX3   ## suffix ## Map;  \
    typedef Eigen::Map< MatrixX4   ## suffix >  MatrixX4   ## suffix ## Map;  \
    \
    typedef Eigen::Map< Matrix2    ## suffix const >  Matrix2    ## suffix ## ConstMap;  \
    typedef Eigen::Map< Matrix3    ## suffix const >  Matrix3    ## suffix ## ConstMap;  \
    typedef Eigen::Map< Matrix4    ## suffix const >  Matrix4    ## suffix ## ConstMap;  \
    typedef Eigen::Map< Vector2    ## suffix const >  Vector2    ## suffix ## ConstMap;  \
    typedef Eigen::Map< Vector3    ## suffix const >  Vector3    ## suffix ## ConstMap;  \
    typedef Eigen::Map< Vector4    ## suffix const >  Vector4    ## suffix ## ConstMap;  \
    typedef Eigen::Map< RowVector2 ## suffix const >  RowVector2 ## suffix ## ConstMap;  \
    typedef Eigen::Map< RowVector3 ## suffix const >  RowVector3 ## suffix ## ConstMap;  \
    typedef Eigen::Map< RowVector4 ## suffix const >  RowVector4 ## suffix ## ConstMap;  \
    typedef Eigen::Map< Matrix2X   ## suffix const >  Matrix2X   ## suffix ## ConstMap;  \
    typedef Eigen::Map< Matrix3X   ## suffix const >  Matrix3X   ## suffix ## ConstMap;  \
    typedef Eigen::Map< Matrix4X   ## suffix const >  Matrix4X   ## suffix ## ConstMap;  \
    typedef Eigen::Map< MatrixX2   ## suffix const >  MatrixX2   ## suffix ## ConstMap;  \
    typedef Eigen::Map< MatrixX3   ## suffix const >  MatrixX3   ## suffix ## ConstMap;  \
    typedef Eigen::Map< MatrixX4   ## suffix const >  MatrixX4   ## suffix ## ConstMap;

THEA_DECL_MATRIX_MAP_TYPEDEFS()
THEA_DECL_MATRIX_MAP_TYPEDEFS(f)
THEA_DECL_MATRIX_MAP_TYPEDEFS(d)
THEA_DECL_MATRIX_MAP_TYPEDEFS(cf)
THEA_DECL_MATRIX_MAP_TYPEDEFS(cd)
THEA_DECL_MATRIX_MAP_TYPEDEFS(i)

#undef THEA_DECL_MATRIX_MAP_TYPEDEFS

#define THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(suffix)                       \
    typedef Eigen::Map< MatrixX    ## suffix >  MatrixX    ## suffix ## Map;  \
    typedef Eigen::Map< VectorX    ## suffix >  VectorX    ## suffix ## Map;  \
    typedef Eigen::Map< RowVectorX ## suffix >  RowVectorX ## suffix ## Map;  \
    \
    typedef Eigen::Map< MatrixX    ## suffix const >  MatrixX    ## suffix ## ConstMap;  \
    typedef Eigen::Map< VectorX    ## suffix const >  VectorX    ## suffix ## ConstMap;  \
    typedef Eigen::Map< RowVectorX ## suffix const >  RowVectorX ## suffix ## ConstMap;

THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(f)
THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(d)
THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(cf)
THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(cd)
THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS(i)

#undef THEA_DECL_RESIZABLE_MATRIX_MAP_TYPEDEFS

/** Alias for Eigen::Map< Matrix<...> >. */
template <int Rows, int Cols, typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Rows,
          int MaxColsAtCompileTime = Cols>
using MatrixMap = Eigen::Map< Matrix<Rows, Cols, T, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> >;

/** Alias for Eigen::Map< Matrix<...> const >. */
template <int Rows, int Cols, typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Rows,
          int MaxColsAtCompileTime = Cols>
using MatrixConstMap = Eigen::Map< Matrix<Rows, Cols, T, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> const >;

/** Alias for Eigen::Map< Vector<...> >. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Size>
using VectorMap = Eigen::Map< Vector<Size, T, Options, MaxRowsAtCompileTime> >;

/** Alias for Eigen::Map< Vector<...> const >. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Size>
using VectorConstMap = Eigen::Map< Vector<Size, T, Options, MaxRowsAtCompileTime> const >;

/** Alias for Eigen::Map< RowVector<...> >. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Size>
using RowVectorMap = Eigen::Map< RowVector<Size, T, Options, MaxColsAtCompileTime> >;

/** Alias for Eigen::Map< RowVector<...> const >. */
template <int Size, typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Size>
using RowVectorConstMap = Eigen::Map< RowVector<Size, T, Options, MaxColsAtCompileTime> const >;

/** Alias for Eigen::Map< MatrixX<...> >. */
template <typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Eigen::Dynamic,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using MatrixXMap = Eigen::Map< MatrixX<T, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> >;

/** Alias for Eigen::Map< MatrixX<...> const >. */
template <typename T = Real,
          int Options = DEFAULT_MATRIX_LAYOUT,
          int MaxRowsAtCompileTime = Eigen::Dynamic,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using MatrixXConstMap = Eigen::Map< MatrixX<T, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> const >;

/** Alias for Eigen::Map< VectorX<...> >. */
template <typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Eigen::Dynamic>
using VectorXMap = Eigen::Map< VectorX<T, Options, MaxRowsAtCompileTime> >;

/** Alias for Eigen::Map< VectorX<...> const >. */
template <typename T = Real,
          int Options = MatrixLayout::COLUMN_MAJOR,
          int MaxRowsAtCompileTime = Eigen::Dynamic>
using VectorXConstMap = Eigen::Map< VectorX<T, Options, MaxRowsAtCompileTime> const >;

/** Alias for Eigen::Map< RowVectorX<...> >. */
template <typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using RowVectorXMap = Eigen::Map< RowVectorX<T, Options, MaxColsAtCompileTime> >;

/** Alias for Eigen::Map< RowVectorX<...> const >. */
template <typename T = Real,
          int Options = MatrixLayout::ROW_MAJOR,
          int MaxColsAtCompileTime = Eigen::Dynamic>
using RowVectorXConstMap = Eigen::Map< RowVectorX<T, Options, MaxColsAtCompileTime> const >;

//=============================================================================================================================
// Typedef other useful Eigen classes.
//=============================================================================================================================

/** Alias for Eigen::Quaternion<...>. */
template <typename T = Real> using Quaternion = Eigen::Quaternion<T>;

} // namespace Thea

#include "MatrixUtil.hpp"

#endif
