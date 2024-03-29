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
// First version: 2011
//
//============================================================================

#ifndef __Thea_MatrixFormat_hpp__
#define __Thea_MatrixFormat_hpp__

#include "Common.hpp"

// Resolve a conflict between X.h (which #defines Success) and Eigen (which has an enum called Success). Fortunately both
// headers currently assign the same value to Success (0) so the #undef should be all that's needed.
// https://stackoverflow.com/questions/22400905/eigen-and-cimg-compatibility-issues
#ifdef Success
#  undef Success
#endif
#include <Eigen/Core>

namespace Thea {

/** %Matrix layouts (enum class). */
struct THEA_API MatrixLayout
{
  /** Supported values. */
  enum Value
  {
    ROW_MAJOR     =  (int)Eigen::RowMajor,  ///< Row-major layout.
    COLUMN_MAJOR  =  (int)Eigen::ColMajor,  ///< Column-major layout.
  };

  THEA_ENUM_CLASS_BODY(MatrixLayout)
};

} // namespace Thea

#endif
