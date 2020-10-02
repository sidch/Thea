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

#ifndef __Thea_Algorithms_BiasedNormalizedCuts_hpp__
#define __Thea_Algorithms_BiasedNormalizedCuts_hpp__

#include "../Common.hpp"
#include "../CompressedSparseMatrix.hpp"
#include "IEigenSolver.hpp"
#include <complex>

namespace Thea {
namespace Algorithms {

/**
 */
class THEA_API BiasedNormalizedCuts
{
  public:
    template <typename MatrixT>
    static void compute(MatrixT const & adjacency_matrix, Array<intx> const & vertices_of_interest, double correlation,
                        Array<double> & vertex_indicators)
    {
    }

}; // class BiasedNormalizedCuts

} // namespace Algorithms
} // namespace Thea

#endif
