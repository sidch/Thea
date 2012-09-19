//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#include "LinearLeastSquares.hpp"
#include "../Matrix.hpp"
#include "NNLS/nnls.h"
#include "SVD.hpp"

namespace Thea {
namespace Algorithms {

LinearLeastSquares::LinearLeastSquares(long ndim_)
: ndim(ndim_), num_objectives(0), has_solution(true)
{
  alwaysAssertM(ndim_ > 0, "LinearLeastSquares: Number of dimensions must be positive");
  solution.resize((array_size_t)ndim_, 0);
}

void
LinearLeastSquares::addObjective(double const * coeffs, double constant)
{
  alwaysAssertM(coeffs, "LinearLeastSquares: Pointer to objective coefficients cannot be null");

  a_values.insert(a_values.end(), coeffs, coeffs + ndim);
  b.push_back(constant);
  num_objectives++;
}

void
LinearLeastSquares::clearObjectives()
{
  a_values.clear();
  num_objectives = 0;
}

bool
LinearLeastSquares::solve(Constraint constraint, double tolerance)
{
  if (num_objectives < ndim)
    THEA_WARNING << "LinearLeastSquares: Fewer objectives than dimensions -- the solution will not be unique";

  if (num_objectives <= 0)
  {
    // Arbitrary solution, set everything to zero (no real need to do this, but whatever...)
    for (array_size_t i = 0; i < solution.size(); ++i)
      solution[i] = 0;

    has_solution = true;
  }
  else
  {
    switch (constraint)
    {
      case Constraint::NON_NEGATIVE:
      {
        array_size_t m = (array_size_t)num_objectives;
        array_size_t n = (array_size_t)ndim;

        // NNLS requires a COLUMN-MAJOR (Fortran-style) matrix. Plus the values will be overwritten, so we need a copy anyway
        TheaArray<double> nnls_a(a_values.size());
        for (array_size_t i = 0; i < m; i++)
          for (array_size_t j = 0; j < n; j++)
            nnls_a[i + j * m] = a_values[i * n + j];

        TheaArray<double>  nnls_b(b);  // the values will be overwritten, so we need to make a copy
        double             rnorm;
        TheaArray<double>  w(n);
        TheaArray<double>  zz(m);
        TheaArray<int>     index(n);
        int                mode;

        int mda = (int)m, im = (int)m, in = (int)n;
        nnls_c(&nnls_a[0], &mda, &im, &in, &nnls_b[0], &solution[0], &rnorm, &w[0], &zz[0], &index[0], &mode);

        if (mode == 1)
          has_solution = true;
        else
        {
          switch (mode)
          {
            // Should never be 2 since we've checked for this above, but recheck all the same
            case 2:  THEA_DEBUG << "LinearLeastSquares: NNLS error (bad problem dimensions)"; break;
            case 3:  THEA_DEBUG << "LinearLeastSquares: NNLS error (iteration count exceeded)"; break;
            default: THEA_DEBUG << "LinearLeastSquares: Unknown NNLS error";
          }

          has_solution = false;
        }

        break;
      }

      case Constraint::UNCONSTRAINED:
      {
        // Create a matrix wrapping the objective coefficients
        Matrix<double, MatrixLayout::ROW_MAJOR> a(&a_values[0], num_objectives, ndim);

        Matrix<double> pseudo_inverse;
        has_solution = SVD::pseudoInverse(a, pseudo_inverse, tolerance);
        if (has_solution)
          pseudo_inverse.postmulVector(&b[0], &solution[0]);

        break;
      }

      default:
        throw Error("LinearLeastSquares: Unsupported constraint");
    }
  }

  return has_solution;
}

} // namespace Algorithms
} // namespace Thea
