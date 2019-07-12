//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#include "LogisticRegression.hpp"
#include "StdLinearSolver.hpp"
#include "../Math.hpp"
#include "../MatrixWrapper.hpp"
#include <limits>

namespace Thea {
namespace Algorithms {

LogisticRegression::LogisticRegression(intx ndims_)
: ndims(ndims_)
{
  alwaysAssertM(ndims_ >= 1, "LogisticRegression: Number of dimensions must be at least 1");
}

LogisticRegression::~LogisticRegression()
{
}

void
LogisticRegression::addObservation(double const * x, double y)
{
  //                  y   =  1 / (1 + exp(-a - b.x))
  // =>  log(1 / y - 1)   =  -a - b.x
  // =>         a + b.x   =  log(1 / y - 1)

  alwaysAssertM(x, "LogisticRegression: Observed vector x cannot be null");

  llsq_coeffs.push_back(1);
  llsq_coeffs.insert(llsq_coeffs.end(), x, x + ndims);

  //     c   =   log(1 / y - 1)
  // =>  y   =   1 / (1 + exp(c))
  static double const MIN_Y = 1.0 / std::min(std::numeric_limits<double>::max(), 1.0 + 1.0e+30);
  static double const MAX_Y = 1.0 / (1.0 + std::max(std::numeric_limits<double>::epsilon(), 1.0e-30));

  // std::cout << "MIN_Y = " << MIN_Y << ", MAX_Y = " << MAX_Y << std::endl;

  y = Math::clamp(y, MIN_Y, MAX_Y);
  llsq_consts.push_back(std::log(1.0 / y - 1.0));

  // std::cout << "LogisticRegression: Adding LLSQ objective [";
  // for (size_t i = llsq_coeffs.size() - (size_t)ndims - 1; i < llsq_coeffs.size(); ++i) std::cout << ' ' << llsq_coeffs[i];
  // std::cout << " ] = " << llsq_consts.back() << std::endl;
}

void
LogisticRegression::clearObservations()
{
  llsq_coeffs.clear();
  llsq_consts.clear();
}

bool
LogisticRegression::solve(double tolerance)
{
  if (numObservations() <= 0)
  {
    THEA_WARNING << "LogisticRegression: Solving empty problem";
    has_solution = true;
    solution.resize(0);
  }
  else
  {
    StdLinearSolver llsq(StdLinearSolver::Method::DEFAULT, StdLinearSolver::Constraint::UNCONSTRAINED);
    llsq.setTolerance(tolerance);

    typedef Eigen::Map< MatrixX<double, MatrixLayout::ROW_MAJOR> > M;
    M a(&llsq_coeffs[0], numObservations(), ndims);
    has_solution = llsq.solve(MatrixWrapper<M>(&a), &llsq_consts[0]);

    if (has_solution)
      solution = VectorXdMap(const_cast<double *>(llsq.getSolution()), ndims);
    else
      solution.resize(0);
  }

  return has_solution;
}

} // namespace Algorithms
} // namespace Thea
