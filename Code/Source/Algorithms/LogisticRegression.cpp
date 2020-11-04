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
// First version: 2012
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
    MatrixWrapper<M> a_wrap(&a);
    has_solution = llsq.solve(&a_wrap, &llsq_consts[0]);

    if (has_solution)
      solution = VectorXd::Map(const_cast<double *>(llsq.getSolution()), ndims);
    else
      solution.resize(0);
  }

  return has_solution;
}

} // namespace Algorithms
} // namespace Thea
