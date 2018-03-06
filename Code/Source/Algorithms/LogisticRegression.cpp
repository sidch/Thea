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
#include "FastCopy.hpp"
#include "LinearLeastSquares.hpp"
#include "../Math.hpp"
#include <limits>

namespace Thea {
namespace Algorithms {

LogisticRegression::LogisticRegression(long ndim_)
: llsq(new LinearLeastSquares(ndim_ + 1))
{
  alwaysAssertM(ndim_ >= 1, "LogisticRegression: Number of dimensions must be at least 1");
}

LogisticRegression::~LogisticRegression()
{
  delete llsq;
}

long
LogisticRegression::numDimensions() const
{
  return llsq->numDimensions() - 1;
}

long
LogisticRegression::numObservations() const
{
  return llsq->numObjectives();
}

void
LogisticRegression::addObservation(double const * x, double y)
{
  //                  y   =  1 / (1 + exp(-a - b.x))
  // =>  log(1 / y - 1)   =  -a - b.x
  // =>         a + b.x   =  log(1 / y - 1)

  alwaysAssertM(x, "LogisticRegression: Observed vector x cannot be null");

  long ndims = numDimensions();
  llsq_coeffs.resize((size_t)(ndims + 1));

  llsq_coeffs[0] = 1;
  fastCopy(x, x + ndims, &llsq_coeffs[0] + 1);

  //     c   =   log(1 / y - 1)
  // =>  y   =   1 / (1 + exp(c))
  static double const MIN_Y = 1.0 / std::min(std::numeric_limits<double>::max(), 1.0 + 1.0e+30);
  static double const MAX_Y = 1.0 / (1.0 + std::max(std::numeric_limits<double>::epsilon(), 1.0e-30));

  // std::cout << "MIN_Y = " << MIN_Y << ", MAX_Y = " << MAX_Y << std::endl;

  y = Math::clamp(y, MIN_Y, MAX_Y);
  double llsq_constant = std::log(1.0 / y - 1.0);

  // std::cout << "LogisticRegression: Adding LLSQ objective [";
  // for (size_t i = 0; i < llsq_coeffs.size(); ++i) std::cout << ' ' << llsq_coeffs[i];
  // std::cout << " ] = " << llsq_constant << std::endl;

  llsq->addObjective(&llsq_coeffs[0], llsq_constant);
}

void
LogisticRegression::clearObservations()
{
  llsq->clearObjectives();
}

bool
LogisticRegression::solve(double tolerance)
{
  return llsq->solve(LinearLeastSquares::Constraint::UNCONSTRAINED, tolerance);
}

bool
LogisticRegression::hasSolution() const
{
  return llsq->hasSolution();
}

TheaArray<double> const &
LogisticRegression::getSolution() const
{
  return llsq->getSolution();
}

} // namespace Algorithms
} // namespace Thea
