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

#include "OPTPPFunction.hpp"
#include "../../Algorithms/AnalyticD1ScalarFunction.hpp"
#include "../../Algorithms/AnalyticD2ScalarFunction.hpp"
#include "../../AddressableMatrix.hpp"
#include <NLF.h>
#include <Opt.h>

namespace Thea {
namespace Algorithms {

class THEA_OPTPP_DLL_LOCAL OPTPPSymMatWrapper : public AddressableMatrix<NEWMAT::Real>
{
  public:
    OPTPPSymMatWrapper(NEWMAT::SymmetricMatrix * m_) : m(m_) {}

    long numRows() const { return m->Nrows(); }
    long numColumns() const { return m->Ncols(); }

    NEWMAT::Real const & get(long row, long col) const
    {
      // For const matrices newmat returns by value, so we const-cast it get a reference
      return const_cast<NEWMAT::SymmetricMatrix *>(m)->element((int)row, (int)col);
    }

    NEWMAT::Real & getMutable(long row, long col) { return m->element((int)row, (int)col); }

    void set(long row, long col, NEWMAT::Real const & value) { m->element((int)row, (int)col) = value; }

    void fill(NEWMAT::Real const & value) { *m = value; }

  private:
    NEWMAT::SymmetricMatrix * m;

}; // class SymMatWrapper

void
nlf0(int ndim, NEWMAT::ColumnVector const & x, OPTPP::real & fx, int & result, void * data)
{
  ScalarFunction const * func = static_cast<ScalarFunction *>(data);
  fx = static_cast<OPTPP::real>(func->valueAt(x.data()));
  result = OPTPP::NLPFunction;
}

void
nlf1(int mode, int ndim, NEWMAT::ColumnVector const & x, OPTPP::real & fx, NEWMAT::ColumnVector & gx, int & result, void * data)
{
  result = OPTPP::NLPNoOp;

  if (mode & OPTPP::NLPFunction)
    nlf0(ndim, x, fx, result, data);

  if (mode & OPTPP::NLPGradient)
  {
    ScalarFunction const * func = static_cast<ScalarFunction *>(data);
    AnalyticD1ScalarFunction const * func1 = dynamic_cast<AnalyticD1ScalarFunction const *>(func);
    alwaysAssertM(func1, "OPTPPNumericalOptimizer: Function does not have analytical first derivative");

    gx.ReSize(ndim);
    func1->gradientAt(x.data(), gx.data());
    result |= OPTPP::NLPGradient;
  }
}

void
nlf2(int mode, int ndim, NEWMAT::ColumnVector const & x, OPTPP::real & fx, NEWMAT::ColumnVector & gx,
     NEWMAT::SymmetricMatrix & Hx, int & result, void * data)
{
  result = OPTPP::NLPNoOp;

  if (mode & OPTPP::NLPFunction || mode & OPTPP::NLPGradient)
    nlf1(mode, ndim, x, fx, gx, result, data);

  if (mode & OPTPP::NLPHessian)
  {
    ScalarFunction const * func = static_cast<ScalarFunction *>(data);
    AnalyticD2ScalarFunction const * func2 = dynamic_cast<AnalyticD2ScalarFunction const *>(func);
    alwaysAssertM(func2, "OPTPPNumericalOptimizer: Function does not have analytical second derivative");

    Hx.ReSize(ndim);
    OPTPPSymMatWrapper Hwrapper(&Hx);
    func2->hessianAt(x.data(), Hwrapper);
    result |= OPTPP::NLPHessian;
  }
}

void
nlfInit(int ndim, NEWMAT::ColumnVector & x)
{
  // No-op for now, since I can't figure out a way to set a function-specific initial value here.
}

OPTPP::NLPBase *
toOPTPPFunction(ScalarFunction const & f, double const * init_pt)
{
  // Ignore init_pt for now, since I can't figure out how to bind it to the init function.

  if (dynamic_cast<AnalyticD2ScalarFunction const *>(&f))
    return new OPTPP::NLF2(f.numDimensions(), nlf2, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
  else if (dynamic_cast<AnalyticD1ScalarFunction const *>(&f))
    return new OPTPP::NLF1(f.numDimensions(), nlf1, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
  else
    return new OPTPP::NLF0(f.numDimensions(), nlf0, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
}

} // namespace Algorithms
} // namespace Thea
