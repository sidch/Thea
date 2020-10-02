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
// First version: 2009
//
//============================================================================

#include "OPTPPFunction.hpp"
#include "../../Algorithms/IAnalyticD1ScalarFunction.hpp"
#include "../../Algorithms/IAnalyticD2ScalarFunction.hpp"
#include "../../IAddressableMatrix.hpp"
#include <NLF.h>
#include <Opt.h>

namespace Thea {
namespace Algorithms {

class THEA_OPTPP_DLL_LOCAL OPTPPSymMatWrapper : public virtual IAddressableMatrix<NEWMAT::Real>
{
  public:
    OPTPPSymMatWrapper(NEWMAT::SymmetricMatrix * m_) : m(m_) {}

    intx rows() const { return m->Nrows(); }
    intx cols() const { return m->Ncols(); }

    NEWMAT::Real const & get(intx row, intx col) const
    {
      // For const matrices newmat returns by value, so we const-cast it get a reference
      return const_cast<NEWMAT::SymmetricMatrix *>(m)->element((int)row, (int)col);
    }

    NEWMAT::Real & getMutable(intx row, intx col) { return m->element((int)row, (int)col); }

    void set(intx row, intx col, NEWMAT::Real const & value) { m->element((int)row, (int)col) = value; }

  private:
    NEWMAT::SymmetricMatrix * m;

}; // class SymMatWrapper

void
nlf0(int ndim, NEWMAT::ColumnVector const & x, OPTPP::real & fx, int & result, void * data)
{
  IScalarFunction const * func = static_cast<IScalarFunction *>(data);
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
    IScalarFunction const * func = static_cast<IScalarFunction *>(data);
    IAnalyticD1ScalarFunction const * func1 = dynamic_cast<IAnalyticD1ScalarFunction const *>(func);
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
    IScalarFunction const * func = static_cast<IScalarFunction *>(data);
    IAnalyticD2ScalarFunction const * func2 = dynamic_cast<IAnalyticD2ScalarFunction const *>(func);
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
toOPTPPFunction(IScalarFunction const & f, double const * init_pt)
{
  // Ignore init_pt for now, since I can't figure out how to bind it to the init function.

  if (dynamic_cast<IAnalyticD2ScalarFunction const *>(&f))
    return new OPTPP::NLF2(f.dims(), nlf2, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
  else if (dynamic_cast<IAnalyticD1ScalarFunction const *>(&f))
    return new OPTPP::NLF1(f.dims(), nlf1, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
  else
    return new OPTPP::NLF0(f.dims(), nlf0, nlfInit, const_cast<void *>(static_cast<void const *>(&f)));
}

} // namespace Algorithms
} // namespace Thea
