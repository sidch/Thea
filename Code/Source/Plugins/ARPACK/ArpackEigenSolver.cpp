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

#include "ArpackEigenSolver.hpp"
#include <algorithm>
#include <cstring>

namespace Thea {
namespace Algorithms {

ArpackEigenSolver::ArpackEigenSolver(std::string const & name_)
: NamedObject(name_)
{}

ArpackEigenSolver::~ArpackEigenSolver()
{}

int64
ArpackEigenSolver::solve(IMatrix<float64> const * m, int8 compute_eigenvectors, int64 num_requested_eigenpairs,
                         IOptions const * options)
{
  if (!m || !Math::isSquare(*m))
  {
    THEA_ERROR << getName() << ": Operator matrix must be square and non-null";
    return -1;
  }

  eigenvalues [0].clear(); eigenvalues [1].clear();
  eigenvectors[0].clear(); eigenvectors[1].clear();

  ndims = m->rows();
  if (ndims <= 0)
  {
    THEA_WARNING << getName() << ": Operator matrix has zero dimensions";
    return 0;
  }

  if (num_requested_eigenpairs < 0)
    num_requested_eigenpairs = ndims;

  char     const * DEFAULT_WHICH = "LM";
  int32    const   DEFAULT_NCV = std::min(2 * num_requested_eigenpairs + 1, ndims - 1);
  float64  const   DEFAULT_TOL = std::max(1e-7, 5e-9 * ndims);
  int32    const   DEFAULT_MAXIT = 100 * num_requested_eigenpairs;
  int32    const   DEFAULT_SIGMA = 0;

  std::string which_s = options ? options->getString("which", DEFAULT_WHICH) : DEFAULT_WHICH;
  char which[3];
  std::strncpy(which, which_s.c_str(), 2);  // make a copy since ARPACK++ doesn't const-protect this string
  which[2] = 0;

  int32 ncv    =  options ? (int32)options->getInteger("ncv",   DEFAULT_NCV  ) : DEFAULT_NCV;
  float64 tol  =  options ?        options->getFloat  ("tol",   DEFAULT_TOL  ) : DEFAULT_TOL;
  int32 maxit  =  options ? (int32)options->getInteger("maxit", DEFAULT_MAXIT) : DEFAULT_MAXIT;

  float64 sigma = 0;
  int8 shift_invert = options && options->hasOption("sigma");
  if (shift_invert)
    sigma = options->getFloat("sigma", DEFAULT_SIGMA);

  THEA_DEBUG << getName() << ": nev = " << num_requested_eigenpairs << ", which = " << which << ", ncv = " << ncv
             << ", tol = " << tol << ", maxit = " << maxit;

  if (m->asSparse() && m->asSparse()->asCompressed())
    return solveSparse(*m->asSparse()->asCompressed(), num_requested_eigenpairs, shift_invert, sigma, which, ncv, tol, maxit);
  else if (m->asAddressable() && m->asAddressable()->asDense())
    return solveDense(*m->asAddressable()->asDense(), num_requested_eigenpairs, shift_invert, sigma, which, ncv, tol, maxit);
  else
  {
    THEA_ERROR << getName() << ": Unsupported matrix format";
    return -1;
  }
}

int8
ArpackEigenSolver::getEigenvalue(int64 i, float64 * re, float64 * im) const
{
  if (i < 0 || i >= (int64)eigenvalues[0].size())
  {
    THEA_ERROR << getName() << ": Index of eigenvalue out of range";
    return false;
  }

  if (i < 0 || i >= (int64)eigenvalues[0].size())
  {
    THEA_ERROR << getName() << ": Index of eigenvalue out of range";
    return false;
  }

  if (re) *re = eigenvalues[0][(size_t)i];
  if (im) *im = eigenvalues[1][(size_t)i];

  return true;
}

int8
ArpackEigenSolver::getEigenvector(int64 i, float64 const ** re, float64 const ** im) const
{
  if (i < 0 || i >= (int64)eigenvalues[0].size())
  {
    THEA_ERROR << getName() << ": Index of eigenvector out of range";
    return false;
  }

  if (i >= (int64)eigenvectors[0].size())  // eigenvectors were not computed
    return false;

  if (re) *re = eigenvectors[0][(size_t)i].data();
  if (im) *im = eigenvectors[1][(size_t)i].data();

  return true;
}

int8
ArpackEigenSolver::getRelativeError(int64 i, float64 * error) const
{
  THEA_ERROR << getName() << ": Relative errors not available";
  return false;
}

ArpackEigenSolverFactory::~ArpackEigenSolverFactory()
{
  destroyAllEigenSolvers();
}

IEigenSolver *
ArpackEigenSolverFactory::createEigenSolver(char const * name)
{
  ArpackEigenSolver * es = new ArpackEigenSolver(name);
  eigen_solvers.insert(es);
  return es;
}

int8
ArpackEigenSolverFactory::destroyEigenSolver(IEigenSolver * eigen_solver)
{
  eigen_solvers.erase(eigen_solver);
  delete eigen_solver;

  return true;
}

int8
ArpackEigenSolverFactory::destroyAllEigenSolvers()
{
  for (EigenSolverSet::iterator ei = eigen_solvers.begin(); ei != eigen_solvers.end(); ++ei)
    delete *ei;

  eigen_solvers.clear();

  return true;
}

} // namespace Algorithms
} // namespace Thea
