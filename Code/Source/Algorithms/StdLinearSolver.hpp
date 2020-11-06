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

#ifndef __Thea_Algorithms_StdLinearSolver_hpp__
#define __Thea_Algorithms_StdLinearSolver_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../MatVec.hpp"
#include "../SparseMatVec.hpp"
#include "ILinearSolver.hpp"
#include "NNLS/nnls.h"
#include <Eigen/QR>
#include <Eigen/SVD>

namespace Thea {
namespace Algorithms {

// Forward declarations
namespace StdLinearSolverInternal { class StdLinearSolverImpl; }

/**
 * Solve dense and sparse linear systems of the form Ax = b for x. This class implements the ILinearSolver interface to provide a
 * variety of built-in algorithms. Other solvers may be provided by plugins implementing the ILinearSolver interface.
 */
class THEA_API StdLinearSolver : public virtual ILinearSolver, public NamedObject
{
  public:
    THEA_DECL_SMART_POINTERS(StdLinearSolver)

    /** The constraints to be imposed on the solution (enum class). */
    struct THEA_API Constraint
    {
      /** Supported values. */
      enum Value
      {
        UNCONSTRAINED,  ///< The solution is not constrained.
        NON_NEGATIVE    ///< The solution must have non-negative elements only.
      };

      THEA_ENUM_CLASS_BODY(Constraint)

      THEA_ENUM_CLASS_STRINGS_BEGIN(Constraint)
        THEA_ENUM_CLASS_STRING(UNCONSTRAINED,  "unconstrained")
        THEA_ENUM_CLASS_STRING(NON_NEGATIVE,   "non-negative")
      THEA_ENUM_CLASS_STRINGS_END(Constraint)
    };

    /** Solution methods (enum class). */
    struct THEA_API Method
    {
      /** Supported values. */
      enum Value
      {
        DEFAULT,                             ///< Automatically pick a solution method.

        // Dense solvers
        HOUSEHOLDER_QR,                      ///< Eigen's HouseholderQR solver.
        COL_PIV_HOUSEHOLDER_QR,              ///< Eigen's ColPivHouseholderQR solver.
        FULL_PIV_HOUSEHOLDER_QR,             ///< Eigen's FullPivHouseholderQR solver.
        COMPLETE_ORTHOGONAL_DECOMPOSITION,   ///< Eigen's CompleteOrthogonalDecomposition solver.
        BDCSVD,                              ///< Eigen's BDCSVD solver.
        NNLS,                                ///< NNLS non-negative least squares solver.

        // Sparse solvers
        SIMPLICIALT_LLT,                     ///< Eigen's SimplicialLLT solver.
        SIMPLICIALT_LDLT,                    ///< Eigen's SimplicialLDLT solver.
        SPARSE_LU,                           ///< Eigen's SparseLU solver.
        SPARSE_QR,                           ///< Eigen's SparseQR solver.
        CONJUGATE_GRADIENT,                  ///< Eigen's ConjugateGradient solver.
        LEAST_SQUARES_CONJUGATE_GRADIENT,    ///< Eigen's LeastSquaresConjugateGradient solver.
        BICGSTAB,                            ///< Eigen's BiCGSTAB solver.
      };

      THEA_ENUM_CLASS_BODY(Method)

      THEA_ENUM_CLASS_STRINGS_BEGIN(Method)
        THEA_ENUM_CLASS_STRING(DEFAULT,                            "default")

        THEA_ENUM_CLASS_STRING(HOUSEHOLDER_QR,                     "HouseholderQR")
        THEA_ENUM_CLASS_STRING(COL_PIV_HOUSEHOLDER_QR,             "ColPivHouseholderQR")
        THEA_ENUM_CLASS_STRING(FULL_PIV_HOUSEHOLDER_QR,            "FullPivHouseholderQR")
        THEA_ENUM_CLASS_STRING(COMPLETE_ORTHOGONAL_DECOMPOSITION,  "CompleteOrthogonalDecomposition")
        THEA_ENUM_CLASS_STRING(BDCSVD,                             "BDSVD")
        THEA_ENUM_CLASS_STRING(NNLS,                               "NNLS")

        THEA_ENUM_CLASS_STRING(SIMPLICIALT_LLT,                    "SimplicialLLT")
        THEA_ENUM_CLASS_STRING(SIMPLICIALT_LDLT,                   "SimplicialLDLT")
        THEA_ENUM_CLASS_STRING(SPARSE_LU,                          "SparseLU")
        THEA_ENUM_CLASS_STRING(SPARSE_QR,                          "SparseQR")
        THEA_ENUM_CLASS_STRING(CONJUGATE_GRADIENT,                 "ConjugateGradient")
        THEA_ENUM_CLASS_STRING(LEAST_SQUARES_CONJUGATE_GRADIENT,   "LeastSquaresConjugateGradient")
        THEA_ENUM_CLASS_STRING(BICGSTAB,                           "BiCGSTAB")
      THEA_ENUM_CLASS_STRINGS_END(Method)
    };

    /** Constructor. */
    StdLinearSolver(Method method_ = Method::DEFAULT, Constraint constraint_ = Constraint::UNCONSTRAINED);

    /** Destructor. */
    ~StdLinearSolver();

    /** Get the solution method. */
    Method getMethod() const;

    /** Get the solution constraint. */
    Constraint getConstraint() const;

    /* Get the solution tolerance/threshold. A negative number implies the default tolerance. */
    double getTolerance() const;

    /* Get the maximum number of solver iterations, if the solver is iterative. */
    intx maxIterations() const;

    /** Set the solution method. */
    void setMethod(Method method_);

    /** Set the solution constraint. */
    void setConstraint(Constraint constraint_);

    /* Set the solution tolerance/threshold. */
    void setTolerance(float64 tol);

    /* Get the maximum number of solver iterations, if the solver is iterative. A negative number implies the default value. */
    void setMaxIterations(intx max_iters_);

    /** Solve the linear system Ax = b for a dense double-precision matrix A. */
    bool solve(Eigen::Ref< MatrixXd > const & a, double const * b, IOptions const * options = nullptr);

    /** Solve the linear system Ax = b for a sparse double-precision matrix A. */
    bool solve(Eigen::Ref< SparseMatrix<double> > const & a, double const * b, IOptions const * options = nullptr);

    // Functions from ILinearSolver
    int8 THEA_ICALL solve(IMatrix<float64> const * a, float64 const * b, IOptions const * options = nullptr);
    int64 THEA_ICALL dims() const;
    int8 THEA_ICALL hasSolution() const;
    float64 const * THEA_ICALL getSolution() const;
    int8 THEA_ICALL getSquaredError(float64 * err) const;

  private:
    StdLinearSolverInternal::StdLinearSolverImpl * impl;  ///< Contains base function implementations using PIMPL idiom.

}; // class StdLinearSolver

} // namespace Algorithms
} // namespace Thea

#endif
