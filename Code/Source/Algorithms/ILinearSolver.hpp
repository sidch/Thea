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

#ifndef __Thea_Algorithms_ILinearSolver_hpp__
#define __Thea_Algorithms_ILinearSolver_hpp__

#include "../Common.hpp"
#include "../IMatrix.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "../Options.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for solvers of systems of linear equations. The matrices must currently be real.
 *
 * If the matrix is square and of full rank, the solver should compute the unique solution. Else, if the matrix is rectangular,
 * it should try to return an approximate solution (typically one that minimizes the squared error). In the latter case, classes
 * implementing the ILinearSolver interface may be considered alternatives to the default linear least-squares solvers provided
 * by StdLinearSolver.
 *
 * To create an instance of a ILinearSolver, one typically loads the plugin for the relevant implementation and calls
 * ILinearSolverFactory::createLinearSolver().
 */
class THEA_API ILinearSolver : public virtual INamedObject
{
  public:
    /** Destructor. */
    virtual ~ILinearSolver() = 0;

    /**
     * Solve the system of linear equations. A set of options may be specified to control the solution process.
     *
     * @param a The coefficient matrix A in the system Ax = b.
     * @param b The constant vector b in the system Ax = b.
     * @param options Backend-specific options. See the documentation of each derived class wrapping a particular backend.
     *
     * @return True if the system was successfully solved, else false. (The same value is returned by subsequent calls to
     *   hasSolution().)
     */
    virtual int8 THEA_ICALL solve(IMatrix<float64> const * a, float64 const * b, IOptions const * options = nullptr) = 0;

    /** Get the size of the solution, which is also the number of columns of the coefficient matrix A. */
    virtual int64 THEA_ICALL dims() const = 0;

    /** Was the linear system successfully solved by the last call to solve()? */
    virtual int8 THEA_ICALL hasSolution() const = 0;

    /** Get the solution vector x of the linear system Ax = b. Valid only if hasSolution() returns true. */
    virtual float64 const * THEA_ICALL getSolution() const = 0;

    /**
     * If the squared error || Ax - b ||^2 was computed during the solution process, put it in \a err and return true. Else
     * return false.
     */
    virtual int8 THEA_ICALL getSquaredError(float64 * err) const = 0;

}; // class ILinearSolver

inline ILinearSolver::~ILinearSolver() {}

/** Interface for a linear solver factory. Should be implemented and registered by each actual linear solver. */
class THEA_API ILinearSolverFactory
{
  public:
    /** Destructor. */
    virtual ~ILinearSolverFactory() = 0;

    /** Create a linear solver with the given name. The linear solver must be destroyed using destroyLinearSolver(). */
    virtual ILinearSolver * THEA_ICALL createLinearSolver(char const * name) = 0;

    /** Destroy a linear solver created with createLinearSolver(). */
    virtual int8 THEA_ICALL destroyLinearSolver(ILinearSolver * linear_solver) = 0;

}; // class ILinearSolverFactory

inline ILinearSolverFactory::~ILinearSolverFactory() {}

/** Manages available linear solver factories. */
class THEA_API LinearSolverManager
{
  public:
    /**
     * Install a factory for a particular linear solver type. The factory pointer should not be null.
     *
     * @return True if the factory was successfully installed, false if a factory of the specified type (with case-insensitive
     *   matching) is already installed.
     */
    bool installFactory(std::string const & type, ILinearSolverFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for linear solver of a given type. An error is thrown if no such factory has been installed. */
    ILinearSolverFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, ILinearSolverFactory *> FactoryMap;  ///< Maps linear solver types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each linear solver type.

}; // class LinearSolverManager

} // namespace Algorithms
} // namespace Thea

#endif
