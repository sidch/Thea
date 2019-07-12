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

#ifndef __Thea_Algorithms_LinearSolver_hpp__
#define __Thea_Algorithms_LinearSolver_hpp__

#include "../Common.hpp"
#include "../AbstractMatrix.hpp"
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
 * implementing the LinearSolver interface may be considered alternatives to the default linear least-squares solvers provided
 * by StdLinearSolver.
 *
 * To create an instance of a LinearSolver, one typically loads the plugin for the relevant implementation and calls
 * LinearSolverFactory::createLinearSolver().
 */
class THEA_API LinearSolver : public virtual AbstractNamedObject
{
  public:
    /** Destructor. */
    virtual ~LinearSolver() = 0;

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
    virtual int8 solve(AbstractMatrix<float64> const & a, float64 const * b, AbstractOptions const * options = NULL) = 0;

    /** Get the size of the solution, which is also the number of columns of the coefficient matrix A. */
    virtual int64 dims() const = 0;

    /** Was the linear system successfully solved by the last call to solve()? */
    virtual int8 hasSolution() const = 0;

    /** Get the solution vector x of the linear system Ax = b. Valid only if hasSolution() returns true. */
    virtual float64 const * getSolution() const = 0;

    /**
     * If the squared error || Ax - b ||^2 was computed during the solution process, put it in \a err and return true. Else
     * return false.
     */
    virtual int8 getSquaredError(float64 & err) const = 0;

}; // class LinearSolver

// Pure virtual destructor should have implementation
inline LinearSolver::~LinearSolver() {}

/** An interface for a linear solver factory. Should be implemented and registered by each actual linear solver. */
class THEA_API LinearSolverFactory
{
  public:
    /** Destructor. */
    virtual ~LinearSolverFactory() {}

    /** Create a linear solver with the given name. The linear solver must be destroyed using destroyLinearSolver(). */
    virtual LinearSolver * createLinearSolver(char const * name) = 0;

    /** Destroy a linear solver created with createLinearSolver(). */
    virtual void destroyLinearSolver(LinearSolver * linear_solver) = 0;

}; // class LinearSolverFactory

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
    bool installFactory(std::string const & type, LinearSolverFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for linear solver of a given type. An error is thrown if no such factory has been installed. */
    LinearSolverFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, LinearSolverFactory *> FactoryMap;  ///< Maps linear solver types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each linear solver type.

}; // class LinearSolverManager

} // namespace Algorithms
} // namespace Thea

#endif
