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

#ifndef __Thea_Algorithms_EigenSolver_hpp__
#define __Thea_Algorithms_EigenSolver_hpp__

#include "../Common.hpp"
#include "../AbstractMatrix.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "../Options.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Eigenvalue solver interface. The operator matrix must currently be real (and float64-precision), though complex
 * eigenvalues/eigenvectors can be returned.
 *
 * To create an instance of an EigenSolver, one typically loads the plugin for the relevant implementation and calls
 * EigenSolverFactory::createEigenSolver().
 */
class THEA_API EigenSolver : public virtual AbstractNamedObject
{
  public:
    /** Destructor. */
    virtual ~EigenSolver() = 0;

    /**
     * Find the eigenvalues (and optionally eigenvectors) of a matrix, optionally requesting a particular number of eigenpairs.
     * A set of additional options may be specified to control the solution process.
     *
     * @param m The operator matrix for the eigenproblem. Must be square.
     * @param compute_eigenvectors If false, only the eigenvalues will be computed, not the eigenvectors. getEigenvalue() will
     *   fail on subsequent calls.
     * @param num_requested_eigenpairs The number of eigenpairs to look for, which need <b>not</b> be the number
     *   returned. A negative value implies that there is no upper bound on the number of eigenpairs.
     * @param options Backend-specific options. See the documentation of each derived class wrapping a particular backend.
     *
     * @return The number of eigenpairs found, or a negative number on error.
     */
    virtual int64 solve(AbstractMatrix<float64> const & m, int8 compute_eigenvectors = true,
                        int64 num_requested_eigenpairs = -1, AbstractOptions const * options = NULL) = 0;

    /** Get the size of each eigenvector, which is also the number of rows (or columns) of the input operator matrix. */
    virtual int64 dims() const = 0;

    /** Get the number of eigenpairs computed by the last call to solve(). */
    virtual int64 numEigenpairs() const = 0;

    /**
     * Get the \a i'th eigenvalue computed by the last call to solve(), as a complex number with real and imaginary parts.
     *
     * @return True on success, false on error.
     */
    virtual int8 getEigenvalue(int64 i, float64 & re, float64 & im) const = 0;

    /**
     * Get the \a i'th eigenvector computed by the last call to solve(), if available. The vector is returned as two arrays of
     * floating-point numbers, one for the real parts of the coordinates and one for the imaginary parts of the coordinates.
     * Each array has dims() entries.
     *
     * @return True on success, false on error.
     */
    virtual int8 getEigenvector(int64 i, float64 const * & re, float64 const * & im) const = 0;

    /** Check if the solver has stored the relative error of each eigenpair. */
    virtual int8 hasRelativeErrors() const = 0;

    /**
     * Get the relative error of the \a i'th eigenpair computed by the last call to solve(), if available.
     *
     * @return True on success, false on error.
     */
    virtual int8 getRelativeError(int64 i, float64 & error) const = 0;

}; // class EigenSolver

// Pure virtual destructor should have implementation
inline EigenSolver::~EigenSolver() {}

/** An interface for a eigensolver factory. Should be implemented and registered by each actual eigensolver. */
class THEA_API EigenSolverFactory
{
  public:
    /** Destructor. */
    virtual ~EigenSolverFactory() {}

    /** Create a eigensolver with the given name. The eigensolver must be destroyed using destroyEigenSolver(). */
    virtual EigenSolver * createEigenSolver(char const * name) = 0;

    /** Destroy a eigensolver created with createEigenSolver(). */
    virtual void destroyEigenSolver(EigenSolver * eigen_solver) = 0;

}; // class EigenSolverFactory

/** Manages available eigensolver factories. */
class THEA_API EigenSolverManager
{
  public:
    /**
     * Install a factory for a particular eigensolver type. The factory pointer should not be null.
     *
     * @return True if the factory was successfully installed, false if a factory of the specified type (with case-insensitive
     *   matching) is already installed.
     */
    bool installFactory(std::string const & type, EigenSolverFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for eigensolver of a given type. An error is thrown if no such factory has been installed. */
    EigenSolverFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, EigenSolverFactory *> FactoryMap;  ///< Maps eigensolver types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each eigensolver type.

}; // class EigenSolverManager

} // namespace Algorithms
} // namespace Thea

#endif
