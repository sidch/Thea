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

#ifndef __Thea_Algorithms_IEigenSolver_hpp__
#define __Thea_Algorithms_IEigenSolver_hpp__

#include "../Common.hpp"
#include "../IMatrix.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "../Options.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Eigenvalue solver interface. The operator matrix must currently be real (and float64-precision), though complex
 * eigenvalues/eigenvectors can be returned.
 *
 * To create an instance of an IEigenSolver, one typically loads the plugin for the relevant implementation and calls
 * IEigenSolverFactory::createEigenSolver().
 */
class THEA_API IEigenSolver : public virtual INamedObject
{
  public:
    /** Destructor. */
    virtual ~IEigenSolver() = 0;

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
    virtual int64 THEA_ICALL solve(IMatrix<float64> const * m, int8 compute_eigenvectors = true,
                                   int64 num_requested_eigenpairs = -1, IOptions const * options = nullptr) = 0;

    /** Get the size of each eigenvector, which is also the number of rows (or columns) of the input operator matrix. */
    virtual int64 THEA_ICALL dims() const = 0;

    /** Get the number of eigenpairs computed by the last call to solve(). */
    virtual int64 THEA_ICALL numEigenpairs() const = 0;

    /**
     * Get the \a i'th eigenvalue computed by the last call to solve(), as a complex number with real and imaginary parts.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getEigenvalue(int64 i, float64 * re, float64 * im) const = 0;

    /**
     * Get the \a i'th eigenvector computed by the last call to solve(), if available. The vector is returned as two arrays of
     * floating-point numbers, one for the real parts of the coordinates and one for the imaginary parts of the coordinates.
     * Each array has dims() entries.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getEigenvector(int64 i, float64 const ** re, float64 const ** im) const = 0;

    /** Check if the solver has stored the relative error of each eigenpair. */
    virtual int8 THEA_ICALL hasRelativeErrors() const = 0;

    /**
     * Get the relative error of the \a i'th eigenpair computed by the last call to solve(), if available.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getRelativeError(int64 i, float64 * error) const = 0;

}; // class IEigenSolver

inline IEigenSolver::~IEigenSolver() {}

/** Interface for a eigensolver factory. Should be implemented and registered by each actual eigensolver. */
class THEA_API IEigenSolverFactory
{
  public:
    /** Destructor. */
    virtual ~IEigenSolverFactory() = 0;

    /** Create a eigensolver with the given name. The eigensolver must be destroyed using destroyEigenSolver(). */
    virtual IEigenSolver * THEA_ICALL createEigenSolver(char const * name) = 0;

    /** Destroy a eigensolver created with createEigenSolver(). */
    virtual int8 THEA_ICALL destroyEigenSolver(IEigenSolver * eigen_solver) = 0;

}; // class IEigenSolverFactory

inline IEigenSolverFactory::~IEigenSolverFactory() {}

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
    bool installFactory(std::string const & type, IEigenSolverFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for eigensolver of a given type. An error is thrown if no such factory has been installed. */
    IEigenSolverFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, IEigenSolverFactory *> FactoryMap;  ///< Maps eigensolver types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each eigensolver type.

}; // class EigenSolverManager

} // namespace Algorithms
} // namespace Thea

#endif
