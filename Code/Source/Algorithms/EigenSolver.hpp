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
#include "../Array.hpp"
#include "../CompressedSparseMatrix.hpp"
#include "../Map.hpp"
#include "../Matrix.hpp"
#include "../MatrixFormat.hpp"
#include "../MatrixUtil.hpp"
#include "../MatrixWrapper.hpp"
#include "../NamedObject.hpp"
#include "../Options.hpp"
#include <complex>

namespace Thea {
namespace Algorithms {

/**
 * Eigenvalue solver interface. The operator matrix must currently be real, though complex eigenvalues/eigenvectors can be
 * returned. The matrix should be a (row-major or column-major) Matrix, CompressedRowMatrix, or CompressedColumnMatrix. Several
 * other matrix types may be converted to these types by the relevant constructors.
 *
 * To create an instance of an EigenSolver, one typically loads the plugin for the relevant implementation and calls
 * EigenSolverFactory::createEigenSolver().
 *
 * @note A derived class typically needs to implement only the solve() function, and optionally getPreferredFormat(). To find
 * the %eigenvalues and %eigenvectors of the operator matrix <b>A</b>, the solve() function reads the operator matrix <b>A</b>
 * from #matrix. The computed %eigenvalues, %eigenvectors and relative errors are placed in #eigenvalues, #eigenvectors and
 * #relative_errors respectively.
 *
 * FIXME: This classes currently passes STL classes across DLL boundaries. It should be an abstract base class.
 */
class THEA_API EigenSolver : private Noncopyable, public virtual NamedObject
{
  public:
    /** A complex double-precision eigenvalue. */
    typedef std::complex<double>  Eigenvalue;

    /** An eigenvector (array of double precision complex numbers). */
    typedef TheaArray< std::complex<double> >  Eigenvector;

    /** Destructor. */
    virtual ~EigenSolver() {}

    /** Set the operator matrix for the eigenproblem. */
    template <typename MatrixT> void setMatrix(MatrixT const & src)
    {
      alwaysAssertM(MatrixUtil::isSquare(src), "EigenSolver: Operator matrix is not square");
      matrix.setMatrix(src, getPreferredFormat(MatrixUtil::getFormat(src)));
    }

    /**
     * Solve the eigenproblem, optionally requesting a particular number of eigenpairs. A set of options may be specified to
     * control the solution process.
     *
     * @param num_requested_eigenpairs The number of eigenpairs to look for, which need <b>not</b> be the number
     *   returned. A negative value implies that there is no upper bound on the number of eigenpairs.
     * @param options Backend-specific options. See the documentation of each derived class wrapping a particular backend.
     *
     * @return The number of eigenpairs found.
     */
    virtual long solve(int num_requested_eigenpairs = -1, Options const & options = Options()) = 0;

    /** Get the set of eigenvalues computed by the last call to solve(). */
    TheaArray<Eigenvalue> const & getEigenvalues() const { return eigenvalues; }

    /** Get the set of eigenvectors computed by the last call to solve(). */
    TheaArray<Eigenvector> const & getEigenvectors() const { return eigenvectors; }

    /** Get the set of relative errors of the eigenpairs computed by the last call to solve(). */
    TheaArray<double> const & getRelativeErrors() const { return relative_errors; }

  protected:
    /** Get the preferred storage format for a matrix in the given input format. */
    virtual MatrixFormat getPreferredFormat(MatrixFormat input_format) { return input_format; }

    MatrixWrapper<double>   matrix;           ///< Operator matrix.
    TheaArray<Eigenvalue>   eigenvalues;      ///< Set of eigenvalues.
    TheaArray<Eigenvector>  eigenvectors;     ///< Set of eigenvectors.
    TheaArray<double>       relative_errors;  ///< Relative error for each eigenvalue.

}; // class EigenSolver

/** An interface for a eigensolver factory. Should be implemented and registered by each actual eigensolver. */
class THEA_API EigenSolverFactory
{
  public:
    /** Destructor. */
    virtual ~EigenSolverFactory() {}

    /** Create a eigensolver with the given name. The eigensolver must be destroyed using destroyEigenSolver(). */
    virtual EigenSolver * createEigenSolver(std::string const & name) = 0;

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
    typedef TheaMap<std::string, EigenSolverFactory *> FactoryMap;  ///< Maps eigensolver types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each eigensolver type.

}; // class EigenSolverManager

} // namespace Algorithms
} // namespace Thea

#endif
