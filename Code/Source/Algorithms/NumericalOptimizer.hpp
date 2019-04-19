//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_NumericalOptimizer_hpp__
#define __Thea_Algorithms_NumericalOptimizer_hpp__

#include "../Common.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "../Options.hpp"
#include "ScalarConstraint.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class ScalarFunction;

/**
 * Interface for numerical optimizers, that minimize an objective function f:R^n-->R under given constraints.
 *
 * To create an instance of a NumericalOptimizer, one typically loads the plugin for the relevant implementation and calls
 * NumericalOptimizerFactory::createNumericalOptimizer().
 */
class THEA_API NumericalOptimizer : public virtual AbstractNamedObject
{
  public:
    /** Destructor. */
    virtual ~NumericalOptimizer() {}

    /** Get the dimensionality of the problem. */
    virtual long dims() const = 0;

    /** Add a linear constraint to the optimization. */
    virtual void addConstraint(ScalarConstraint::ConstPtr constraint);

    /** Clear all constraints. */
    virtual void clearConstraints() = 0;

    /**
     * Minimize an objective function. A set of options may be specified to control the optimization process.
     *
     * @param objective The objective function to minimize.
     * @param hint A guessed value for the initial solution, if non-null. May be ignored by the optimizer.
     * @param options Backend-specific options. See the documentation of each derived class wrapping a particular backend.
     *
     * @return True if a local minimum was successfully found, else false. (The same value is returned by successive calls to
     *   hasSolution().)
     */
    virtual bool minimize(ScalarFunction const & objective, double const * hint = NULL,
                          AbstractOptions const * options = NULL) = 0;

    /** Was a local minimum found by the last call to minimize()? */
    virtual bool hasSolution() const = 0;

    /** Get the solution vector of the optimization problem. Valid only if hasSolution() returns true. */
    virtual double const * getSolution() const = 0;

}; // class NumericalOptimizer

/** An interface for a numerical optimizer factory. Should be implemented and registered by each actual optimizer. */
class THEA_API NumericalOptimizerFactory
{
  public:
    /** Destructor. */
    virtual ~NumericalOptimizerFactory() {}

    /**
     * Create a numerical optimizer with the given name. The numerical optimizer must be destroyed using
     * destroyNumericalOptimizer().
     */
    virtual NumericalOptimizer * createNumericalOptimizer(char const * name) = 0;

    /** Destroy a numerical optimizer created with createNumericalOptimizer(). */
    virtual void destroyNumericalOptimizer(NumericalOptimizer * optimizer) = 0;

}; // class NumericalOptimizerFactory

/** Manages available numerical optimizer factories. */
class THEA_API NumericalOptimizerManager
{
  public:
    /**
     * Install a factory for a particular numerical optimizer type. The factory pointer should not be null.
     *
     * @return True if the factory was successfully installed, false if a factory of the specified type (with case-insensitive
     *   matching) is already installed.
     */
    bool installFactory(std::string const & type, NumericalOptimizerFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for numerical optimizer of a given type. An error is thrown if no such factory has been installed. */
    NumericalOptimizerFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, NumericalOptimizerFactory *> FactoryMap;  /**< Maps numerical optimizer types to factory
                                                                                instances. */

    FactoryMap installed_factories;  ///< Set of installed factories, one for each numerical optimizer type.

}; // class NumericalOptimizerManager

} // namespace Algorithms
} // namespace Thea

#endif
