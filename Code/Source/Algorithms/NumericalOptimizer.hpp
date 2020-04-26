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
// First version: 2011
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
    virtual ~NumericalOptimizer() = 0;

    /** Get the dimensionality of the problem. */
    virtual int64 dims() const = 0;

    /** Add a linear constraint to the optimization. The constraint object must exist as long as this object does. */
    virtual void addConstraint(ScalarConstraint const * constraint) = 0;

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
    virtual int8 minimize(ScalarFunction const & objective, float64 const * hint = nullptr,
                          AbstractOptions const * options = nullptr) = 0;

    /** Was a local minimum found by the last call to minimize()? */
    virtual int8 hasSolution() const = 0;

    /** Get the solution vector of the optimization problem. Valid only if hasSolution() returns true. */
    virtual float64 const * getSolution() const = 0;

}; // class NumericalOptimizer

// Pure virtual destructor should have implementation
inline NumericalOptimizer::~NumericalOptimizer() {}

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
