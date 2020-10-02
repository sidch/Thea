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
#include "IScalarConstraint.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class IScalarFunction;

/**
 * Interface for numerical optimizers, that minimize an objective function f:R^n-->R under given constraints.
 *
 * To create an instance of a INumericalOptimizer, one typically loads the plugin for the relevant implementation and calls
 * INumericalOptimizerFactory::createNumericalOptimizer().
 */
class THEA_API INumericalOptimizer : public virtual INamedObject
{
  public:
    /** Destructor. */
    virtual ~INumericalOptimizer() = 0;

    /** Get the dimensionality of the problem. */
    virtual int64 THEA_ICALL dims() const = 0;

    /** Add a linear constraint to the optimization. The constraint object must exist as long as this object does. */
    virtual void THEA_ICALL addConstraint(IScalarConstraint const * constraint) = 0;

    /** Clear all constraints. */
    virtual void THEA_ICALL clearConstraints() = 0;

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
    virtual int8 THEA_ICALL minimize(IScalarFunction const * objective, float64 const * hint = nullptr,
                          IOptions const * options = nullptr) = 0;

    /** Was a local minimum found by the last call to minimize()? */
    virtual int8 THEA_ICALL hasSolution() const = 0;

    /** Get the solution vector of the optimization problem. Valid only if hasSolution() returns true. */
    virtual float64 const * THEA_ICALL getSolution() const = 0;

}; // class INumericalOptimizer

inline INumericalOptimizer::~INumericalOptimizer() {}

/** Interface for a numerical optimizer factory. Should be implemented and registered by each actual optimizer. */
class THEA_API INumericalOptimizerFactory
{
  public:
    /** Destructor. */
    virtual ~INumericalOptimizerFactory() = 0;

    /**
     * Create a numerical optimizer with the given name. The numerical optimizer must be destroyed using
     * destroyNumericalOptimizer().
     */
    virtual INumericalOptimizer * THEA_ICALL createNumericalOptimizer(char const * name) = 0;

    /** Destroy a numerical optimizer created with createNumericalOptimizer(). */
    virtual void THEA_ICALL destroyNumericalOptimizer(INumericalOptimizer * optimizer) = 0;

}; // class INumericalOptimizerFactory

inline INumericalOptimizerFactory::~INumericalOptimizerFactory() {}

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
    bool installFactory(std::string const & type, INumericalOptimizerFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for numerical optimizer of a given type. An error is thrown if no such factory has been installed. */
    INumericalOptimizerFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, INumericalOptimizerFactory *> FactoryMap;  /**< Maps numerical optimizer types to factory
                                                                                instances. */

    FactoryMap installed_factories;  ///< Set of installed factories, one for each numerical optimizer type.

}; // class NumericalOptimizerManager

} // namespace Algorithms
} // namespace Thea

#endif
