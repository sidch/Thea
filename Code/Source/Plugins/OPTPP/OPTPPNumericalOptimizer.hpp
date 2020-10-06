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

#ifndef __Thea_Algorithms_OPTPPNumericalOptimizer_hpp__
#define __Thea_Algorithms_OPTPPNumericalOptimizer_hpp__

#include "OPTPPCommon.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/INumericalOptimizer.hpp"

namespace Thea {
namespace Algorithms {

/**
 * OPT++-based nonlinear optimizer.
 *
 * @see http://software.sandia.gov/opt++/
 */
class THEA_OPTPP_DLL_LOCAL OPTPPNumericalOptimizer : public INumericalOptimizer
{
  private:
    typedef INumericalOptimizer BaseType;

  public:
    /** Constructor. */
    OPTPPNumericalOptimizer(std::string const & name_);

    /**
     * @copydoc INumericalOptimizer::minimize()
     *
     * Valid options for the OPT++ backend are:
     * - <b>method</b>: Solution method to use
     *   - <i>Type:</i> <code>std::string</code> in {"BaNewton", "BaQNewton", "BCEllipsoid", "BCNewton", "BCQNewton", "CG",
     *     "DHNIPS", "FDNewton", "FDNIPS", "GSS", "LBFGS", "Newton", "NIPS", "PDS", "QNewton", "QNIPS"}
     *   - <i>Default</i>: "CG"
     *   - <i>Note</i>: See the OPT++ documentation to decide which method is suitable in which case. Some methods require
     *     the objective function to have analytic first and possibly second derivatives.
     *
     * @todo This function is currently a no-op.
     */
    bool THEA_ICALL minimize(IScalarFunction const & objective, double const * hint = nullptr, Options const & options = Options());

}; // class OPTPPNumericalOptimizer

/** Factory for creating OPT++ numerical optimizers. */
class THEA_OPTPP_DLL_LOCAL OPTPPNumericalOptimizerFactory : public INumericalOptimizerFactory
{
  public:
    /** Destructor. */
    ~OPTPPNumericalOptimizerFactory();

    INumericalOptimizer * THEA_ICALL createNumericalOptimizer(char const * name);
    void THEA_ICALL destroyNumericalOptimizer(INumericalOptimizer * optimizer);

    /** Destroy all numerical optimizers created with this factory. */
    void destroyAllNumericalOptimizers();

  private:
    typedef Set<INumericalOptimizer *> NumericalOptimizerSet;  ///< Set of numerical optimizers.

    NumericalOptimizerSet optimizers;  ///< All numerical optimizers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
