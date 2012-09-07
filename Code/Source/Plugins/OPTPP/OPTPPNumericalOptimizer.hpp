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

#ifndef __Thea_Algorithms_OPTPPNumericalOptimizer_hpp__
#define __Thea_Algorithms_OPTPPNumericalOptimizer_hpp__

#include "OPTPPCommon.hpp"
#include "../../Set.hpp"
#include "../../Algorithms/NumericalOptimizer.hpp"

namespace Thea {
namespace Algorithms {

/**
 * OPT++-based nonlinear optimizer.
 *
 * @see http://software.sandia.gov/opt++/
 */
class THEA_OPTPP_DLL_LOCAL OPTPPNumericalOptimizer : public NumericalOptimizer
{
  private:
    typedef NumericalOptimizer BaseType;

  public:
    /** Constructor. */
    OPTPPNumericalOptimizer(std::string const & name_);

    /**
     * {@inheritDoc}
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
    bool minimize(ScalarFunction const & objective, double const * hint = NULL, Options const & options = Options());

}; // class OPTPPNumericalOptimizer

/** Factory for creating OPT++ numerical optimizers. */
class THEA_OPTPP_DLL_LOCAL OPTPPNumericalOptimizerFactory : public NumericalOptimizerFactory
{
  public:
    /** Destructor. */
    ~OPTPPNumericalOptimizerFactory();

    NumericalOptimizer * createNumericalOptimizer(std::string const & name);
    void destroyNumericalOptimizer(NumericalOptimizer * optimizer);

    /** Destroy all numerical optimizers created with this factory. */
    void destroyAllNumericalOptimizers();

  private:
    typedef TheaSet<NumericalOptimizer *> NumericalOptimizerSet;  ///< Set of numerical optimizers.

    NumericalOptimizerSet optimizers;  ///< All numerical optimizers created by this factory.
};

} // namespace Algorithms
} // namespace Thea

#endif
