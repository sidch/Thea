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

#ifndef __Thea_Algorithms_CSPARSEPlugin_hpp__
#define __Thea_Algorithms_CSPARSEPlugin_hpp__

#include "../../Plugin.hpp"
#include "CSPARSECommon.hpp"

namespace Thea {
namespace Algorithms {

// Forward declaration
class CSPARSELinearSolverFactory;

/** A CSPARSE-based plugin for solving sparse systems of linear equations. */
class THEA_CSPARSE_DLL_LOCAL CSPARSEPlugin : public Plugin
{
  public:
    /** Default constructor. */
    CSPARSEPlugin(FactoryRegistry * registry_);

    /** Destructor. */
    ~CSPARSEPlugin();

    char const * getName() const;
    void install();
    void startup();
    void shutdown();
    void uninstall();

  private:
    FactoryRegistry * registry;
    CSPARSELinearSolverFactory * factory;
    bool started;

}; // class CSPARSEPlugin

} // namespace Algorithms
} // namespace Thea

#endif
