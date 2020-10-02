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

#ifndef __Thea_Algorithms_CsparsePlugin_hpp__
#define __Thea_Algorithms_CsparsePlugin_hpp__

#include "../../IPlugin.hpp"
#include "CsparseCommon.hpp"

namespace Thea {
namespace Algorithms {

// Forward declaration
class CsparseLinearSolverFactory;

/** A CSPARSE-based plugin for solving sparse systems of linear equations. */
class THEA_CSPARSE_DLL_LOCAL CsparsePlugin : public IPlugin
{
  public:
    /** Default constructor. */
    CsparsePlugin(IFactoryRegistry * registry_);

    /** Destructor. */
    ~CsparsePlugin();

    char const * getName() const;
    void install();
    void startup();
    void shutdown();
    void uninstall();

  private:
    IFactoryRegistry * registry;
    CsparseLinearSolverFactory * factory;
    bool started;

}; // class CsparsePlugin

} // namespace Algorithms
} // namespace Thea

#endif
