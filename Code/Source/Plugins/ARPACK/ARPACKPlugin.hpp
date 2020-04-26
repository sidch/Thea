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

#ifndef __Thea_Algorithms_ARPACKPlugin_hpp__
#define __Thea_Algorithms_ARPACKPlugin_hpp__

#include "../../Plugin.hpp"
#include "ARPACKCommon.hpp"

namespace Thea {
namespace Algorithms {

// Forward declaration
class ARPACKEigenSolverFactory;

/** An ARPACK-based eigensolver plugin. */
class THEA_ARPACK_DLL_LOCAL ARPACKPlugin : public Plugin
{
  public:
    /** Constructor. */
    ARPACKPlugin(FactoryRegistry * registry_);

    /** Destructor. */
    ~ARPACKPlugin();

    char const * getName() const;
    void install();
    void startup();
    void shutdown();
    void uninstall();

  private:
    FactoryRegistry * registry;
    ARPACKEigenSolverFactory * factory;
    bool started;

}; // class ARPACKPlugin

} // namespace Algorithms
} // namespace Thea

#endif
