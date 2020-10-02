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

#ifndef __Thea_Algorithms_ArpackPlugin_hpp__
#define __Thea_Algorithms_ArpackPlugin_hpp__

#include "../../IPlugin.hpp"
#include "ArpackCommon.hpp"

namespace Thea {
namespace Algorithms {

// Forward declaration
class ArpackEigenSolverFactory;

/** An ARPACK-based eigensolver plugin. */
class THEA_ARPACK_DLL_LOCAL ArpackPlugin : public IPlugin
{
  public:
    /** Constructor. */
    ArpackPlugin(IFactoryRegistry * registry_);

    /** Destructor. */
    ~ArpackPlugin();

    char const * getName() const;
    void install();
    void startup();
    void shutdown();
    void uninstall();

  private:
    IFactoryRegistry * registry;
    ArpackEigenSolverFactory * factory;
    bool started;

}; // class ArpackPlugin

} // namespace Algorithms
} // namespace Thea

#endif
