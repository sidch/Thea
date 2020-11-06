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

#ifndef __Thea_Graphics_Gl_GlPlugin_hpp__
#define __Thea_Graphics_Gl_GlPlugin_hpp__

#include "../../IPlugin.hpp"
#include "GlCommon.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

// Forward declaration
class GlRenderSystemFactory;

/** An OpenGL rendering plugin. */
class THEA_GL_DLL_LOCAL GlPlugin : public virtual IPlugin
{
  public:
    /** Constructor. */
    GlPlugin(IFactoryRegistry * registry_);

    /** Destructor. */
    ~GlPlugin();

    char const * THEA_ICALL getName() const;
    int8 THEA_ICALL setName(char const * s) { return false;  /* name is read-only */ }

    void THEA_ICALL install();
    void THEA_ICALL startup();
    void THEA_ICALL shutdown();
    void THEA_ICALL uninstall();

  private:
    IFactoryRegistry * registry;
    GlRenderSystemFactory * factory;
    bool started;

}; // class GlPlugin

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
