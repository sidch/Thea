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

#ifndef __Thea_Graphics_GLPlugin_hpp__
#define __Thea_Graphics_GLPlugin_hpp__

#include "../../Plugin.hpp"
#include "GLCommon.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Forward declaration
class GLRenderSystemFactory;

/** An OpenGL rendering plugin. */
class THEA_GL_DLL_LOCAL GLPlugin : public Plugin
{
  public:
    /** Constructor. */
    GLPlugin(FactoryRegistry * registry_);

    /** Destructor. */
    ~GLPlugin();

    char const * getName() const;
    void install();
    void startup();
    void shutdown();
    void uninstall();

  private:
    FactoryRegistry * registry;
    GLRenderSystemFactory * factory;
    bool started;

}; // class GLPlugin

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
