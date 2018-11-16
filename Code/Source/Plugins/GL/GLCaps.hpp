//============================================================================
//
// From the G3D project. Copyright (C) 2000-2008 Morgan McGuire.
//
// For the full G3D license see LICENSE.txt in the documentation.
//
//============================================================================

#ifndef __Thea_Graphics_GL_GLCaps_hpp__
#define __Thea_Graphics_GL_GLCaps_hpp__

#include "../../Graphics/Texture.hpp"
#include "GLCommon.hpp"
#include "GLHeaders.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

#define THEA_GL_SUPPORTS(extension)  (GLEW_##extension)
#define THEA_WGL_SUPPORTS(extension) (WGLEW_##extension)
#define THEA_GLX_SUPPORTS(extension) (GLXEW_##extension)

/**
 * Maintains information about the OpenGL installation, and loads extensions on startup.
 *
 * The <code>hasBug_</code> methods detect specific common bugs on graphics cards.  They can be used to switch to fallback
 * rendering paths.
 *
 * This class is absolutely <b>NOT threadsafe</b>.
 */
class THEA_GL_DLL_LOCAL GLCaps
{
  public:
    /** Common OpenGL vendors. */
    enum Vendor { ATI, NVIDIA, MESA, ARB };

    /** Maximum number of texture coordinates supported by GLRenderSystem; used to preallocate some static arrays. */
    enum { MAX_TEXTURE_UNITS = 8 };

    /**
     * Loads OpenGL extensions (such as glBindBufferARB). Call this once at the beginning of the program, <b>after</b> a video
     * device is created. This is called for you if you use GLRenderSystem. It is safe to call this function twice as long
     * as the calls are from the same thread.
     */
    static void init();

    /** Is the headless rendering context, if any, currently active? */
    static bool isHeadless();

    /**
     * Check if the OpenGL implementation supports a set of extensions. The names of the extensions should be separated by
     * spaces and include the initial GL/WGL/GLX... prefix, e.g. <code>"GL_ARB_multitexture WGLEW_ARB_pbuffer"</code>.
     */
    static bool supports(std::string const & extensions);

    /** Returns true if the given texture format is supported on this device for Textures. */
    static bool supportsTexture(Texture::Format const * fmt);

    /** Returns true if the given texture format is supported on this device for RenderBuffers. */
    static bool supportsRenderBuffer(Texture::Format const * fmt);

    /** Get the OpenGL version. */
    static std::string const & glVersion();

    /** Get the driver version. */
    static std::string const & driverVersion();

    /** Get the name of the OpenGL vendor. */
    static std::string const & vendor();

    /** Get the enum representing the OpenGL vendor. */
    static Vendor enumVendor();

    /** Get the name of the rendering device. */
    static std::string const & renderer();

    /** Check if either GL_EXT_stencil_two_side or GL_ATI_separate_stencil is supported. */
    static bool supportsTwoSidedStencil();

    /** Get the number of texture coordinate sets. */
    static int numTextureCoords() { return _numTextureCoords; }

    /** Some devices (such as NVIDIA cards) support more textures than texture matrices */
    static int numTextures() { return _numTextures; }

    /** Get the number of available texture units. */
    static int numTextureUnits() { return _numTextureUnits; }

    /**
     * Returns true if cube map support has a specific known bug on this card. Returns false otherwise, or if cube maps are not
     * supported at all on this card.
     *
     * Call after OpenGL is intialized.  Will render on the backbuffer but not make the backbuffer visible. Safe to call
     * multiple times; the result is memoized.
     *
     * On some Radeon Mobility cards (up to Mobility 9200), glMultiTexCoord3fvARB and glVertex4fv together create incorrect
     * texture lookups from cube maps. Using glVertex3fv or glTexCoord with glActiveTextureARB avoids this problem, as does
     * using normal map generation.
     */
    static bool hasBug_glMultiTexCoord3fvARB();

    /**
     * Returns true if cube map support has a specific known bug on this card that prevents correct normal map coordinate
     * generation. Returns false otherwise, or if cube maps are not supported at all on this card.
     *
     * Call after OpenGL is intialized.  Will render on the backbuffer but not make the backbuffer visible. Safe to call
     * multiple times; the result is memoized.
     *
     * Radeon Mobility 7500 has been shown to have a bug where not only does hasBug_glMultiTexCoord3fvARB() exist, but normal
     * maps can't work around the problem.
     */
    static bool hasBug_normalMapTexGen();

    /**
     * Check if the red and blue channels are occasionally flipped when auto-generating mipmaps. The Radeon Mobility 7500 is
     * known to have this bug. This has proven to be a reliable test.
     *
     * If this bug is detected, GLTexture switches to RGBA8 formats for RGB8 data.
     */
    static bool hasBug_redBlueMipmapSwap();

    /**
     * Returns true if SGIS auto mip-map generation occasionally produces buggy results (usually, pieces of other textures in
     * the low-level mipmaps).
     *
     * Radeon Mobility 9200 has this bug for some drivers.
     *
     * If this bug is detected, GLTexture reverts to software mipmap generation.
     */
    static bool hasBug_mipmapGeneration();

    /**
     * Some graphics cards (e.g. Radeon Mobility 7500) support the VBO extension but it is slower than main memory in most cases
     * due to poor cache behavior. This method performs a speed test the first time it is invoked and identifies those cards.
     */
    static bool hasBug_slowVBO();

    /** Print a (multi-line) description of the OpenGL system to an output stream. */
    static void describeSystem(std::ostream & os);

    /** Get a description of the OpenGL system as a (multi-line) string. */
    static std::string describeSystem();

  private:
    /** Create a headless rendering context, if possible. */
    static bool createHeadlessContext();

    /** Destroy the headless rendering context, if any. */
    static void destroyHeadlessContext();

    /** Load all OpenGL extensions and get system info. */
    static void loadExtensions();

    /** Detect the OpenGL vendor. */
    static Vendor computeVendor();

    /**
     * Returns the version string for the video driver.
     *
     * @cite Based in part on code by Ted Peck tpeck@roundwave.com http://www.codeproject.com/dll/ShowVer.asp
     */
    static std::string getDriverVersion();

    // Tests for hasBug_glMultiTexCoord3fvARB and hasBug_glNormalMapTexGenARB
    static void checkBug_cubeMapBugs();
    static void checkBug_redBlueMipmapSwap();
    static void checkBug_mipmapGeneration();
    static void checkBug_slowVBO();

    /** Runs all of the checkBug_ methods. Called from loadExtensions(). */
    static void checkAllBugs();

    /** True when init has been called. */
    static bool _initialized;

    static bool has_headless_context;  ///< Has a headless context been created?
    static GLContext headless_context;  ///< Handle to headless context, if any.

#if defined(THEA_GL_OSMESA)
    static char * headless_buffer;
#elif defined(THEA_WINDOWS)
#elif defined(THEA_LINUX) || defined(THEA_BSD)
    static Display * headless_display;
    static Pixmap headless_pixmap;
    static GLXPixmap headless_glx_pixmap;
#elif defined(THEA_MAC)
#endif

    /** True when loadExtensions has already been called. */
    static bool _loadedExtensions;

    /** True if this is GL 2.0 or greater, which mandates certain extensions.*/
    static bool _hasGLMajorVersion2;

    static int  _numTextureCoords;
    static int  _numTextures;
    static int  _numTextureUnits;

    /** True when checkAllBugs has been called. */
    static bool _checkedForBugs;

    static bool bug_glMultiTexCoord3fvARB;
    static bool bug_normalMapTexGen;
    static bool bug_redBlueMipmapSwap;
    static bool bug_mipmapGeneration;
    static bool bug_slowVBO;

};  // class GLCaps

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
