//============================================================================
//
// From the G3D project.Copyright (C) 2000-2008 Morgan McGuire.
//
// For the full G3D license see LICENSE.txt in the documentation.
//
//============================================================================

#include "GLCaps.hpp"
#include "../../FilePath.hpp"
#include "../../Map.hpp"
#include "../../Math.hpp"
#include <cstring>
#include <sstream>

#ifdef THEA_WINDOWS
#  include <windows.h>
#endif

namespace Thea {
namespace Graphics {
namespace GL {

#ifdef THEA_WINDOWS

namespace GLInternal {

bool registryKeyExists(std::string const & key);
bool registryValueExists(std::string const & key, std::string const & value);
bool registryReadString(std::string const & key, std::string const & value, std::string & data);

} // namespace GLInternal

#endif // THEA_WINDOWS

// Global init flags for GLCaps. Because this is an integer constant (equal to zero), we can safely assume that it will be
// initialized before this translation unit is entered.
bool GLCaps::_loadedExtensions   = false;
bool GLCaps::_initialized        = false;
bool GLCaps::_checkedForBugs     = false;
bool GLCaps::_hasGLMajorVersion2 = false;

bool GLCaps::has_headless_context = false;
GLContext GLCaps::headless_context;

#if defined(THEA_GL_OSMESA)
  char * GLCaps::headless_buffer = NULL;
#elif defined(THEA_WINDOWS)
#elif defined(THEA_LINUX) || defined(THEA_BSD)
  Display * GLCaps::headless_display = NULL;
  Pixmap GLCaps::headless_pixmap;
  GLXPixmap GLCaps::headless_glx_pixmap;
#elif defined(THEA_OSX)
#endif

int GLCaps::_numTextureCoords = 0;
int GLCaps::_numTextures      = 0;
int GLCaps::_numTextureUnits  = 0;

bool GLCaps::bug_glMultiTexCoord3fvARB = false;
bool GLCaps::bug_normalMapTexGen       = false;
bool GLCaps::bug_redBlueMipmapSwap     = false;
bool GLCaps::bug_mipmapGeneration      = false;
bool GLCaps::bug_slowVBO               = false;

// Cache of values supplied to supportsImageFormat. Works on pointers since there is no way for users to construct their own
// ImageFormats.
static TheaMap<Texture::Format const *, bool> _supportedTextureFormat;
static TheaMap<Texture::Format const *, bool> _supportedRenderBufferFormat;

GLCaps::Vendor
GLCaps::computeVendor()
{
  std::string s = vendor();

  if (s == "ATI Technologies Inc.")
    return ATI;
  else if (s == "NVIDIA Corporation")
    return NVIDIA;
  else if ((s == "Brian Paul") || (s == "Mesa project: www.mesa3d.org"))
    return MESA;
  else
    return ARB;
}

GLCaps::Vendor
GLCaps::enumVendor()
{
  return computeVendor();
}

#ifdef THEA_WINDOWS
// Used by the Windows version of getDriverVersion().
// Based on code by Ted Peck tpeck@roundwave.com http://www.codeproject.com/dll/ShowVer.asp
struct VS_VERSIONINFO
{
  WORD                wLength;
  WORD                wValueLength;
  WORD                wType;
  WCHAR               szKey[1];
  WORD                Padding1[1];
  VS_FIXEDFILEINFO    Value;
  WORD                Padding2[1];
  WORD                Children[1];
};
#endif

std::string
GLCaps::getDriverVersion()
{
  if (computeVendor() == MESA)
  {
    // Mesa includes the driver version in the renderer version e.g., "1.5 Mesa 6.4.2"
    static std::string _glVersion = (char *)glGetString(GL_VERSION);
    std::size_t i = _glVersion.rfind(' ');
    if (i == std::string::npos)
      return "Unknown (bad MESA driver string)";
    else
      return _glVersion.substr(i + 1, _glVersion.length() - i);
  }

#ifdef THEA_WINDOWS
  // Locate the driver on Windows and get the version this systems expects Windows 2000/XP/Vista
  std::string videoDriverKey;
  bool canCheckVideoDevice = GLInternal::registryKeyExists("HKEY_LOCAL_MACHINE\\HARDWARE\\DEVICEMAP\\VIDEO");
  if (canCheckVideoDevice)
  {
    // Find the driver expected to load
    std::string videoDeviceKey = "HKEY_LOCAL_MACHINE\\HARDWARE\\DEVICEMAP\\VIDEO";
    std::string videoDeviceValue = "\\Device\\Video";
    int videoDeviceNum = 0;

    while (GLInternal::registryValueExists(videoDeviceKey, format("%s%d", videoDeviceValue.c_str(), videoDeviceNum)))
      ++videoDeviceNum;

    // Find the key where the installed driver lives
    std::string installedDriversKey;
    GLInternal::registryReadString(videoDeviceKey, format("%s%d", videoDeviceValue.c_str(), videoDeviceNum - 1),
                                   installedDriversKey);

    // Find and remove the "\Registry\Machine\" part of the key
    std::size_t subKeyIndex = installedDriversKey.find('\\', 1);
    subKeyIndex = installedDriversKey.find('\\', subKeyIndex + 1);

    installedDriversKey.erase(0, subKeyIndex);

    // Read the list of driver files this display driver installed/loads. This is a multi-string value, but we only care about
    // the first entry so reading one string is fine
    std::string videoDrivers;
    GLInternal::registryReadString("HKEY_LOCAL_MACHINE" + installedDriversKey, "InstalledDisplayDrivers", videoDrivers);

    if (videoDrivers.find(',', 0) != std::string::npos)
      videoDrivers = videoDrivers.substr(0, videoDrivers.find(',', 0));

    char systemDirectory[512] = "";
    GetSystemDirectoryA(systemDirectory, sizeof(systemDirectory));

    std::string driverFileName = format("%s\\%s.dll", systemDirectory, videoDrivers.c_str());

    DWORD dummy;
    int size = GetFileVersionInfoSizeA((LPCSTR)driverFileName.c_str(), &dummy);
    if (size == 0)
      return "Unknown (Can't find driver)";

    void * buffer = new uint8[size];
    if (GetFileVersionInfoA((LPCSTR)driverFileName.c_str(), 0, size, buffer) == 0)
    {
      delete[] (uint8*)buffer;
      return "Unknown (Can't find driver)";
    }

    // Interpret the VS_VERSIONINFO header pseudo-struct
    VS_VERSIONINFO * pVS = (VS_VERSIONINFO *)buffer;
    debugAssertM(!wcscmp(pVS->szKey, L"VS_VERSION_INFO"), "GLCaps: Version info not found");

    uint8 * pVt = (uint8 *) &pVS->szKey[wcslen(pVS->szKey) + 1];

#   define roundoffs(a,b,r)   (((uint8 *)(b) - (uint8 *)(a) + ((r) - 1)) & ~((r) - 1))
#   define roundpos(b, a, r)  (((uint8 *)(a)) + roundoffs(a, b, r))

    VS_FIXEDFILEINFO * pValue = (VS_FIXEDFILEINFO *) roundpos(pVt, pVS, 4);

#   undef roundoffs
#   undef roundpos

    std::string result = "Unknown (Can't find driver)";

    if (pVS->wValueLength)
    {
      result = format("%d.%d.%d.%d",
                      pValue->dwProductVersionMS >> 16,
                      pValue->dwProductVersionMS & 0xFFFF,
                      pValue->dwProductVersionLS >> 16,
                      pValue->dwProductVersionLS & 0xFFFF);
    }

    delete[] (uint8 *)buffer;

    return result;
  }
  else
    return "Unknown (Can't find driver)";
#else
    return "Unknown";
#endif
}

void
GLCaps::init()
{
  if (!_initialized)
  {
    if (!glGetCurrentContext())
    {
      if (!createHeadlessContext())
      {
        throw Error("GLCaps: Error creating headless OpenGL context");
      }
    }

    loadExtensions();
    THEA_CHECK_GL_OK

    checkAllBugs();
    THEA_CHECK_GL_OK
  }
}

bool
GLCaps::createHeadlessContext()
{
  if (has_headless_context)
    return true;

#if defined(THEA_GL_OSMESA)

#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  // specify Z, stencil, accum sizes
  headless_context = OSMesaCreateContextExt(OSMESA_RGBA, 16, 0, 0, NULL);
#else
  headless_context = OSMesaCreateContext(OSMESA_RGBA, NULL);
#endif

  if (!headless_context)
  {
    THEA_ERROR << "GLCaps: Could not create OSMesa context";
    return false;
  }

  // Bind to a small dummy pixmap instead of binding to a window
  static int const DUMMY_FB_WIDTH   =  32;
  static int const DUMMY_FB_HEIGHT  =  32;
  headless_buffer = new char[DUMMY_FB_WIDTH * DUMMY_FB_HEIGHT];
  if (!headless_buffer)
  {
    THEA_ERROR << "GLCaps: Could not allocate dummy OSMesa offscreen buffer";
    return false;
  }

  if (!OSMesaMakeCurrent(headless_context, headless_buffer, GL_UNSIGNED_BYTE, DUMMY_FB_WIDTH, DUMMY_FB_HEIGHT))
  {
    THEA_ERROR << "GLCaps: Could not make new OSMesa context current";
    return false;
  }

  has_headless_context = true;

#elif defined(THEA_WINDOWS)

  // TODO

#elif defined(THEA_LINUX) || defined(THEA_BSD)

  headless_display = XOpenDisplay(0);
  if (!headless_display)
  {
    THEA_ERROR << "GLCaps: Could not open X display";
    return false;
  }

  static int attribs[] = {
    GLX_RGBA,
    GLX_RED_SIZE,    8,
    GLX_GREEN_SIZE,  8,
    GLX_BLUE_SIZE,   8,
    GLX_ALPHA_SIZE,  8,
    GLX_DEPTH_SIZE, 24,
    None
  };

  XVisualInfo * vi = glXChooseVisual(headless_display, DefaultScreen(headless_display), attribs);
  if (!vi)
  {
    THEA_ERROR << "GLCaps: Could not choose X visual";
    return false;
  }

  headless_context = glXCreateContext(headless_display, vi, 0, GL_TRUE);
  if (!headless_context)
  {
    THEA_ERROR << "GLCaps: Could not create X context";
    return false;
  }

  // Bind to a small dummy pixmap instead of binding to a window
  static int const DUMMY_FB_WIDTH   =  32;
  static int const DUMMY_FB_HEIGHT  =  32;
  headless_pixmap = XCreatePixmap(headless_display, DefaultRootWindow(headless_display), DUMMY_FB_WIDTH, DUMMY_FB_HEIGHT, 24);
  headless_glx_pixmap = glXCreateGLXPixmap(headless_display, vi, headless_pixmap);
  if (!glXMakeCurrent(headless_display, headless_glx_pixmap, headless_context))
  {
    THEA_ERROR << "GLCaps: Could not make new X context current";
    return false;
  }

  has_headless_context = true;

#elif defined(THEA_OSX)

  CGLPixelFormatAttribute attribs[] = {
    kCGLPFAColorSize,     (CGLPixelFormatAttribute)24,
    kCGLPFAAlphaSize,     (CGLPixelFormatAttribute) 8,
    kCGLPFADepthSize,     (CGLPixelFormatAttribute)24,
    kCGLPFAAccelerated,
    kCGLPFASampleBuffers, (CGLPixelFormatAttribute) 1,
    kCGLPFASamples,       (CGLPixelFormatAttribute) 4,
    kCGLPFAAllowOfflineRenderers,
    (CGLPixelFormatAttribute)0
  };

  CGLError err;
  CGLPixelFormatObj pix;
  GLint npix;
  if ((err = CGLChoosePixelFormat(attribs, &pix, &npix)) != kCGLNoError)
  {
    THEA_ERROR << "GLCaps: Could not choose CGL pixel format (error: " << CGLErrorString(err) << ')';
    return false;
  }

  err = CGLCreateContext(pix, 0, &headless_context);
  if (err != kCGLNoError || !headless_context)
  {
    THEA_ERROR << "GLCaps: Could not create CGL context (error: " << CGLErrorString(err) << ')';
    return false;
  }

  if ((err = CGLSetCurrentContext(headless_context)) != kCGLNoError)
  {
    THEA_ERROR << "GLCaps: Could not make new CGL context current (error: " << CGLErrorString(err) << ')';
    return false;
  }

  has_headless_context = true;

#endif

  return has_headless_context;
}

void
GLCaps::destroyHeadlessContext()
{
  if (!has_headless_context)
    return;

#if defined(THEA_GL_OSMESA)

  delete [] headless_buffer;
  OSMesaDestroyContext(headless_context);

#elif defined(THEA_WINDOWS)

  // TODO

#elif defined(THEA_LINUX) || defined(THEA_BSD)

  glXDestroyContext(headless_display, headless_context);
  glXDestroyPixmap(headless_display, headless_glx_pixmap);
  XFreePixmap(headless_display, headless_pixmap);

#elif defined(THEA_OSX)

  CGLDestroyContext(headless_context);

#endif

  has_headless_context = false;
}

bool
GLCaps::isHeadless()
{
  return has_headless_context && glGetCurrentContext() == headless_context;
}

void
GLCaps::loadExtensions()
{
  alwaysAssertM(glGetString(GL_RENDERER) != NULL, "Error initializing OpenGL, please check your GL installation");

  if (_loadedExtensions)
    return;

  alwaysAssertM(!_initialized, "Attempt to initialize OpenGL twice");
  alwaysAssertM(glGetCurrentContext(), "Unable to load OpenGL extensions without a current context.");

  GLenum err = glewInit();
  if (err != GLEW_OK)
    throw Error(format("Couldn't load extensions via GLEW (%s)", glewGetErrorString(err)));

  _loadedExtensions = true;

  // Initialize statically cached strings
  vendor();
  renderer();
  glVersion();
  driverVersion();

  // Initialize cached GL major version pulled from glVersion() for extensions made into 2.0 core
  std::string const glver = glVersion();
  _hasGLMajorVersion2 = beginsWith(glver, "2.");

  // GL_ARB_texture_cube_map doesn't work on Radeon Mobility
  // GL Renderer:    MOBILITY RADEON 9000 DDR x86/SSE2
  // GL Version:     1.3.4204 WinXP Release
  // Driver version: 6.14.10.6430

  // GL Vendor:      ATI Technologies Inc.
  // GL Renderer:    MOBILITY RADEON 7500 DDR x86/SSE2
  // GL Version:     1.3.3842 WinXP Release
  // Driver version: 6.14.10.6371

  if ((beginsWith(renderer(), "MOBILITY RADEON") || beginsWith(renderer(), "ATI MOBILITY RADEON"))
   && beginsWith(driverVersion(), "6.14.10.6"))
  {
      THEA_WARNING << "This ATI Radeon Mobility card has a known bug with cube maps. Put cube map texture coordinates in the"
                      " normals and use ARB_NORMAL_MAP to work around.";
  }

  if (THEA_GL_SUPPORTS(ARB_multitexture))
  {
    // Don't use more texture units than allowed at compile time
    GLint num_texunits;
    glGetIntegerv(GL_MAX_TEXTURE_UNITS_ARB, &num_texunits);
    _numTextureUnits = std::min(num_texunits, (int)MAX_TEXTURE_UNITS);
  }
  else
    _numTextureUnits = 1;

  // NVIDIA cards with GL_NV_fragment_program have different numbers of texture coords, units, and textures
  if (THEA_GL_SUPPORTS(NV_fragment_program))
  {
// To make this compile in the absence of the extension
#ifdef GL_MAX_TEXTURE_COORDS_NV
    glGetIntegerv(GL_MAX_TEXTURE_COORDS_NV, &_numTextureCoords);
    _numTextureCoords = Math::clamp(_numTextureCoords, _numTextureUnits, (int)MAX_TEXTURE_UNITS);

    glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS_NV, &_numTextures);
    _numTextures = Math::clamp(_numTextures, _numTextureUnits, (int)MAX_TEXTURE_UNITS);
#endif
  }
  else
  {
    _numTextureCoords = _numTextureUnits;
    _numTextures      = _numTextureUnits;
  }

  if (!THEA_GL_SUPPORTS(ARB_multitexture))
  {
    // No multitexture
    THEA_DEBUG << "No GL_ARB_multitexture support: forcing number of texture units to no more than 1";

    _numTextureCoords = std::max(1, _numTextureCoords);
    _numTextures      = std::max(1, _numTextures);
    _numTextureUnits  = std::max(1, _numTextureUnits);
  }

  THEA_CHECK_GL_OK

  _initialized = true;
}

void
GLCaps::checkAllBugs()
{
  if (_checkedForBugs)
    return;

  alwaysAssertM(_loadedExtensions, "Cannot check for OpenGL bugs before extensions are loaded.");

  bool is_headless = isHeadless();
  GLuint fb, color_tex, depth_rb;
  if (is_headless)
  {
    static GLsizei const FB_SIZE = 32;

    glGenTextures(1, &color_tex);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, FB_SIZE, FB_SIZE, 0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);
    THEA_CHECK_GL_OK

    glGenRenderbuffersEXT(1, &depth_rb);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depth_rb);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, FB_SIZE, FB_SIZE);
    THEA_CHECK_GL_OK

    glGenFramebuffersEXT(1, &fb);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, color_tex, 0);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depth_rb);
    THEA_CHECK_GL_OK

    GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
      throw Error(format("GLCaps: Could not create offscreen framebuffer for checking bugs (error code %d)", (int)status));

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fb);
    glViewport(0, 0, FB_SIZE, FB_SIZE);
    THEA_CHECK_GL_OK
  }

  checkBug_cubeMapBugs();
  checkBug_redBlueMipmapSwap();
  checkBug_mipmapGeneration();
  checkBug_slowVBO();

  if (is_headless)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteRenderbuffersEXT(1, &depth_rb);
    glDeleteFramebuffersEXT(1, &color_tex);
    glDeleteTextures(1, &color_tex);
  }

  _checkedForBugs = true;
}

bool
GLCaps::hasBug_glMultiTexCoord3fvARB()
{
  alwaysAssertM(_initialized, "GLCaps has not been initialized.");
  return bug_glMultiTexCoord3fvARB;
}

bool
GLCaps::hasBug_normalMapTexGen()
{
  alwaysAssertM(_initialized, "GLCaps has not been initialized.");
  return bug_normalMapTexGen;
}

bool
GLCaps::hasBug_redBlueMipmapSwap()
{
  alwaysAssertM(_initialized, "GLCaps has not been initialized.");
  return bug_redBlueMipmapSwap;
}

bool
GLCaps::hasBug_mipmapGeneration()
{
  alwaysAssertM(_initialized, "GLCaps has not been initialized.");
  return bug_mipmapGeneration;
}

bool GLCaps::hasBug_slowVBO()
{
  alwaysAssertM(_initialized, "GLCaps has not been initialized.");
  return bug_slowVBO;
}

bool
GLCaps::supports(std::string const & extensions)
{
  return glewIsSupported(extensions.c_str()) != GL_FALSE;
}

bool
GLCaps::supportsTexture(Texture::Format const * fmt)
{
  // First, check if we've already tested this format
  TheaMap<Texture::Format const *, bool>::iterator existing = _supportedTextureFormat.find(fmt);
  if (existing == _supportedTextureFormat.end())
  {
    bool supportsFormat = false;

    if (!(fmt->floatingPoint && !THEA_GL_SUPPORTS(ARB_texture_float)))
    {
      // Allocate some space for making a dummy texture
      uint8 bytes[8 * 8 * 4];

      glPushAttrib(GL_TEXTURE_BIT);
      {
        // Clear the error bit
        glGetError();

        // See if we can create a texture in this format
        unsigned int id;
        glGenTextures(1, &id);
        glBindTexture(GL_TEXTURE_2D, id);

        // Clear the old error flag
        glGetError();
        // 2D texture, level of detail 0 (normal), internal format, x size from image, y size from image,
        // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
        glTexImage2D(GL_TEXTURE_2D, 0, fmt->openGLFormat, 8, 8, 0, fmt->openGLBaseFormat, GL_UNSIGNED_BYTE, bytes);

        supportsFormat = (glGetError() == GL_NO_ERROR);

        glBindTexture(GL_TEXTURE_2D, 0);
        glDeleteTextures(1, &id);
      }
      glPopAttrib();
    }

    _supportedTextureFormat[fmt] = supportsFormat;
    return supportsFormat;
  }
  else
    return existing->second;
}

bool
GLCaps::supportsRenderBuffer(const Texture::Format* fmt)
{
  // First, check if we've already tested this format
  TheaMap<Texture::Format const *, bool>::iterator existing = _supportedRenderBufferFormat.find(fmt);
  if (existing == _supportedRenderBufferFormat.end())
  {
    bool supportsFormat = false;

    if (THEA_GL_SUPPORTS(EXT_framebuffer_object) && !(fmt->floatingPoint && !THEA_GL_SUPPORTS(ARB_texture_float)))
    {
      glPushAttrib(GL_COLOR_BUFFER_BIT);
      {
        // Clear the error bit
        glGetError();

        // See if we can create a render buffer in this format
        unsigned int id;
        glGenRenderbuffersEXT (1, &id);

        // Clear the old error flag
        glGetError();

        glBindRenderbufferEXT (GL_RENDERBUFFER_EXT, id);
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, fmt->openGLFormat, 8, 8);

        supportsFormat = (glGetError() == GL_NO_ERROR);

        glBindRenderbufferEXT (GL_RENDERBUFFER_EXT, 0);
        glDeleteRenderbuffersEXT(1, &id);
      }
      glPopAttrib();
    }

    _supportedRenderBufferFormat[fmt] = supportsFormat;
    return supportsFormat;
  }
  else
    return existing->second;
}

std::string const &
GLCaps::glVersion()
{
  alwaysAssertM(_loadedExtensions, "Cannot call GLCaps::glVersion before GLCaps::init().");
  static std::string _glVersion = (char *)glGetString(GL_VERSION);
  return _glVersion;
}

std::string const &
GLCaps::driverVersion()
{
  alwaysAssertM(_loadedExtensions, "Cannot call GLCaps::driverVersion before GLCaps::init().");
  static std::string _driverVersion = getDriverVersion().c_str();
  return _driverVersion;
}

std::string const &
GLCaps::vendor()
{
  alwaysAssertM(_loadedExtensions, "Cannot call GLCaps::vendor before GLCaps::init().");
  static std::string _driverVendor = (char *)glGetString(GL_VENDOR);
  return _driverVendor;
}

std::string const & GLCaps::renderer()
{
  alwaysAssertM(_loadedExtensions, "Cannot call GLCaps::renderer before GLCaps::init().");
  static std::string _glRenderer = (char *)glGetString(GL_RENDERER);
  return _glRenderer;
}

bool
GLCaps::supportsTwoSidedStencil()
{
  return THEA_GL_SUPPORTS(ATI_separate_stencil) || THEA_GL_SUPPORTS(EXT_stencil_two_side);
}

void
GLCaps::checkBug_cubeMapBugs()
{
  if (!THEA_GL_SUPPORTS(ARB_texture_cube_map))
  {
    // No cube map == no bug
    bug_glMultiTexCoord3fvARB = false;
    bug_normalMapTexGen = false;
    return;
  }

  // Save current GL state
  unsigned int id;
  glGenTextures(1, &id);

  THEA_CHECK_GL_OK

  glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (!isHeadless())
    {
      glDrawBuffer(GL_FRONT);
      glReadBuffer(GL_FRONT);

      // Normally, GL_FRONT is guaranteed to be present regardless of what sort of output device we're using. Except in Qt 5,
      // QOpenGLWidget does not have a front buffer and renders everything offscreen.
      if (glGetError() != GL_NO_ERROR)
      {
        glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
        glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
      }
    }

    glClearColor(0, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);

    THEA_CHECK_GL_OK

    GLenum target[] = {
        GL_TEXTURE_CUBE_MAP_POSITIVE_X_ARB,
        GL_TEXTURE_CUBE_MAP_NEGATIVE_X_ARB,
        GL_TEXTURE_CUBE_MAP_POSITIVE_Y_ARB,
        GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_ARB,
        GL_TEXTURE_CUBE_MAP_POSITIVE_Z_ARB,
        GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_ARB};

    // Face colors
    unsigned char color[6];

    // Create a cube map
    glActiveTextureARB(GL_TEXTURE0_ARB);
    glBindTexture(GL_TEXTURE_CUBE_MAP_ARB, id);
    glEnable(GL_TEXTURE_CUBE_MAP_ARB);

    THEA_CHECK_GL_OK

    {
      static int const N = 16;
      unsigned int image[N * N];
      for (int f = 0; f < 6; ++f)
      {
        color[f] = f * 40;

        // Fill each face with a different color
        std::memset(image, color[f], N * N * sizeof(unsigned int));

        // 2D texture, level of detail 0 (normal), internal format, x size from image, y size from image,
        // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
        glTexImage2D(target[f], 0, GL_RGBA, N, N, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
      }

      glTexParameteri(GL_TEXTURE_CUBE_MAP_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_CUBE_MAP_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_CUBE_MAP_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_CUBE_MAP_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_CUBE_MAP_ARB, GL_TEXTURE_WRAP_R, GL_CLAMP);
    }

    THEA_CHECK_GL_OK

    // Set orthogonal projection
    float viewport[4];
    glGetFloatv(GL_VIEWPORT, viewport);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(viewport[0], viewport[0] + viewport[2], viewport[1] + viewport[3], viewport[1], -1, 10);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glActiveTextureARB(GL_TEXTURE0_ARB + 0);
    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();

    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
    glDisable(GL_TEXTURE_GEN_R);

    // Render one sample from each cube map face
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glColor3f(1, 1, 1);

    static float const corner[] = {
       1.000000, -1.000000,  1.000000,
       1.000000, -1.000000, -1.000000,
       1.000000,  1.000000, -1.000000,
       1.000000,  1.000000,  1.000000,

      -1.000000,  1.000000,  1.000000,
      -1.000000,  1.000000, -1.000000,
      -1.000000, -1.000000, -1.000000,
      -1.000000, -1.000000,  1.000000,

       1.000000,  1.000000,  1.000000,
       1.000000,  1.000000, -1.000000,
      -1.000000,  1.000000, -1.000000,
      -1.000000,  1.000000,  1.000000,

       1.000000, -1.000000,  1.000000,
      -1.000000, -1.000000,  1.000000,
      -1.000000, -1.000000, -1.000000,
       1.000000, -1.000000, -1.000000,

      -1.000000, -1.000000,  1.000000,
       1.000000, -1.000000,  1.000000,
       1.000000,  1.000000,  1.000000,
      -1.000000,  1.000000,  1.000000,

      -1.000000,  1.000000, -1.000000,
       1.000000,  1.000000, -1.000000,
       1.000000, -1.000000, -1.000000,
      -1.000000, -1.000000, -1.000000
    };

    for (int i = 0; i < 2; ++i)
    {
      // First time through, use multitex coord
      if (i == 1)
      {
        // Second time through, use normal map generation
        glActiveTextureARB(GL_TEXTURE0_ARB + 0);
        glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP_ARB);
        glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP_ARB);
        glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_NORMAL_MAP_ARB);
        glEnable(GL_TEXTURE_GEN_S);
        glEnable(GL_TEXTURE_GEN_T);
        glEnable(GL_TEXTURE_GEN_R);
      }

      glBegin(GL_QUADS);
        for (int f = 0; f < 6; ++f) {
          const float s = 10.0f;

          glMultiTexCoord3fvARB(GL_TEXTURE0_ARB, corner + 12 * f + 0);
          glNormal3fv(corner + 12 * f + 0);
          glVertex4f(f * s, 0, -1, 1);

          glMultiTexCoord3fvARB(GL_TEXTURE0_ARB, corner + 12 * f + 3);
          glNormal3fv(corner + 12 * f + 3);
          glVertex4f(f * s, s, -1, 1);

          glMultiTexCoord3fvARB(GL_TEXTURE0_ARB, corner + 12 * f + 6);
          glNormal3fv(corner + 12 * f + 6);
          glVertex4f((f + 1) * s, s, -1, 1);

          glMultiTexCoord3fvARB(GL_TEXTURE0_ARB, corner + 12 * f + 9);
          glNormal3fv(corner + 12 * f + 9);
          glVertex4f((f + 1) * s, 0, -1, 1);
        }
      glEnd();

      // Read back results
      unsigned int readback[60];
      glReadPixels(0, (int)(viewport[3] - 5), 60, 1, GL_RGBA, GL_UNSIGNED_BYTE, readback);

      // Test result for errors
      bool texbug = false;
      for (int f = 0; f < 6; ++f)
      {
        if ((readback[f * 10 + 5] & 0xFF) != color[f])
        {
          texbug = true;
          break;
        }
      }

      if (i == 0)
          bug_glMultiTexCoord3fvARB = texbug;
      else
          bug_normalMapTexGen = texbug;

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

  glPopAttrib();

  glDeleteTextures(1, &id);

  THEA_CHECK_GL_OK
}

void
GLCaps::checkBug_redBlueMipmapSwap()
{
  if (!THEA_GL_SUPPORTS(SGIS_generate_mipmap))
  {
    // No auto-mipmaps == no bug
    bug_redBlueMipmapSwap = false;
    return;
  }

  glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLuint id;
    glGenTextures(1, &id);

    glBindTexture(GL_TEXTURE_2D, id);

    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);

    int N = 4 * 4 * 3;
    uint8 * bytes = new uint8[N];

    std::memset(bytes, 0, 4 *4 *3);
    for (int i = 0; i < N; i += 3)
      bytes[i] = 0xFF;

    // 2D texture, level of detail 0 (normal), internal format, x size from image, y size from image, border 0 (normal), rgb
    // color data, unsigned byte data, and finally the data itself.
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, 4, 4, 0, GL_RGB, GL_UNSIGNED_BYTE, bytes);

    // Read the data back.
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, bytes);

    // Verify that the data made the round trip correctly
    bug_redBlueMipmapSwap = (bytes[0] != 0xFF) || (bytes[1] != 0x00) || (bytes[2] != 0x00);

    delete[] bytes;
    glDeleteTextures(1, &id);
  glPopAttrib();

  THEA_CHECK_GL_OK
}

void
GLCaps::checkBug_mipmapGeneration()
{
  std::string const & r = renderer();

  // The mip-maps are arbitrarily corrupted; we have not yet generated a reliable test for this case.

  bug_mipmapGeneration = THEA_GL_SUPPORTS(SGIS_generate_mipmap) && (beginsWith(r, "MOBILITY RADEON 90")
                                                                     || beginsWith(r, "MOBILITY RADEON 57")
                                                                     || beginsWith(r, "Intel 845G")
                                                                     || beginsWith(r, "Intel 854G"));

  THEA_CHECK_GL_OK
}

void
GLCaps::checkBug_slowVBO()
{
  if (!THEA_GL_SUPPORTS(ARB_vertex_buffer_object))
  {
    // No VBO == no bug
    bug_slowVBO = false;
    return;
  }

  std::string const & r = renderer();
  bug_slowVBO = beginsWith(r, "MOBILITY RADEON 7500");

  return;

  // Note: Test incomplete.
}

void
GLCaps::describeSystem(std::ostream & os)
{
  GLint max_tex_size;
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_tex_size);

  os << "GPU {\n"
     << "  Chipset              =  " << renderer() << '\n'
     << "  Vendor               =  " << vendor() << '\n'
     << "  Driver               =  " << driverVersion() << '\n'
     << "  OpenGL version       =  " << glVersion() << '\n'
     << "  Textures             =  " << numTextures() << '\n'
     << "  Texture coordinates  =  " << numTextureCoords() << '\n'
     << "  Texture units        =  " << numTextureUnits() << '\n'
     << "  Texture size         =  " << max_tex_size << '\n'
     << "}" << std::endl;
}

std::string
GLCaps::describeSystem()
{
  std::ostringstream os;
  describeSystem(os);
  return os.str();
}

#ifdef THEA_WINDOWS

namespace GLInternal {

HKEY
getRootKeyFromString(char const * str, size_t length)
{
  if (str)
  {
    if       ( std::strncmp(str, "HKEY_CLASSES_ROOT", length) == 0 )
      return HKEY_CLASSES_ROOT;
    else if  ( std::strncmp(str, "HKEY_CURRENT_CONFIG", length) == 0 )
      return HKEY_CURRENT_CONFIG;
    else if  ( std::strncmp(str, "HKEY_CURRENT_USER", length) == 0 )
      return HKEY_CURRENT_USER;
    else if  ( std::strncmp(str, "HKEY_LOCAL_MACHINE", length) == 0 )
      return HKEY_LOCAL_MACHINE;
    else if  ( std::strncmp(str, "HKEY_PERFORMANCE_DATA", length) == 0 )
      return HKEY_PERFORMANCE_DATA;
    else if  ( std::strncmp(str, "HKEY_PERFORMANCE_NLSTEXT", length) == 0 )
      return HKEY_PERFORMANCE_NLSTEXT;
    else if  ( std::strncmp(str, "HKEY_PERFORMANCE_TEXT", length) == 0 )
      return HKEY_PERFORMANCE_TEXT;
    else if  ( std::strncmp(str, "HKEY_CLASSES_ROOT", length) == 0 )
      return HKEY_CLASSES_ROOT;
    else
      return NULL;
  }
  else
    return NULL;
}

bool
registryKeyExists(std::string const & key)
{
  size_t pos = key.find('\\', 0);
  if (pos == std::string::npos)
    return false;

  HKEY hkey = getRootKeyFromString(key.c_str(), pos);
  if (hkey == NULL)
    return false;

  HKEY openKey;
  int32 result = RegOpenKeyExA(hkey, (key.c_str() + pos + 1), 0, KEY_READ, &openKey);
  debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't open registry key");

  if (result == ERROR_SUCCESS)
  {
    RegCloseKey(openKey);
    return true;
  }
  else
    return false;
}

bool
registryValueExists(std::string const & key, std::string const & value)
{
  size_t pos = key.find('\\', 0);
  if (pos == std::string::npos)
    return false;

  HKEY hkey = getRootKeyFromString(key.c_str(), pos);
  if ( hkey == NULL )
    return false;

  HKEY openKey;
  int32 result = RegOpenKeyExA(hkey, (key.c_str() + pos + 1), 0, KEY_READ, &openKey);
  debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't open registry key");

  if (result == ERROR_SUCCESS)
  {
    uint32 dataSize = 0;
    result = RegQueryValueExA(openKey, value.c_str(), NULL, NULL, NULL, reinterpret_cast<LPDWORD>(&dataSize));
    debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't query registry value");
    RegCloseKey(openKey);
  }

  return (result == ERROR_SUCCESS);
}

bool
registryReadString(std::string const & key, std::string const & value, std::string & data)
{
  size_t pos = key.find('\\', 0);
  if (pos == std::string::npos)
    return false;

  HKEY hkey = getRootKeyFromString(key.c_str(), pos);
  if (hkey == NULL)
    return false;

  HKEY openKey;
  int32 result = RegOpenKeyExA(hkey, (key.c_str() + pos + 1), 0, KEY_READ, &openKey);
  debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't open registry key");

  if (result == ERROR_SUCCESS)
  {
    uint32 dataSize = 0;
    result = RegQueryValueExA(openKey, value.c_str(), NULL, NULL, NULL, reinterpret_cast<LPDWORD>(&dataSize));
    debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't query registry value");

    if (result == ERROR_SUCCESS)
    {
      // increment datasize to allow for non null-terminated strings in registry
      dataSize += 1;

      char * tmpStr = static_cast<char *>(std::malloc(dataSize));
      std::memset(tmpStr, 0, (size_t)dataSize);
      result = RegQueryValueExA(openKey, value.c_str(), NULL, NULL, reinterpret_cast<LPBYTE>(tmpStr),
                                reinterpret_cast<LPDWORD>(&dataSize));
      debugAssertM(result == ERROR_SUCCESS || result == ERROR_FILE_NOT_FOUND, "Couldn't query registry value");

      if (result == ERROR_SUCCESS)
        data = tmpStr;

      RegCloseKey(openKey);
      std::free(tmpStr);
    }
  }

  return (result == ERROR_SUCCESS);
}

} // namespace GLInternal

#endif // THEA_WINDOWS

} // namespace GL
} // namespace Graphics
} // namespace Thea
