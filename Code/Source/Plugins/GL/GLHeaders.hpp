//============================================================================
//
// From the G3D project. Copyright (C) 2000-2008 Morgan McGuire.
//
// For the full G3D license see LICENSE.txt in the documentation.
//
//============================================================================

#ifndef __Thea_Graphics_GL_GLHeaders_hpp__
#define __Thea_Graphics_GL_GLHeaders_hpp__

#include "GLCommon.hpp"
#include "../../FilePath.hpp"

// Main includes

#ifndef GLEW_STATIC
#  define GLEW_STATIC
#endif

#include "glew.h"

#if defined(THEA_GL_OSMESA)

// glew.h undefs a bunch of visibility macros after it's done with them, but we need to restore them for osmesa.h
#  if defined(_WIN32)
#    ifndef APIENTRY
#      if defined(__MINGW32__)
#        define APIENTRY __stdcall
#      elif (_MSC_VER >= 800) || defined(_STDCALL_SUPPORTED) || defined(__BORLANDC__)
#        define APIENTRY __stdcall
#      else
#        define APIENTRY
#      endif
#    endif
#    ifndef GLAPI
#      if defined(__MINGW32__)
#        define GLAPI extern
#      endif
#    endif
     /* <wingdi.h> and <winnt.h> */
#    ifndef WINGDIAPI
#      define WINGDIAPI __declspec(dllimport)
#    endif
#    ifndef GLAPI
#      if defined(__MINGW32__)
#        define GLAPI extern
#      else
#        define GLAPI WINGDIAPI
#      endif
#    endif

#    ifndef GLAPIENTRY
#      define GLAPIENTRY APIENTRY
#    endif

#  else /* _UNIX */

#    define APIENTRY

     /* <glu.h> */
#    ifndef GLAPI
#      define GLAPI extern
#    endif
#    ifndef GLAPIENTRY
#      define GLAPIENTRY
#    endif

#  endif /* _WIN32 */

#  include <GL/osmesa.h>

#endif

#if defined(THEA_MAC) && !defined(THEA_GL_OSMESA)
#  include <OpenGL/OpenGL.h>  // required to pull in CGL, which OpenGL/gl.h does not
#endif

// GUI compatibility. Requires a display device.
#ifndef THEA_GL_OSMESA
#  ifdef THEA_WIN32
#    include "wglew.h"
#  elif defined(THEA_LINUX) || defined(THEA_BSD)
#    include "glxew.h"
#  endif
#endif

namespace Thea {
namespace Graphics {

/** OpenGL rendering backend. */
namespace GL {

/** RAII idiom for saving and restoring OpenGL server attribute states. */
struct THEA_GL_DLL_LOCAL GLScope
{
  /** Constructor. Creates a new scope for a set of attributes. */
  GLScope(GLbitfield attribs_to_save) { glPushAttrib(attribs_to_save); }

  /** Destructor. Restores saved attributes from the stack. */
  ~GLScope() { glPopAttrib(); }
};

/** RAII idiom for saving and restoring OpenGL client attribute states. */
struct THEA_GL_DLL_LOCAL GLClientScope
{
  /** Constructor. Creates a new scope for a set of attributes. */
  GLClientScope(GLbitfield attribs_to_save) { glPushClientAttrib(attribs_to_save); }

  /** Destructor. Restores saved attributes from the stack. */
  ~GLClientScope() { glPopClientAttrib(); }
};

/** Replacement for deprecated gluErrorString(). */
inline GLubyte const *
theaGLErrorString(GLuint error)
{
  switch (error)
  {
#define THEA_GL_ERROR_STRING_BRANCH(p) case(p): return (GLubyte const *)#p;
    THEA_GL_ERROR_STRING_BRANCH(GL_NO_ERROR)
    THEA_GL_ERROR_STRING_BRANCH(GL_INVALID_ENUM)
    THEA_GL_ERROR_STRING_BRANCH(GL_INVALID_VALUE)
    THEA_GL_ERROR_STRING_BRANCH(GL_INVALID_OPERATION)
    THEA_GL_ERROR_STRING_BRANCH(GL_STACK_OVERFLOW)
    THEA_GL_ERROR_STRING_BRANCH(GL_STACK_UNDERFLOW)
    THEA_GL_ERROR_STRING_BRANCH(GL_OUT_OF_MEMORY)
    THEA_GL_ERROR_STRING_BRANCH(GL_TABLE_TOO_LARGE)
    default: break;
#undef THEA_GL_ERROR_STRING_BRANCH
  }

  if ((error >= GLU_NURBS_ERROR1) && (error <= GLU_NURBS_ERROR37))
    return (GLubyte const *)"NURBS error";
  if ((error >= GLU_TESS_ERROR1) && (error <= GLU_TESS_ERROR8))
    return (GLubyte const *)"tessellation error";

  return (GLubyte const *)"unknown error";
}

/**
 * Check if the OpenGL state is error-free.
 *
 * @todo Make this a no-op in release mode?
 */
#define THEA_CHECK_GL_OK \
{ \
  GLenum err_code; \
  GLubyte const * err_string; \
\
  if ((err_code = glGetError()) != GL_NO_ERROR) \
  { \
    err_string = theaGLErrorString(err_code); \
    throw Error(Thea::format("%s:%ld: OpenGL error: %s", \
                Thea::FilePath::objectName(__FILE__).c_str(), (long)__LINE__, err_string)); \
  } \
}

/** A convenient shortcut for getting a single integer state. */
inline THEA_GL_DLL_LOCAL GLint
glGetInteger(GLenum which)
{
  // We allocate an array in case the caller accidentally invoked on a value that will return more than just one integer
  GLint result[32];
  glGetIntegerv(which, result);
  return result[0];
}

/**
 * Get a handle to the current OpenGL context, if available. On platforms where we can't return a context, this function always
 * returns true.
 */
#if defined(THEA_GL_OSMESA)

  typedef OSMesaContext GLContext;
  inline THEA_GL_DLL_LOCAL GLContext glGetCurrentContext()
  {
    return OSMesaGetCurrentContext();
  }

#elif defined(THEA_WINDOWS)

  typedef HGLRC GLContext;
  inline THEA_GL_DLL_LOCAL GLContext glGetCurrentContext()
  {
    return wglGetCurrentContext();
  }

#elif defined(THEA_LINUX) || defined(THEA_BSD)

  typedef GLXContext GLContext;
  inline THEA_GL_DLL_LOCAL GLContext glGetCurrentContext()
  {
    return glXGetCurrentContext();
  }

#elif defined(THEA_MAC)

  typedef CGLContextObj GLContext;
  inline THEA_GL_DLL_LOCAL GLContext glGetCurrentContext()
  {
      return CGLGetCurrentContext();
  }

#endif

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
