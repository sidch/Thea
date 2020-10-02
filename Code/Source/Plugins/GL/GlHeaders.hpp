//============================================================================
//
// From the G3D project. Copyright (C) 2000-2008 Morgan McGuire.
//
// For the full G3D license see LICENSE.txt in the documentation.
//
//============================================================================

#ifndef __Thea_Graphics_Gl_GlHeaders_hpp__
#define __Thea_Graphics_Gl_GlHeaders_hpp__

#include "GlCommon.hpp"
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
namespace Gl {

/** RAII idiom for saving and restoring OpenGL server attribute states. */
struct THEA_GL_DLL_LOCAL GlScope
{
  /** Constructor. Creates a new scope for a set of attributes. */
  GlScope(GLbitfield attribs_to_save) { glPushAttrib(attribs_to_save); }

  /** Destructor. Restores saved attributes from the stack. */
  ~GlScope() { glPopAttrib(); }
};

/** RAII idiom for saving and restoring OpenGL client attribute states. */
struct THEA_GL_DLL_LOCAL GlClientScope
{
  /** Constructor. Creates a new scope for a set of attributes. */
  GlClientScope(GLbitfield attribs_to_save) { glPushClientAttrib(attribs_to_save); }

  /** Destructor. Restores saved attributes from the stack. */
  ~GlClientScope() { glPopClientAttrib(); }
};

/** Replacement for deprecated gluErrorString(). */
inline char const *
theaGlErrorString(GLuint error)
{
  switch (error)
  {
#define THEA_GL_ERROR_STRING_BRANCH(p) case(p): return #p;
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
    return "NURBS error";
  if ((error >= GLU_TESS_ERROR1) && (error <= GLU_TESS_ERROR8))
    return "tessellation error";

  return "unknown error";
}

/**
 * Check if the OpenGL state is error-free and throw an error if not.
 *
 * @see THEA_GL_OK_OR_RETURN(0)
 */
#define THEA_CHECK_GL_OK \
{ \
  GLenum err_code; \
  char const * err_string; \
\
  if ((err_code = glGetError()) != GL_NO_ERROR) \
  { \
    err_string = Thea::Graphics::Gl::theaGlErrorString(err_code); \
    throw Error(Thea::format("%s:%ld: OpenGL error: %s", \
                Thea::FilePath::objectName(__FILE__).c_str(), (long)__LINE__, err_string)); \
  } \
}

/**
 * Check if the OpenGL state is error-free and return zero/false if not.
 *
 * @see THEA_CHECK_GL_OK
 */
#define THEA_GL_OK_OR_RETURN(retval) \
{ \
  GLenum err_code; \
  if ((err_code = glGetError()) != GL_NO_ERROR) \
  { \
    THEA_ERROR << Thea::FilePath::objectName(__FILE__) << ':' << __LINE__ << ": OpenGL error: " \
               << Thea::Graphics::Gl::theaGlErrorString(err_code); \
    return (retval); \
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

  typedef OSMesaContext GlContext;
  inline THEA_GL_DLL_LOCAL GlContext glGetCurrentContext()
  {
    return OSMesaGetCurrentContext();
  }

#elif defined(THEA_WINDOWS)

  typedef HGLRC GlContext;
  inline THEA_GL_DLL_LOCAL GlContext glGetCurrentContext()
  {
    return wglGetCurrentContext();
  }

#elif defined(THEA_LINUX) || defined(THEA_BSD)

  typedef GLXContext GlContext;
  inline THEA_GL_DLL_LOCAL GlContext glGetCurrentContext()
  {
    return glXGetCurrentContext();
  }

#elif defined(THEA_MAC)

  typedef CGLContextObj GlContext;
  inline THEA_GL_DLL_LOCAL GlContext glGetCurrentContext()
  {
    return CGLGetCurrentContext();
  }

#endif

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
