#===============================================================================================================================
#
# Build script for the Thea OpenGL rendering plugin.
#
# Copyright (C) 2010, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginGL CXX C)

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

# Dependency: OpenGL
IF(WITH_OSMESA)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-osmesa)
    SET(OSMesa_ROOT ${THEA_DEPS_ROOT}/installed-osmesa)
  ELSE()
    SET(OSMesa_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  FIND_PACKAGE(OSMesa REQUIRED)

  IF(OSMesa_FOUND AND OSMesa_GLU_FOUND)
    SET(OPENGL_FOUND TRUE)
    SET(OPENGL_GLU_FOUND TRUE)

    # Replace GL libraries with OSMesa ones, if available
    SET(OPENGL_INCLUDE_DIRS ${OSMesa_INCLUDE_DIRS} ${OSMesa_GLU_INCLUDE_DIRS})
    SET(OPENGL_CFLAGS "${OSMesa_CFLAGS} ${OSMesa_GLU_CFLAGS}")
    SET(OPENGL_LIBRARIES ${OSMesa_GLU_LIBRARIES} ${OSMesa_LIBRARIES})
    SET(OPENGL_LIBRARY_DIRS ${OSMesa_LIBRARY_DIRS} ${OSMesa_GLU_LIBRARY_DIRS})
    SET(OPENGL_LDFLAGS "${OSMesa_LDFLAGS} ${OSMesa_GLU_LDFLAGS}")
    ADD_DEFINITIONS(-DTHEA_GL_OSMESA=1 -DGLEW_OSMESA=1)
  ELSE()
    MESSAGE(FATAL_ERROR "OSMesa GL or GLU library not found")
  ENDIF()
ELSE()
  FIND_PACKAGE(OpenGL)
ENDIF()

# If we can't find OpenGL, don't build this plugin
IF(OPENGL_FOUND AND OPENGL_GLU_FOUND)

# Additional platform-specific libraries
IF(WIN32)
  SET(TheaPluginGL_PLATFORM_LIBRARIES "version.lib")
ENDIF()

ADD_DEFINITIONS(-DGLEW_STATIC)
IF(APPLE)
  ADD_DEFINITIONS(-DGL_SILENCE_DEPRECATION -DCL_SILENCE_DEPRECATION)
ENDIF()

# Shared library flags
ADD_DEFINITIONS(-DTHEA_GL_DLL -DTHEA_GL_DLL_EXPORTS)
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENGL_CFLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENGL_CFLAGS}")

# Include directories
INCLUDE_DIRECTORIES(BEFORE ${OPENGL_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${OPENGL_LIBRARY_DIRS})

# Source file lists
FILE(GLOB TheaPluginGLSources
     ${SourceRoot}/Plugins/GL/*.cpp
     ${SourceRoot}/Plugins/GL/glew.c)

# Libraries to link to
SET(TheaPluginGLLibraries
    Thea
    ${OPENGL_LIBRARIES}
    ${Thea_DEPS_LIBRARIES}
    ${TheaPluginGL_PLATFORM_LIBRARIES})

# Build products
ADD_LIBRARY(TheaPluginGL SHARED ${TheaPluginGLSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginGL ${TheaPluginGLLibraries})
SET_TARGET_PROPERTIES(TheaPluginGL PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS} ${OPENGL_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaPluginGL "${TheaPluginGLLibraries}")
ENDIF()

# Install rules
SET_TARGET_PROPERTIES(TheaPluginGL
                        PROPERTIES
                          INSTALL_RPATH_USE_LINK_PATH TRUE
                          INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

INSTALL(TARGETS TheaPluginGL DESTINATION lib)
INSTALL(TARGETS TheaPluginGL DESTINATION share/Thea/Build/Output/lib)  # for tools to use

ELSE(OPENGL_FOUND AND OPENGL_GLU_FOUND)

MESSAGE(STATUS "OPENGL not found: GL plugin will not be built")

ENDIF(OPENGL_FOUND AND OPENGL_GLU_FOUND)
