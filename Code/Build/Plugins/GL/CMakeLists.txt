#===============================================================================================================================
#
# Build script for the Thea OpenGL rendering plugin.
#
# Copyright (C) 2010, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginGL)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# See cmake --help-policy CMP0003 for details on this one
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# If you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Postfix for debug builds
SET(CMAKE_DEBUG_POSTFIX "d")

# Project root path
GET_FILENAME_COMPONENT(ProjectRoot ../../.. ABSOLUTE)

# Path for build products
SET(OutputRoot ${ProjectRoot}/Build/Output)

# Path to put executables in
SET(EXECUTABLE_OUTPUT_PATH ${OutputRoot}/bin)

# Path to put libraries in
SET(LIBRARY_OUTPUT_PATH ${OutputRoot}/lib)

# Path for customized CMake modules
IF(NOT CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH ${ProjectRoot}/Build/Common/CMake/Modules)
ENDIF()
GET_FILENAME_COMPONENT(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ABSOLUTE)

# Path to root folder for source code
SET(SourceRoot ${ProjectRoot}/Source)

# Path to folder with installations of the dependencies
IF(NOT THEA_INSTALLATIONS_ROOT)
  SET(THEA_INSTALLATIONS_ROOT ${CMAKE_INSTALL_PREFIX})
ENDIF()
SET(THEA_INSTALLATIONS_ROOT ${THEA_INSTALLATIONS_ROOT} CACHE PATH "Path to folder with installations of dependencies")

# Locate dependencies
SET(Thea_FIND_Boost  TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# -- TODO: Spin off into separate offscreen rendering plugin --
# # Dependency: OSMesa
# IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-osmesa)
#   SET(OSMesa_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-osmesa)
# ELSE()
#   SET(OSMesa_ROOT ${THEA_INSTALLATIONS_ROOT})
# ENDIF()
# FIND_PACKAGE(OSMesa REQUIRED)
#
# IF(NOT OSMesa_GLU_FOUND)
#   MESSAGE(FATAL_ERROR "OSMesa-compatible GLU was not found")
# ENDIF()
# SET(OpenGL_INCLUDE_DIRS ${OSMesa_INCLUDE_DIRS} ${OSMesa_GLU_INCLUDE_DIRS})
# SET(OpenGL_CFLAGS ${OSMesa_CFLAGS} ${OSMesa_GLU_CFLAGS})
# SET(OpenGL_LIBRARIES ${OSMesa_LIBRARIES} ${OSMesa_GLU_LIBRARIES})
# SET(OpenGL_LIBRARY_DIRS ${OSMesa_LIBRARY_DIRS} ${OSMesa_GLU_LIBRARY_DIRS})
# SET(OpenGL_LDFLAGS ${OSMesa_LDFLAGS} ${OSMesa_GLU_LDFLAGS})
#
# ADD_DEFINITIONS(-DTHEA_GL_OSMESA)
# ADD_DEFINITIONS(-DTHEA_DISABLE_DISPLAY)

# Dependency: OpenGL
FIND_PACKAGE(OpenGL REQUIRED)
IF(NOT OPENGL_GLU_FOUND)
  MESSAGE(FATAL_ERROR "GLU was not found")
ENDIF()
SET(OpenGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
SET(OpenGL_CFLAGS )
SET(OpenGL_LIBRARIES    ${OPENGL_LIBRARIES})
SET(OpenGL_LIBRARY_DIRS )
SET(OpenGL_LDFLAGS )
IF(APPLE)
  ADD_DEFINITIONS(-DTHEA_GL_USE_FRAMEWORK)
ENDIF()
MESSAGE(STATUS "Found OpenGL: headers at ${OpenGL_INCLUDE_DIRS}, libraries at ${OpenGL_LIBRARIES}")

# Additional platform-specific libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(PLATFORM_LIBRARIES "-framework Carbon")
ENDIF()

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(${CGAL_DEBUG_CFLAGS} ${OpenGL_CFLAGS} -Wall -g2 -fno-strict-aliasing)  # -fvisibility=hidden
  ELSE()
    ADD_DEFINITIONS(${CGAL_RELEASE_CFLAGS} ${OpenGL_CFLAGS} -Wall -O2 -fno-strict-aliasing -DNDEBUG)  # -fvisibility=hidden
  ENDIF()
ENDIF()
ADD_DEFINITIONS(-DGLEW_STATIC)

# Shared library flags
ADD_DEFINITIONS(-DTHEA_GL_DLL -DTHEA_GL_DLL_EXPORTS)
ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)

# "extern template" support
IF(NOT DEFINED THEA_EXTERN_TEMPLATES)
  SET(THEA_EXTERN_TEMPLATES FALSE)
ENDIF()
SET(THEA_EXTERN_TEMPLATES ${THEA_EXTERN_TEMPLATES} CACHE BOOL "Use extern templates?")

IF(THEA_EXTERN_TEMPLATES)
  MESSAGE(STATUS "Compiler support for 'extern template' required")
  ADD_DEFINITIONS(-DTHEA_EXTERN_TEMPLATES)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE
                    ${Boost_INCLUDE_DIRS}
                    ${OpenGL_INCLUDE_DIRS})

# Source file lists
FILE(GLOB TheaPluginGLSources
     ${SourceRoot}/Plugins/GL/*.cpp
     ${SourceRoot}/Plugins/GL/glew.c)

# Libraries to link to
SET(TheaPluginGLLibraries
    Thea
    ${OpenGL_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_LIBRARY(TheaPluginGL SHARED ${TheaPluginGLSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginGL ${TheaPluginGLLibraries})

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
