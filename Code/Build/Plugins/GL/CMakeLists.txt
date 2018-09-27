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
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

# See cmake --help-policy CMP0042 for details on this one
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 NEW)
ENDIF(POLICY CMP0042)

# See cmake --help-policy CMP0074 for details on this one
IF(POLICY CMP0074)
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF(POLICY CMP0074)

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
IF(NOT THEA_DEPS_ROOT)
  SET(THEA_DEPS_ROOT ${CMAKE_INSTALL_PREFIX})
ENDIF()
SET(THEA_DEPS_ROOT ${THEA_DEPS_ROOT} CACHE PATH "Path to folder with installations of dependencies")

# Locate dependencies
SET(Thea_FIND_Boost  TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Dependency: OSMesa
IF(WITH_OSMESA)
  FIND_PACKAGE(OSMesa REQUIRED)
  IF(OSMesa_FOUND AND OSMesa_GLU_FOUND)
    SET(OPENGL_FOUND TRUE)
    SET(OPENGL_GLU_FOUND TRUE)

    # Replace GL libraries with OSMesa ones, if available
    SET(OPENGL_INCLUDE_DIRS ${OSMesa_INCLUDE_DIRS} ${OSMesa_GLU_INCLUDE_DIRS})
    SET(OPENGL_CFLAGS ${OSMesa_CFLAGS} ${OSMesa_GLU_CFLAGS})
    SET(OPENGL_LIBRARIES ${OSMesa_GLU_LIBRARIES} ${OSMesa_LIBRARIES})
    SET(OPENGL_LIBRARY_DIRS ${OSMesa_LIBRARY_DIRS} ${OSMesa_GLU_LIBRARY_DIRS})
    SET(OPENGL_LDFLAGS ${OSMesa_LDFLAGS} ${OSMesa_GLU_LDFLAGS})
    ADD_DEFINITIONS(-DTHEA_GL_OSMESA=1)
  ELSE()
    MESSAGE(FATAL_ERROR "OSMesa GL or GLU library not found")
  ENDIF()
ELSE()
  FIND_PACKAGE(OpenGL)
ENDIF()

# If we can't find OpenGL, don't build this plugin
IF(OPENGL_FOUND AND OPENGL_GLU_FOUND)

# Additional platform-specific libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(PLATFORM_LIBRARIES "-framework Carbon")
ELSEIF(WIN32)
  SET(PLATFORM_LIBRARIES "version.lib")
ENDIF()

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  STRING(REPLACE ";" " " EXTRA_DEBUG_CFLAGS "${CGAL_DEBUG_CFLAGS}")
  STRING(REPLACE ";" " " EXTRA_RELEASE_CFLAGS "${CGAL_RELEASE_CFLAGS}")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENGL_CFLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENGL_CFLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

ELSEIF(MSVC)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
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
                    ${OPENGL_INCLUDE_DIRS})

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
INSTALL(TARGETS TheaPluginGL DESTINATION share/Thea/Build/Output/lib)  # for tools to use

ELSE(OPENGL_FOUND AND OPENGL_GLU_FOUND)

MESSAGE(STATUS "OPENGL not found: GL plugin will not be built")

ENDIF(OPENGL_FOUND AND OPENGL_GLU_FOUND)
