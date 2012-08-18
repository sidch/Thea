#===============================================================================================================================
#
# Build script for the Thea ARPACK eigensolver plugin.
#
# Copyright (C) Siddhartha Chaudhuri, 2010
#
#===============================================================================================================================

PROJECT(TheaPluginARPACK)

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
SET(Thea_FIND_G3D    TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Dependency: ARPACK
IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-arpack)
  SET(ARPACK_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-arpack)
ELSE()
  SET(ARPACK_ROOT ${THEA_INSTALLATIONS_ROOT})
ENDIF()
FIND_PACKAGE(ARPACK REQUIRED)

# Dependency: ARPACK++
IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-arpack++)
  SET(ARPACKPP_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-arpack++)
ELSE()
  SET(ARPACKPP_ROOT ${THEA_INSTALLATIONS_ROOT})
ENDIF()
FIND_PACKAGE(ARPACKPP REQUIRED)

# Dependency: SuperLU (optional)
IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-superlu)
  SET(SuperLU_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-superlu)
ELSE()
  SET(SuperLU_ROOT ${THEA_INSTALLATIONS_ROOT})
ENDIF()
SET(SuperLU_LIBS_ONLY TRUE)
FIND_PACKAGE(SuperLU)

IF(SuperLU_FOUND)
  ADD_DEFINITIONS(-DTHEA_HAVE_SUPERLU)
ELSE()
  SET(SuperLU_INCLUDE_DIRS "")
  SET(SuperLU_LIBRARIES "")
ENDIF()

# Additional platform-specific libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(PLATFORM_LIBRARIES "-framework Carbon")
ENDIF()

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(${CGAL_DEBUG_CFLAGS} -Wall -g2 -fno-strict-aliasing)  # -fvisibility=hidden
  ELSE()
    ADD_DEFINITIONS(${CGAL_RELEASE_CFLAGS} -Wall -O2 -fno-strict-aliasing -DNDEBUG)  # -fvisibility=hidden
  ENDIF()
ENDIF()

# Shared library flags
ADD_DEFINITIONS(-DTHEA_ARPACK_DLL -DTHEA_ARPACK_DLL_EXPORTS)
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
                    ${G3D_INCLUDE_DIRS}
                    ${ARPACKPP_INCLUDE_DIRS})

# Source file lists
FILE(GLOB TheaPluginARPACKSources ${SourceRoot}/Plugins/ARPACK/*.cpp)

# Libraries to link to
SET(TheaPluginARPACKLibraries
    Thea
    ${G3D_LIBRARIES}
    ${SuperLU_LIBRARIES}
    ${ARPACK_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_LIBRARY(TheaPluginARPACK SHARED ${TheaPluginARPACKSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginARPACK ${TheaPluginARPACKLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaPluginARPACK "${TheaPluginARPACKLibraries}")
ENDIF()

# Install rules
SET_TARGET_PROPERTIES(TheaPluginARPACK
                        PROPERTIES
                          INSTALL_RPATH_USE_LINK_PATH TRUE
                          INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

INSTALL(TARGETS TheaPluginARPACK DESTINATION lib)
