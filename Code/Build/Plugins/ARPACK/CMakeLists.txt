#===============================================================================================================================
#
# Build script for the Thea ARPACK eigensolver plugin.
#
# Copyright (C) 2010, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginARPACK)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# See cmake --help-policy CMP0003 for details on this one
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

# See cmake --help-policy CMP0042 for details on this one
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 OLD)
ENDIF(POLICY CMP0042)

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

# Dependency: ARPACK
IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-arpack)
  SET(ARPACK_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-arpack)
ELSE()
  SET(ARPACK_ROOT ${THEA_INSTALLATIONS_ROOT})
ENDIF()
FIND_PACKAGE(ARPACK)

# Dependency: ARPACK++
IF(ARPACK_FOUND)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-arpack++)
    SET(ARPACKPP_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-arpack++)
  ELSE()
    SET(ARPACKPP_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  FIND_PACKAGE(ARPACKPP)
ENDIF(ARPACK_FOUND)

# If we can't find ARPACK or ARPACK++, don't build this plugin
IF(ARPACK_FOUND AND ARPACKPP_FOUND)

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
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  STRING(REPLACE ";" " " EXTRA_DEBUG_CFLAGS "${CGAL_DEBUG_CFLAGS}")
  STRING(REPLACE ";" " " EXTRA_RELEASE_CFLAGS "${CGAL_RELEASE_CFLAGS}")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

ELSEIF(MSVC)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ENDIF()

# Shared library flags
ADD_DEFINITIONS(-DTHEA_ARPACK_DLL -DTHEA_ARPACK_DLL_EXPORTS)
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
                    ${ARPACKPP_INCLUDE_DIRS})

# Source file lists
FILE(GLOB TheaPluginARPACKSources ${SourceRoot}/Plugins/ARPACK/*.cpp)

# Libraries to link to
SET(TheaPluginARPACKLibraries
    Thea
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

ELSE(ARPACK_FOUND AND ARPACKPP_FOUND)

MESSAGE(STATUS "ARPACK/ARPACK++ not found: ARPACK plugin will not be built")

ENDIF(ARPACK_FOUND AND ARPACKPP_FOUND)
