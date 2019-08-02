#===============================================================================================================================
#
# Build script for the Thea ARPACK eigensolver plugin.
#
# Copyright (C) 2010, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginARPACK)

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

# Dependency: ARPACK
IF(EXISTS ${THEA_DEPS_ROOT}/installed-arpack)
  SET(ARPACK_ROOT ${THEA_DEPS_ROOT}/installed-arpack)
ELSE()
  SET(ARPACK_ROOT ${THEA_DEPS_ROOT})
ENDIF()
FIND_PACKAGE(ARPACK)

# If we can't find ARPACK, don't build this plugin
IF(ARPACK_FOUND)

# Dependency: SuperLU (optional)
IF(EXISTS ${THEA_DEPS_ROOT}/installed-superlu)
  SET(SuperLU_ROOT ${THEA_DEPS_ROOT}/installed-superlu)
ELSE()
  SET(SuperLU_ROOT ${THEA_DEPS_ROOT})
ENDIF()
SET(SuperLU_LIBS_ONLY TRUE)
FIND_PACKAGE(SuperLU)

IF(SuperLU_FOUND)
  ADD_DEFINITIONS(-DTHEA_HAVE_SUPERLU)
ELSE()
  SET(SuperLU_INCLUDE_DIRS "")
  SET(SuperLU_LIBRARIES "")
ENDIF()

# Shared library flags
ADD_DEFINITIONS(-DTHEA_ARPACK_DLL -DTHEA_ARPACK_DLL_EXPORTS)
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# ARPACK++ throws a bunch of harmless (?) warnings which are easier to ignore than to fix
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overloaded-virtual")
ENDIF()

# Source file lists
FILE(GLOB TheaPluginARPACKSources ${SourceRoot}/Plugins/ARPACK/*.cpp)

# Libraries to link to
SET(TheaPluginARPACKLibraries
    Thea
    ${SuperLU_LIBRARIES}
    ${ARPACK_LIBRARIES}
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_LIBRARY(TheaPluginARPACK SHARED ${TheaPluginARPACKSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginARPACK ${TheaPluginARPACKLibraries})
SET_TARGET_PROPERTIES(TheaPluginARPACK PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
INSTALL(TARGETS TheaPluginARPACK DESTINATION share/Thea/Build/Output/lib)  # for tools to use

ELSE(ARPACK_FOUND)

MESSAGE(STATUS "ARPACK not found: ARPACK plugin will not be built")

ENDIF(ARPACK_FOUND)
