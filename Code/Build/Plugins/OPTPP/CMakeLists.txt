#===============================================================================================================================
#
# Build script for the Thea OPT++ nonlinear optimization plugin.
#
# Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginOPTPP CXX C)

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

# Dependency: OPT++
IF(EXISTS ${THEA_DEPS_ROOT}/installed-opt++)
  SET(OPTPP_ROOT ${THEA_DEPS_ROOT}/installed-opt++)
ELSE()
  SET(OPTPP_ROOT ${THEA_DEPS_ROOT})
ENDIF()
FIND_PACKAGE(OPTPP REQUIRED)

# Assume our compiler is up to snuff
ADD_DEFINITIONS("-DHAVE_STD" "-DHAVE_NAMESPACES")

# Shared library flags
ADD_DEFINITIONS(-DTHEA_OPTPP_DLL -DTHEA_OPTPP_DLL_EXPORTS)
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE ${OPTPP_INCLUDE_DIRS})

# Libraries to link to
SET(TheaPluginOPTPPLibraries
    Thea
    ${OPTPP_LIBRARIES}
    ${Thea_DEPS_LIBRARIES})

# Source file lists
FILE(GLOB TheaPluginOPTPPSources
     ${SourceRoot}/Plugins/OPTPP/*.cpp)

# Build products
ADD_LIBRARY(TheaPluginOPTPP SHARED ${TheaPluginOPTPPSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginOPTPP ${TheaPluginOPTPPLibraries})
SET_TARGET_PROPERTIES(TheaPluginOPTPP PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS} ${OPTPP_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaPluginOPTPP "${TheaPluginOPTPPLibraries}")
ENDIF()

# Install rules
SET_TARGET_PROPERTIES(TheaPluginOPTPP
                        PROPERTIES
                          INSTALL_RPATH_USE_LINK_PATH TRUE
                          INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

INSTALL(TARGETS TheaPluginOPTPP DESTINATION lib)
INSTALL(TARGETS TheaPluginOPTPP DESTINATION share/Thea/Build/Output/lib)  # for tools to use
