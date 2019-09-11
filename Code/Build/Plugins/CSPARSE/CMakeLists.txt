#===============================================================================================================================
#
# Build script for the Thea CSPARSE linear solver plugin.
#
# Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaPluginCSPARSE CXX C)

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

# Shared library flags
ADD_DEFINITIONS(-DTHEA_CSPARSE_DLL -DTHEA_CSPARSE_DLL_EXPORTS)
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Source file lists
FILE(GLOB TheaPluginCSPARSESources
     ${SourceRoot}/Plugins/CSPARSE/*.cpp
     ${SourceRoot}/Plugins/CSPARSE/*.c)

# Libraries to link to
SET(TheaPluginCSPARSELibraries
    Thea
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_LIBRARY(TheaPluginCSPARSE SHARED ${TheaPluginCSPARSESources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaPluginCSPARSE ${TheaPluginCSPARSELibraries})
SET_TARGET_PROPERTIES(TheaPluginCSPARSE PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaPluginCSPARSE "${TheaPluginCSPARSELibraries}")
ENDIF()

# Install rules
SET_TARGET_PROPERTIES(TheaPluginCSPARSE
                        PROPERTIES
                          INSTALL_RPATH_USE_LINK_PATH TRUE
                          INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

INSTALL(TARGETS TheaPluginCSPARSE DESTINATION lib)
INSTALL(TARGETS TheaPluginCSPARSE DESTINATION share/Thea/Build/Output/lib)  # for tools to use
