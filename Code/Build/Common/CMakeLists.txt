#===============================================================================================================================
#
# Build script for the Thea library.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#
#===============================================================================================================================

PROJECT(Thea)

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
GET_FILENAME_COMPONENT(ProjectRoot ../.. ABSOLUTE)

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
SET(Thea_FIND_ALL TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(${CGAL_DEBUG_CFLAGS} -Wall -g2 -fno-strict-aliasing)  # -fvisibility=hidden
  ELSE()
    ADD_DEFINITIONS(${CGAL_RELEASE_CFLAGS} -Wall -O2 -fno-strict-aliasing -DNDEBUG)  # -fvisibility=hidden
  ENDIF()
ELSEIF(MSVC)
  IF(NOT CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(/O2)
  ENDIF()
ENDIF()

# Specify that Triangle should be built as a library, use ANSI function declarations and not use a timer
ADD_DEFINITIONS(-DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER)

# Shared library flags
ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_EXPORTS)

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
                    ${FreeImage_INCLUDE_DIRS}
                    ${G3D_INCLUDE_DIRS}
                    ${CGAL_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS}
                    ${CLUTO_INCLUDE_DIRS})

# Source file lists
FILE(GLOB TheaSources
     ${SourceRoot}/*.cpp
     ${SourceRoot}/Algorithms/*.cpp
     ${SourceRoot}/Algorithms/dx/*.c
     ${SourceRoot}/Algorithms/fibheap/*.cpp
     ${SourceRoot}/Algorithms/NNLS/*.c
     ${SourceRoot}/Algorithms/BloomenthalPolygonizer/*.cpp
     ${SourceRoot}/Algorithms/StanHull/*.cpp
     ${SourceRoot}/Algorithms/Triangle/*.c
     ${SourceRoot}/Database/*.cpp
     ${SourceRoot}/Graphics/*.cpp)

# Libraries to link to
SET(TheaLibraries
    ${CLUTO_LIBRARIES}
    ${CGAL_LIBRARIES}
    ${Lib3ds_LIBRARIES}
    ${G3D_LIBRARIES}
    ${FreeImage_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_LIBRARY(Thea SHARED ${TheaSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(Thea ${TheaLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(Thea "${TheaLibraries}")
ENDIF()

# Write a file containing build flags to be used by a calling application
SET(Thea_BUILD_FLAGS )
IF(THEA_EXTERN_TEMPLATES)
  SET(Thea_BUILD_FLAGS "${Thea_BUILD_FLAGS} -DTHEA_EXTERN_TEMPLATES")
ENDIF(THEA_EXTERN_TEMPLATES)
FILE(WRITE "${SourceRoot}/BuildFlags.txt" "${Thea_BUILD_FLAGS}\n")

# Install rules
SET_TARGET_PROPERTIES(Thea PROPERTIES
                           INSTALL_RPATH_USE_LINK_PATH TRUE
                           INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

SET(UseFilesMatching FALSE)
IF(${CMAKE_MAJOR_VERSION} GREATER 2)
  SET(UseFilesMatching TRUE)
ELSEIF(${CMAKE_MAJOR_VERSION} EQUAL 2)
  IF(${CMAKE_MINOR_VERSION} GREATER 4)
    SET(UseFilesMatching TRUE)
  ENDIF()
ENDIF()

INSTALL(TARGETS Thea DESTINATION lib)
IF(UseFilesMatching)  # CMake 2.6 and above
  INSTALL(DIRECTORY ${SourceRoot}/ DESTINATION include/Thea
          FILES_MATCHING
          REGEX "Test$" EXCLUDE
          REGEX "Plugins$" EXCLUDE
          REGEX "Tools$" EXCLUDE
          PATTERN "CVS" EXCLUDE
          PATTERN ".svn" EXCLUDE
          PATTERN "*.hpp"
          PATTERN "*.h"
          PATTERN "BuildFlags.txt")
ELSE()
  # Bit flakier method, since we can't install only the files matching a pattern
  INSTALL(DIRECTORY ${SourceRoot}/ DESTINATION include/Thea
          REGEX "Test$" EXCLUDE
          REGEX "Plugins$" EXCLUDE
          REGEX "Tools$" EXCLUDE
          PATTERN "CVS" EXCLUDE
          PATTERN ".svn" EXCLUDE
          PATTERN "*.cpp" EXCLUDE
          PATTERN "*.c" EXCLUDE)
ENDIF()
