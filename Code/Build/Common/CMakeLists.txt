#===============================================================================================================================
#
# Build script for the Thea library.
#
# Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(Thea)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# See cmake --help-policy CMP0003 for details on this one
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

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
IF(NOT THEA_DEPS_ROOT)
  SET(THEA_DEPS_ROOT ${CMAKE_INSTALL_PREFIX})
ENDIF()
SET(THEA_DEPS_ROOT ${THEA_DEPS_ROOT} CACHE PATH "Path to folder with installations of dependencies")

# Locate dependencies
SET(Thea_FIND_ALL TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  STRING(REPLACE ";" " " EXTRA_DEBUG_CFLAGS "${CGAL_DEBUG_CFLAGS}")
  STRING(REPLACE ";" " " EXTRA_RELEASE_CFLAGS "${CGAL_RELEASE_CFLAGS}")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-strict-aliasing -fPIC")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fno-strict-aliasing -fPIC")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

ELSEIF(MSVC)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ENDIF()

# Additional platform-specific libraries
SET(PLATFORM_LIBRARIES ${CMAKE_DL_LIBS})

# Specify that Triangle should be built as a library, use ANSI function declarations and not use a timer
ADD_DEFINITIONS(-DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER)

# Shared library flags
IF(THEA_DLL)
  MESSAGE(WARNING "
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                        !
  !  WARNING!   WARNING!   WARNING!   WARNING!   WARNING!  !
  !                                                        !
  !  You're trying to build Thea as a shared library!      !
  !                                                        !
  !  This can be VERY DANGEROUS since things can break     !
  !  mysteriously if you mix compilers/build settings (or  !
  !  possibly even with the same compiler and settings).   !
  !                                                        !
  !  Thea uses STL classes and lots of templates in its    !
  !  public interface. This is SUPER-BAD for shared        !
  !  libraries.                                            !
  !                                                        !
  !  Proceed at your own risk. This is NOT recommended.    !
  !                                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ")

  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_EXPORTS)
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
                    ${FreeImage_INCLUDE_DIRS}
                    ${CGAL_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS})
IF(CLUTO_FOUND)
  INCLUDE_DIRECTORIES(BEFORE ${CLUTO_INCLUDE_DIRS})
ENDIF()

# Source file lists
FILE(GLOB TheaSources
     ${SourceRoot}/*.cpp
     ${SourceRoot}/Triangle/*.c
     ${SourceRoot}/Algorithms/*.cpp
     ${SourceRoot}/Algorithms/ImageFeatures/*.cpp
     ${SourceRoot}/Algorithms/MeshFeatures/*.cpp
     ${SourceRoot}/Algorithms/BloomenthalPolygonizer/*.cpp
     ${SourceRoot}/Algorithms/dx/*.c
     ${SourceRoot}/Algorithms/fibheap/*.cpp
     ${SourceRoot}/Algorithms/miniball/*.cpp
     ${SourceRoot}/Algorithms/NNLS/*.c
     ${SourceRoot}/Algorithms/StanHull/*.cpp
     ${SourceRoot}/Graphics/*.cpp)

# Libraries to link to
SET(TheaLibraries
    ${Boost_LIBRARIES}
    ${CGAL_LIBRARIES}
    ${Lib3ds_LIBRARIES}
    ${FreeImage_LIBRARIES}
    ${PLATFORM_LIBRARIES})
IF(CLUTO_FOUND)
  SET(TheaLibraries ${CLUTO_LIBRARIES} ${TheaLibraries})
ENDIF()


# Build products
ADD_LIBRARY(Thea STATIC ${TheaSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(Thea ${TheaLibraries})

# Fix library install names on OS X
IF(THEA_DLL AND APPLE)
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
