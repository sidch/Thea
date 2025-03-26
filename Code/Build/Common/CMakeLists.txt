#===============================================================================================================================
#
# Build script for the Thea library.
#
# Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

PROJECT(TheaLib CXX C)

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

# Source file lists
FILE(GLOB TheaSources
     ${ThirdPartyRoot}/BloomenthalPolygonizer/*.cpp
     ${ThirdPartyRoot}/dx/*.c
     ${ThirdPartyRoot}/fibheap/*.cpp
     ${ThirdPartyRoot}/miniball/*.cpp
     ${ThirdPartyRoot}/NNLS/*.c
     ${ThirdPartyRoot}/StanHull/*.cpp
     ${ThirdPartyRoot}/Triangle/*.c
     ${SourceRoot}/*.cpp
     ${SourceRoot}/Algorithms/*.cpp
     ${SourceRoot}/Algorithms/ImageFeatures/*.cpp
     ${SourceRoot}/Algorithms/SurfaceFeatures/Global/*.cpp
     ${SourceRoot}/Algorithms/SurfaceFeatures/Local/*.cpp
     ${SourceRoot}/Graphics/*.cpp)

# Libraries to link to
SET(TheaLibraries ${Thea_DEPS_LIBRARIES})

# Build products
ADD_LIBRARY(Thea STATIC ${TheaSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(Thea ${TheaLibraries})
SET_TARGET_PROPERTIES(Thea PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(THEA_DLL AND APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(Thea "${TheaLibraries}")
ENDIF()

# Write a file containing build flags to be used by a calling application
SET(Thea_BUILD_FLAGS
    "EXTERN_TEMPLATES=${EXTERN_TEMPLATES}"
    "LITE=${LITE}"
    "WITH_CGAL=${WITH_CGAL}"
    "WITH_CLUTO=${WITH_CLUTO}"
    "WITH_FREEIMAGE=${WITH_FREEIMAGE}"
    "WITH_LIB3DS=${WITH_LIB3DS}"
    "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "CMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}")
STRING(REPLACE ";" "\n" Thea_BUILD_FLAGS "${Thea_BUILD_FLAGS}")
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
          PATTERN ".git" EXCLUDE
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
          PATTERN ".git" EXCLUDE
          PATTERN "*.cpp" EXCLUDE
          PATTERN "*.c" EXCLUDE)
ENDIF()
