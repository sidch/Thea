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

# Mac-specific requirements: don't link to Lion SDK
IF(APPLE)
  # Determine the version of OS X that is installed
  EXECUTE_PROCESS(
    COMMAND sw_vers -productVersion
    ERROR_QUIET
    OUTPUT_VARIABLE SW_CMD
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  STRING(REGEX REPLACE "^(10.[0-9]+).*" "\\1" SW_VER "${SW_CMD}")
  IF(SW_VER STREQUAL "10.7" OR SW_VER STREQUAL "10.8")
    SET(CMAKE_OSX_SYSROOT /Developer/SDKs/MacOSX10.6.sdk CACHE STRING "Path to OS X SDK root" FORCE)
    MESSAGE(STATUS "Setting SDK to 10.6: ${CMAKE_OSX_SYSROOT}")
  ENDIF()
ENDIF()

# Subdirectories containing individual build targets
ADD_SUBDIRECTORY(Common)

IF(NOT APPLE AND NOT WIN32)  # we don't have arpack/superlu working on OS X or Windows yet
  ADD_SUBDIRECTORY(Plugins/ARPACK)
ENDIF()

ADD_SUBDIRECTORY(Plugins/CSPARSE)
ADD_SUBDIRECTORY(Plugins/GL)
ADD_SUBDIRECTORY(Plugins/OPTPP)

# The tests should not be built by default
IF(WITH_TESTS)
  ADD_SUBDIRECTORY(Test EXCLUDE_FROM_ALL)

  SET(WITH_TESTS ${WITH_TESTS} CACHE BOOL "Build test programs?")

  # Help message about building tests
  MESSAGE(STATUS "")
  MESSAGE(STATUS "NOTE: Test programs are not automatically built. To build them, run 'make TheaTests'.")
ELSE()
  MESSAGE(STATUS "")
  MESSAGE(STATUS "NOTE: To build test programs, set the variable WITH_TESTS to true")
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "NOTE: Tools are not automatically built. To build them, run 'make TheaTools'.")
MESSAGE(STATUS "")
ADD_SUBDIRECTORY(Tools EXCLUDE_FROM_ALL)
