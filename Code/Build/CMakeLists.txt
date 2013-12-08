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
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# If you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Postfix for debug builds
SET(CMAKE_DEBUG_POSTFIX "d")

# Set the default build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

# Subdirectories containing individual build targets
ADD_SUBDIRECTORY(Common)

IF(NOT APPLE AND NOT WIN32)  # we don't have arpack/superlu working on OS X or Windows yet
#   ADD_SUBDIRECTORY(Plugins/ARPACK)
ENDIF()

ADD_SUBDIRECTORY(Plugins/CSPARSE)
ADD_SUBDIRECTORY(Plugins/GL)

IF(NOT WIN32)  # we don't have a prebuilt version of OPT++
  ADD_SUBDIRECTORY(Plugins/OPTPP)
ENDIF()

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
