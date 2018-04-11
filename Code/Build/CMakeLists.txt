#===============================================================================================================================
#
# Build script for the Thea library.
#
# Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(Thea)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

# See cmake --help-policy CMP0003 for details on this one
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

# If you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Postfix for debug builds
SET(CMAKE_DEBUG_POSTFIX "d")

# Set the default build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
ENDIF()

# Subdirectories containing individual build targets
ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(Plugins/ARPACK)
ADD_SUBDIRECTORY(Plugins/CSPARSE)
ADD_SUBDIRECTORY(Plugins/GL)

# Omit OPT++ from the default build until we finish this module
IF(WITH_OPTPP)
  ADD_SUBDIRECTORY(Plugins/OPTPP)
ENDIF()

# The tests should not be built by default (some of them are also probably currently broken)
IF(WITH_TESTS)
  ADD_SUBDIRECTORY(Test EXCLUDE_FROM_ALL)

  SET(WITH_TESTS ${WITH_TESTS} CACHE BOOL "Build test programs?")

  # Help message about building tests
  MESSAGE(STATUS "")
  MESSAGE(STATUS "NOTE: Test programs are not automatically built. To build them, run 'make test'.")
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "NOTE: Tools are not automatically built. To build them, run 'make tools'.")
MESSAGE(STATUS "")
ADD_SUBDIRECTORY(Tools EXCLUDE_FROM_ALL)
