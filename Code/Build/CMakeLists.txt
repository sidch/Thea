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

# https://stackoverflow.com/a/42697475
MACRO(install_symlink filepath sympath)
  INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})")
  INSTALL(CODE "MESSAGE(STATUS \"Created symlink: ${sympath} -> ${filepath}\")")
ENDMACRO(install_symlink)

# Subdirectories containing individual build targets
ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(Plugins/ARPACK)
ADD_SUBDIRECTORY(Plugins/CSPARSE)
ADD_SUBDIRECTORY(Plugins/GL)

# Omit OPT++ from the default build until we finish this module
IF(WITH_OPTPP)
  ADD_SUBDIRECTORY(Plugins/OPTPP)
ENDIF()

IF(NOT DEFINED WITH_TESTS)
  SET(WITH_TESTS true)
ENDIF()
IF(WITH_TESTS)
  MESSAGE(STATUS "")
  MESSAGE(STATUS "NOTE: Test programs run when you execute 'make test'.")
  MESSAGE(STATUS "NOTE: Some tests may be out-of-date with the current library.")
  MESSAGE(STATUS "")
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(Test)
ENDIF()

ADD_SUBDIRECTORY(Tools)
