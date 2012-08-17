# Find Bullet includes and library
#
# This module defines
#  Bullet_INCLUDE_DIRS  The location of the Bullet header files
#  Bullet_LIBRARIES     The libraries to link against to use Bullet
#  Bullet_FOUND         If false, do not try to use Bullet
#
# To specify an additional directory to search, set Bullet_ROOT.
#
# If the absence of Bullet is a fatal error, set Bullet_FIND_REQUIRED to TRUE.
#
# Original copyright (c) 2007, Matt Williams
#
# Heavily modified to work on Unix, and with the Meru project conventions, by
# Siddhartha Chaudhuri, 2008.
#
# Redistribution and use is allowed according to the terms of the BSD license.

SET(Bullet_FOUND FALSE)

IF(WIN32) # Windows

  SET(Bullet_DIR ${Bullet_ROOT})
  IF(Bullet_DIR AND EXISTS "${Bullet_DIR}")
    SET(Bullet_INCLUDE_DIRS ${Bullet_DIR}/include)
    SET(Bullet_LIBRARIES
        debug ${Bullet_DIR}/lib/Debug/libbulletcollision_d.lib
        debug ${Bullet_DIR}/lib/Debug/libbulletdynamics_d.lib
        debug ${Bullet_DIR}/lib/Debug/libbulletmath_d.lib
        optimized ${Bullet_DIR}/lib/Release/libbulletcollision.lib
        optimized ${Bullet_DIR}/lib/Release/libbulletdynamics.lib
        optimized ${Bullet_DIR}/lib/Release/libbulletmath.lib)

    SET(Bullet_FOUND TRUE)
  ENDIF(Bullet_DIR AND EXISTS "${Bullet_DIR}")

ELSE(WIN32) # Unix

  SET(Bullet_DIR_DESCRIPTION "root directory of Bullet installation. E.g /usr/local or /opt")

  # Look for an installation, first in the user-specified location and then in the system locations
  FIND_PATH(Bullet_DIR NAMES lib/pkgconfig/bullet.pc PATHS ${Bullet_ROOT} DOC "The ${Bullet_DIR_DESCRIPTION}" NO_DEFAULT_PATH)
  IF(NOT Bullet_DIR)  # now look in system locations
    FIND_PATH(Bullet_DIR NAMES lib/pkgconfig/bullet.pc DOC "The ${Bullet_DIR_DESCRIPTION}")
  ENDIF(NOT Bullet_DIR)

  # Now try to get the include and library path.
  IF(Bullet_DIR)
    SET(Bullet_INCLUDE_DIRS ${Bullet_DIR}/include/bullet)

    IF(EXISTS ${Bullet_INCLUDE_DIRS})
      SET(Bullet_LIBRARY_DIR ${Bullet_DIR}/lib)
      FIND_LIBRARY(Bullet_LIB_DYNAMICS NAMES bulletdynamics PATHS ${Bullet_LIBRARY_DIR} NO_DEFAULT_PATH)
      FIND_LIBRARY(Bullet_LIB_COLLISION NAMES bulletcollision PATHS ${Bullet_LIBRARY_DIR} NO_DEFAULT_PATH)
      FIND_LIBRARY(Bullet_LIB_MATH NAMES bulletmath PATHS ${Bullet_LIBRARY_DIR} NO_DEFAULT_PATH)

      IF(Bullet_LIB_DYNAMICS AND Bullet_LIB_COLLISION AND Bullet_LIB_MATH)
        SET(Bullet_FOUND TRUE)
        SET(Bullet_LIBRARIES ${Bullet_LIB_DYNAMICS} ${Bullet_LIB_COLLISION} ${Bullet_LIB_MATH})
      ENDIF(Bullet_LIB_DYNAMICS AND Bullet_LIB_COLLISION AND Bullet_LIB_MATH)
    ENDIF(EXISTS ${Bullet_INCLUDE_DIRS})
  ENDIF(Bullet_DIR)

ENDIF(WIN32)

IF(Bullet_FOUND)
  IF(NOT Bullet_FIND_QUIETLY)
    MESSAGE(STATUS "Found Bullet at ${Bullet_DIR}")
  ENDIF(NOT Bullet_FIND_QUIETLY)
ELSE(Bullet_FOUND)
  IF(Bullet_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Bullet not found")
  ENDIF(Bullet_FIND_REQUIRED)
ENDIF(Bullet_FOUND)
