# Searches for an installation of the SLEPc library. On success, it sets the following variables:
#
#   SLEPc_FOUND              Set to true to indicate the library was found
#   SLEPc_INCLUDE_DIRS       The directory containing the main SLEPc header file slepc.h
#   SLEPc_LIBRARIES          The libraries needed to use SLEPc (with full path)
#
# To specify an additional directory to search, set SLEPc_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2008
#

# Look for the SLEPc header, first in the user-specified location and then in the system locations
SET(SLEPc_INCLUDE_DOC "The directory containing the SLEPc header file slepc.h")
FIND_PATH(SLEPc_INCLUDE_DIRS NAMES slepc.h PATHS ${SLEPc_ROOT} ${SLEPc_ROOT}/include DOC ${SLEPc_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT SLEPc_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(SLEPc_INCLUDE_DIRS NAMES slepc.h DOC ${SLEPc_INCLUDE_DOC})
ENDIF(NOT SLEPc_INCLUDE_DIRS)

SET(SLEPc_FOUND FALSE)

IF(SLEPc_INCLUDE_DIRS)
  SET(SLEPc_LIBRARY_DIRS ${SLEPc_INCLUDE_DIRS})

  IF("${SLEPc_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(SLEPc_LIBRARY_DIRS ${SLEPc_LIBRARY_DIRS} PATH)
  ENDIF("${SLEPc_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${SLEPc_LIBRARY_DIRS}/lib")
    SET(SLEPc_LIBRARY_DIRS ${SLEPc_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${SLEPc_LIBRARY_DIRS}/lib")

  FIND_LIBRARY(SLEPc_DEBUG_LIBRARY   NAMES slepcd slepc_d libslepcd libslepc_d
               PATH_SUFFIXES "" Debug   PATHS ${SLEPc_LIBRARY_DIRS} NO_DEFAULT_PATH)
  FIND_LIBRARY(SLEPc_RELEASE_LIBRARY NAMES slepc libslepc
               PATH_SUFFIXES "" Release PATHS ${SLEPc_LIBRARY_DIRS} NO_DEFAULT_PATH)

  SET(SLEPc_LIBRARIES)
  IF(SLEPc_DEBUG_LIBRARY AND SLEPc_RELEASE_LIBRARY)
    SET(SLEPc_LIBRARIES debug ${SLEPc_DEBUG_LIBRARY} optimized ${SLEPc_RELEASE_LIBRARY})
  ELSEIF(SLEPc_DEBUG_LIBRARY)
    SET(SLEPc_LIBRARIES ${SLEPc_DEBUG_LIBRARY})
  ELSEIF(SLEPc_RELEASE_LIBRARY)
    SET(SLEPc_LIBRARIES ${SLEPc_RELEASE_LIBRARY})
  ENDIF(SLEPc_DEBUG_LIBRARY AND SLEPc_RELEASE_LIBRARY)

  IF(SLEPc_LIBRARIES)
    SET(SLEPc_FOUND TRUE)
  ENDIF(SLEPc_LIBRARIES)
ENDIF(SLEPc_INCLUDE_DIRS)

IF(SLEPc_FOUND)
  IF(NOT SLEPc_FIND_QUIETLY)
    MESSAGE(STATUS "Found SLEPc: headers at ${SLEPc_INCLUDE_DIRS}, libraries at ${SLEPc_LIBRARIES}")
  ENDIF(NOT SLEPc_FIND_QUIETLY)
ELSE(SLEPc_FOUND)
  IF(SLEPc_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "SLEPc not found")
  ENDIF(SLEPc_FIND_REQUIRED)
ENDIF(SLEPc_FOUND)
