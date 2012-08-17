# Searches for an installation of the ZLIB library. On success, it sets the following variables:
#
#   ZLIB_FOUND              Set to true to indicate the zip library was found
#   ZLIB_INCLUDE_DIRS       The directory containing the header file zlib.h
#   ZLIB_LIBRARIES          The libraries needed to use the ZLIB library
#
# To specify an additional directory to search, set ZLIB_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

# Look for the header, first in the user-specified location and then in the system locations
SET(ZLIB_INCLUDE_DOC "The directory containing the header file zlib.h")
FIND_PATH(ZLIB_INCLUDE_DIRS NAMES zlib.h PATHS ${ZLIB_ROOT} ${ZLIB_ROOT}/include DOC ${ZLIB_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT ZLIB_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(ZLIB_INCLUDE_DIRS NAMES zlib.h DOC ${ZLIB_INCLUDE_DOC})
ENDIF(NOT ZLIB_INCLUDE_DIRS)

SET(ZLIB_FOUND FALSE)

IF(ZLIB_INCLUDE_DIRS)
  SET(ZLIB_LIBRARY_DIRS ${ZLIB_INCLUDE_DIRS})

  IF("${ZLIB_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(ZLIB_LIBRARY_DIRS ${ZLIB_LIBRARY_DIRS} PATH)
  ENDIF("${ZLIB_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${ZLIB_LIBRARY_DIRS}/lib")
    SET(ZLIB_LIBRARY_DIRS ${ZLIB_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${ZLIB_LIBRARY_DIRS}/lib")

  # Find ZLIB libraries
  FIND_LIBRARY(ZLIB_LIBRARIES NAMES z zlib zdll PATHS ${ZLIB_LIBRARY_DIRS} NO_DEFAULT_PATH)
  IF(ZLIB_LIBRARIES)
    SET(ZLIB_FOUND TRUE)
  ENDIF(ZLIB_LIBRARIES)
ENDIF(ZLIB_INCLUDE_DIRS)

IF(ZLIB_FOUND)
  IF(NOT ZLIB_FIND_QUIETLY)
    MESSAGE(STATUS "Found ZLIB: headers at ${ZLIB_INCLUDE_DIRS}, libraries at ${ZLIB_LIBRARY_DIRS}")
  ENDIF(NOT ZLIB_FIND_QUIETLY)
ELSE(ZLIB_FOUND)
  IF(ZLIB_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ZLIB not found")
  ENDIF(ZLIB_FIND_REQUIRED)
ENDIF(ZLIB_FOUND)
