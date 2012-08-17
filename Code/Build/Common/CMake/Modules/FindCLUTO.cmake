# Searches for an installation of the CLUTO libraries. On success, it sets the following variables:
#
#   CLUTO_FOUND              Set to true to indicate the CLUTO library was found
#   CLUTO_INCLUDE_DIRS       The directory containing the main CLUTO header file CLUTO/CLUTO.h
#   CLUTO_LIBRARY_DIRS       The directory containing the CLUTO libraries (for info only, should not be required)
#   CLUTO_LIBRARIES          The libraries needed to use CLUTO (with full paths)
#
# To specify an additional directory to search, set CLUTO_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

# Look for the CLUTO header, first in the user-specified location and then in the system locations
SET(CLUTO_INCLUDE_DOC "The directory containing the CLUTO header file cluto.h")
FIND_PATH(CLUTO_INCLUDE_DIRS NAMES cluto.h PATHS ${CLUTO_ROOT} ${CLUTO_ROOT}/include DOC ${CLUTO_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT CLUTO_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(CLUTO_INCLUDE_DIRS NAMES cluto.h DOC ${CLUTO_INCLUDE_DOC})
ENDIF(NOT CLUTO_INCLUDE_DIRS)

SET(CLUTO_FOUND FALSE)

IF(CLUTO_INCLUDE_DIRS)
  SET(CLUTO_LIBRARY_DIRS ${CLUTO_INCLUDE_DIRS})

  IF("${CLUTO_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(CLUTO_LIBRARY_DIRS ${CLUTO_LIBRARY_DIRS} PATH)
  ENDIF("${CLUTO_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${CLUTO_LIBRARY_DIRS}/lib")
    SET(CLUTO_LIBRARY_DIRS ${CLUTO_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${CLUTO_LIBRARY_DIRS}/lib")

  # Find CLUTO libraries
  FIND_LIBRARY(CLUTO_LIBRARIES NAMES cluto libcluto PATHS ${CLUTO_LIBRARY_DIRS} NO_DEFAULT_PATH)

  IF(CLUTO_LIBRARIES)
    SET(CLUTO_FOUND TRUE)
  ENDIF(CLUTO_LIBRARIES)
ENDIF(CLUTO_INCLUDE_DIRS)

IF(CLUTO_FOUND)
  IF(NOT CLUTO_FIND_QUIETLY)
    MESSAGE(STATUS "Found CLUTO: headers at ${CLUTO_INCLUDE_DIRS}, libraries at ${CLUTO_LIBRARIES}")
  ENDIF(NOT CLUTO_FIND_QUIETLY)
ELSE(CLUTO_FOUND)
  IF(CLUTO_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "CLUTO not found")
  ENDIF(CLUTO_FIND_REQUIRED)
ENDIF(CLUTO_FOUND)
