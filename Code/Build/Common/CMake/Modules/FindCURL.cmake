# Searches for an installation of the CURL library. On success, it sets the following variables:
#
#   CURL_FOUND              Set to true to indicate the library was found
#   CURL_INCLUDE_DIRS       The directory containing the main CURL header file curl/curl.h
#   CURL_LIBRARIES          The libraries needed to use CURL (without the full path)
#
# To specify an additional directory to search, set CURL_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2008
#

# Look for the CURL header, first in the user-specified location and then in the system locations
SET(CURL_INCLUDE_DOC "The directory containing the CURL header file curl/curl.h")
FIND_PATH(CURL_INCLUDE_DIRS NAMES curl/curl.h PATHS ${CURL_ROOT} ${CURL_ROOT}/include DOC ${CURL_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT CURL_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(CURL_INCLUDE_DIRS NAMES curl/curl.h DOC ${CURL_INCLUDE_DOC})
ENDIF(NOT CURL_INCLUDE_DIRS)

SET(CURL_FOUND FALSE)

IF(CURL_INCLUDE_DIRS)
  SET(CURL_LIBRARY_DIRS ${CURL_INCLUDE_DIRS})

  IF("${CURL_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(CURL_LIBRARY_DIRS ${CURL_LIBRARY_DIRS} PATH)
  ENDIF("${CURL_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${CURL_LIBRARY_DIRS}/lib")
    SET(CURL_LIBRARY_DIRS ${CURL_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${CURL_LIBRARY_DIRS}/lib")

  FIND_LIBRARY(CURL_DEBUG_LIBRARY   NAMES curld curl_d libcurld libcurl_d
               PATH_SUFFIXES "" Debug   PATHS ${CURL_LIBRARY_DIRS} NO_DEFAULT_PATH)
  FIND_LIBRARY(CURL_RELEASE_LIBRARY NAMES curl libcurl
               PATH_SUFFIXES "" Release PATHS ${CURL_LIBRARY_DIRS} NO_DEFAULT_PATH)

  SET(CURL_LIBRARIES)
  IF(CURL_DEBUG_LIBRARY AND CURL_RELEASE_LIBRARY)
    SET(CURL_LIBRARIES debug ${CURL_DEBUG_LIBRARY} optimized ${CURL_RELEASE_LIBRARY})
  ELSEIF(CURL_DEBUG_LIBRARY)
    SET(CURL_LIBRARIES ${CURL_DEBUG_LIBRARY})
  ELSEIF(CURL_RELEASE_LIBRARY)
    SET(CURL_LIBRARIES ${CURL_RELEASE_LIBRARY})
  ENDIF(CURL_DEBUG_LIBRARY AND CURL_RELEASE_LIBRARY)

  IF(CURL_LIBRARIES)
    SET(CURL_FOUND TRUE)
  ENDIF(CURL_LIBRARIES)
ENDIF(CURL_INCLUDE_DIRS)

IF(CURL_FOUND)
  IF(NOT CURL_FIND_QUIETLY)
    MESSAGE(STATUS "Found CURL: headers at ${CURL_INCLUDE_DIRS}, libraries at ${CURL_LIBRARIES}")
  ENDIF(NOT CURL_FIND_QUIETLY)
ELSE(CURL_FOUND)
  IF(CURL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "CURL not found")
  ENDIF(CURL_FIND_REQUIRED)
ENDIF(CURL_FOUND)
