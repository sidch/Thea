# Searches for an installation of the LAPACK library. On success, it sets the following variables:
#
#   LAPACK_FOUND      Set to true to indicate the library was found
#   LAPACK_LIBRARIES  All libraries needed to use LAPACK (with full path)
#
# To specify an additional directory to search, set LAPACK_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

SET(LAPACK_FOUND FALSE)

# First look in user-provided root directory, then look in system locations
FIND_LIBRARY(LAPACK_LIBRARIES NAMES lapack liblapack LAPACK libLAPACK PATHS "${LAPACK_ROOT}" "${LAPACK_ROOT}/lib"
             NO_DEFAULT_PATH)
IF(NOT LAPACK_LIBRARIES)
  FIND_LIBRARY(LAPACK_LIBRARIES NAMES lapack liblapack LAPACK libLAPACK)
ENDIF(NOT LAPACK_LIBRARIES)

IF(LAPACK_LIBRARIES)
  SET(LAPACK_FOUND TRUE)
ENDIF(LAPACK_LIBRARIES)

IF(LAPACK_FOUND)
  IF(NOT LAPACK_FIND_QUIETLY)
    MESSAGE(STATUS "Found LAPACK: libraries at ${LAPACK_LIBRARIES}")
  ENDIF(NOT LAPACK_FIND_QUIETLY)
ELSE(LAPACK_FOUND)
  IF(LAPACK_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "LAPACK not found")
  ENDIF(LAPACK_FIND_REQUIRED)
ENDIF(LAPACK_FOUND)
