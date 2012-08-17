# Searches for an installation of the BLAS library. On success, it sets the following variables:
#
#   BLAS_FOUND      Set to true to indicate the library was found
#   BLAS_LIBRARIES  All libraries needed to use BLAS (with full path)
#
# To specify an additional directory to search, set BLAS_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

SET(BLAS_FOUND FALSE)

# First look in user-provided root directory, then look in system locations
FIND_LIBRARY(BLAS_LIBRARIES NAMES lapack liblapack BLAS libBLAS PATHS "${BLAS_ROOT}" "${BLAS_ROOT}/lib"
             NO_DEFAULT_PATH)
IF(NOT BLAS_LIBRARIES)
  FIND_LIBRARY(BLAS_LIBRARIES NAMES lapack liblapack BLAS libBLAS)
ENDIF(NOT BLAS_LIBRARIES)

IF(BLAS_LIBRARIES)
  SET(BLAS_FOUND TRUE)
ENDIF(BLAS_LIBRARIES)

IF(BLAS_FOUND)
  IF(NOT BLAS_FIND_QUIETLY)
    MESSAGE(STATUS "Found BLAS: libraries at ${BLAS_LIBRARIES}")
  ENDIF(NOT BLAS_FIND_QUIETLY)
ELSE(BLAS_FOUND)
  IF(BLAS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "BLAS not found")
  ENDIF(BLAS_FIND_REQUIRED)
ENDIF(BLAS_FOUND)
