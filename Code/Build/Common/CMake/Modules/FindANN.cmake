# Searches for an installation of the ANN (Arya/Mount's Approximate Nearest Neighbors) library. On success, it sets the
# following variables:
#
#   ANN_FOUND              Set to true to indicate the library was found
#   ANN_INCLUDE_DIRS       The directory containing the main ANN header file ANN/ANN.h
#   ANN_LIBRARIES          The libraries needed to use ANN (without the full path)
#
# To specify an additional directory to search, set ANN_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2010
#

# Look for the ANN header, first in the user-specified location and then in the system locations
SET(ANN_INCLUDE_DOC "The directory containing the ANN header file ANN/ANN.h")
FIND_PATH(ANN_INCLUDE_DIRS NAMES ANN/ANN.h PATHS ${ANN_ROOT} ${ANN_ROOT}/include DOC ${ANN_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT ANN_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(ANN_INCLUDE_DIRS NAMES ANN/ANN.h DOC ${ANN_INCLUDE_DOC})
ENDIF(NOT ANN_INCLUDE_DIRS)

SET(ANN_FOUND FALSE)

IF(ANN_INCLUDE_DIRS)
  SET(ANN_LIBRARY_DIRS ${ANN_INCLUDE_DIRS})

  IF("${ANN_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(ANN_LIBRARY_DIRS ${ANN_LIBRARY_DIRS} PATH)
  ENDIF("${ANN_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${ANN_LIBRARY_DIRS}/lib")
    SET(ANN_LIBRARY_DIRS ${ANN_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${ANN_LIBRARY_DIRS}/lib")

  FIND_LIBRARY(ANN_DEBUG_LIBRARY   NAMES ANNd ANN_d libANNd libANN_d
               PATH_SUFFIXES "" Debug   PATHS ${ANN_LIBRARY_DIRS} NO_DEFAULT_PATH)
  FIND_LIBRARY(ANN_RELEASE_LIBRARY NAMES ANN libANN
               PATH_SUFFIXES "" Release PATHS ${ANN_LIBRARY_DIRS} NO_DEFAULT_PATH)

  SET(ANN_LIBRARIES)
  IF(ANN_DEBUG_LIBRARY AND ANN_RELEASE_LIBRARY)
    SET(ANN_LIBRARIES debug ${ANN_DEBUG_LIBRARY} optimized ${ANN_RELEASE_LIBRARY})
  ELSEIF(ANN_DEBUG_LIBRARY)
    SET(ANN_LIBRARIES ${ANN_DEBUG_LIBRARY})
  ELSEIF(ANN_RELEASE_LIBRARY)
    SET(ANN_LIBRARIES ${ANN_RELEASE_LIBRARY})
  ENDIF(ANN_DEBUG_LIBRARY AND ANN_RELEASE_LIBRARY)

  IF(ANN_LIBRARIES)
    SET(ANN_FOUND TRUE)
  ENDIF(ANN_LIBRARIES)
ENDIF(ANN_INCLUDE_DIRS)

IF(ANN_FOUND)
  IF(NOT ANN_FIND_QUIETLY)
    MESSAGE(STATUS "Found ANN: headers at ${ANN_INCLUDE_DIRS}, libraries at ${ANN_LIBRARIES}")
  ENDIF(NOT ANN_FIND_QUIETLY)
ELSE(ANN_FOUND)
  IF(ANN_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ANN not found")
  ENDIF(ANN_FIND_REQUIRED)
ENDIF(ANN_FOUND)
