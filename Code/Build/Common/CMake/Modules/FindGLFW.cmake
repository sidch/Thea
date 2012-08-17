# Searches for an installation of the glfw library. On success, it sets the following variables:
#
#   GLFW_FOUND              Set to true to indicate the glfw library was found
#   GLFW_INCLUDE_DIRS       The directory containing the header file GL/glfw.h
#   GLFW_LIBRARIES          The libraries needed to use the glfw library
#
# To specify an additional directory to search, set GLFW_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

# Look for the header, first in the user-specified location and then in the system locations
SET(GLFW_INCLUDE_DOC "The directory containing the header file GL/glfw.h")
FIND_PATH(GLFW_INCLUDE_DIRS NAMES GL/glfw.h PATHS ${GLFW_ROOT} ${GLFW_ROOT}/include DOC ${GLFW_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT GLFW_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(GLFW_INCLUDE_DIRS NAMES GL/glfw.h DOC ${GLFW_INCLUDE_DOC})
ENDIF(NOT GLFW_INCLUDE_DIRS)

SET(GLFW_FOUND FALSE)

IF(GLFW_INCLUDE_DIRS)
  SET(GLFW_LIBRARY_DIRS ${GLFW_INCLUDE_DIRS})

  IF("${GLFW_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(GLFW_LIBRARY_DIRS ${GLFW_LIBRARY_DIRS} PATH)
  ENDIF("${GLFW_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${GLFW_LIBRARY_DIRS}/lib")
    SET(GLFW_LIBRARY_DIRS ${GLFW_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${GLFW_LIBRARY_DIRS}/lib")

  # Find GLFW libraries
  FIND_LIBRARY(GLFW_LIBRARIES NAMES glfw PATHS ${GLFW_LIBRARY_DIRS} NO_DEFAULT_PATH)

  IF(GLFW_LIBRARIES)
    SET(GLFW_FOUND TRUE)
  ENDIF(GLFW_LIBRARIES)
ENDIF(GLFW_INCLUDE_DIRS)

IF(GLFW_FOUND)
  IF(NOT GLFW_FIND_QUIETLY)
    MESSAGE(STATUS "Found GLFW: headers at ${GLFW_INCLUDE_DIRS}, libraries at ${GLFW_LIBRARIES}")
  ENDIF(NOT GLFW_FIND_QUIETLY)
ELSE(GLFW_FOUND)
  IF(GLFW_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "GLFW library not found")
  ENDIF(GLFW_FIND_REQUIRED)
ENDIF(GLFW_FOUND)
