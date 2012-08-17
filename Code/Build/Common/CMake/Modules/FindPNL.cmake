# Searches for an installation of the Intel Probabilistic Networks library. On success, it sets the following variables:
#
#   PNL_FOUND              Set to true to indicate the PNL library was found
#   PNL_INCLUDE_DIRS       The directory containing the main PNL header file PNL/PNL.h
#   PNL_LIBRARY_DIRS       The directory containing the PNL libraries (for info only, should not be required)
#   PNL_LIBRARIES          The libraries needed to use PNL (with full paths)
#
# To specify an additional directory to search, set PNL_ROOT.
#
# Copyright (C) Siddhartha Chaudhuri, 2010
#

# Look for the PNL headers, first in the user-specified location and then in the system locations
SET(PNL_INCLUDE_DOC "The directory containing the PNL header file pnlGraphicalModel.hpp")
FIND_PATH(PNL_INCLUDE_DIRS NAMES pnlGraphicalModel.hpp PATHS ${PNL_ROOT} ${PNL_ROOT}/include DOC ${PNL_INCLUDE_DOC}
          NO_DEFAULT_PATH)
IF(NOT PNL_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(PNL_INCLUDE_DIRS NAMES pnlGraphicalModel.hpp DOC ${PNL_INCLUDE_DOC})
ENDIF(NOT PNL_INCLUDE_DIRS)

SET(PNL_FOUND FALSE)

IF(PNL_INCLUDE_DIRS)
  SET(PNL_LIBRARY_DIRS ${PNL_INCLUDE_DIRS})
  IF(EXISTS "${PNL_INCLUDE_DIRS}/opencx")
    SET(PNL_INCLUDE_DIRS ${PNL_INCLUDE_DIRS} "${PNL_INCLUDE_DIRS}/opencx")
  ENDIF(EXISTS "${PNL_INCLUDE_DIRS}/opencx")

  IF("${PNL_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(PNL_LIBRARY_DIRS ${PNL_LIBRARY_DIRS} PATH)
  ENDIF("${PNL_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${PNL_LIBRARY_DIRS}/lib")
    SET(PNL_LIBRARY_DIRS ${PNL_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${PNL_LIBRARY_DIRS}/lib")

  # Find PNL libraries
  FIND_LIBRARY(PNL_LIBRARY NAMES pnl PATHS ${PNL_LIBRARY_DIRS} NO_DEFAULT_PATH)

  IF(PNL_LIBRARY)
    SET(PNL_FOUND TRUE)
    SET(PNL_LIBRARIES ${PNL_LIBRARY})

    FIND_LIBRARY(PNL_CXCORE_LIBRARY NAMES cxcore PATHS ${PNL_LIBRARY_DIRS} NO_DEFAULT_PATH)
    IF(PNL_CXCORE_LIBRARY)
      SET(PNL_LIBRARIES ${PNL_LIBRARIES} ${PNL_CXCORE_LIBRARY})
    ENDIF(PNL_CXCORE_LIBRARY)

  ENDIF(PNL_LIBRARY)
ENDIF(PNL_INCLUDE_DIRS)

IF(PNL_FOUND)
  IF(NOT PNL_FIND_QUIETLY)
    MESSAGE(STATUS "Found PNL: headers at ${PNL_INCLUDE_DIRS}, libraries at ${PNL_LIBRARIES}")
  ENDIF(NOT PNL_FIND_QUIETLY)
ELSE(PNL_FOUND)
  IF(PNL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "PNL not found")
  ENDIF(PNL_FIND_REQUIRED)
ENDIF(PNL_FOUND)
