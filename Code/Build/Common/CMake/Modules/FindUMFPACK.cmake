# Searches for an installation of the UMFPACK library. On success, it sets the following variables:
#
#   UMFPACK_FOUND         Set to true to indicate the library was found
#   UMFPACK_INCLUDE_DIRS  Location of UMFPACK header files (not set if UMFPACK_LIBS_ONLY is true -- see below)
#   UMFPACK_LIBRARIES     Location of UMFPACK libraries (with full path)
#
# To specify an additional directory to search, set UMFPACK_ROOT.
# To ignore looking for the headers (e.g. for ARPACK++, which has its own versions of the headers) set UMFPACK_LIBS_ONLY
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

SET(UMFPACK_FOUND FALSE)

IF(UMFPACK_LIBS_ONLY)

  # First look for libumfpack in user-provided root directory, then look in system locations
  FIND_LIBRARY(UMFPACK_LIBRARIES NAMES umfpack libumfpack UMFPACK libUMFPACK PATHS "${UMFPACK_ROOT}" "${UMFPACK_ROOT}/lib"
               NO_DEFAULT_PATH)
  IF(NOT UMFPACK_LIBRARIES)
    FIND_LIBRARY(UMFPACK_LIBRARIES NAMES umfpack libumfpack UMFPACK libUMFPACK)
  ENDIF(NOTUMFPACK_LIBRARIES)

  IF(UMFPACK_LIBRARIES)
    SET(UMFPACK_FOUND TRUE)
  ENDIF(UMFPACK_LIBRARIES)

ELSEIF(UMFPACK_LIBS_ONLY)

  # Look for the UMFPACK header, first in the user-specified location and then in the system locations
  SET(UMFPACK_INCLUDE_DOC "The directory containing the UMFPACK header file umfpack.h")
  FIND_PATH(UMFPACK_INCLUDE_DIRS NAMES umfpack.h PATHS ${UMFPACK_ROOT} ${UMFPACK_ROOT}/include DOC ${UMFPACK_INCLUDE_DOC}
            NO_DEFAULT_PATH)
  IF(NOT UMFPACK_INCLUDE_DIRS)  # now look in system locations
    FIND_PATH(UMFPACK_INCLUDE_DIRS NAMES umfpack.h DOC ${UMFPACK_INCLUDE_DOC})
  ENDIF(NOT UMFPACK_INCLUDE_DIRS)

  IF(UMFPACK_INCLUDE_DIRS)
    SET(UMFPACK_LIBRARY_DIRS ${UMFPACK_INCLUDE_DIRS})

    IF("${UMFPACK_LIBRARY_DIRS}" MATCHES "/include$")
      # Strip off the trailing "/include" in the path.
      GET_FILENAME_COMPONENT(UMFPACK_LIBRARY_DIRS ${UMFPACK_LIBRARY_DIRS} PATH)
    ENDIF("${UMFPACK_LIBRARY_DIRS}" MATCHES "/include$")

    IF(EXISTS "${UMFPACK_LIBRARY_DIRS}/lib")
      SET(UMFPACK_LIBRARY_DIRS ${UMFPACK_LIBRARY_DIRS}/lib)
    ENDIF(EXISTS "${UMFPACK_LIBRARY_DIRS}/lib")

    # Find UMFPACK library
    FIND_LIBRARY(UMFPACK_LIBRARIES NAMES umfpack libumfpack UMFPACK libUMFPACK PATHS ${UMFPACK_LIBRARY_DIRS} NO_DEFAULT_PATH)

    # If we find an AMD installation in the same directory, it's probably required by UMFPACK
    FIND_LIBRARY(UMFPACK_AMD_LIBRARY NAMES amd AMD libamd libAMD PATHS ${UMFPACK_LIBRARY_DIRS} NO_DEFAULT_PATH)
    SET(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${UMFPACK_AMD_LIBRARY})

    IF(UMFPACK_LIBRARIES)
      SET(UMFPACK_FOUND TRUE)
    ENDIF(UMFPACK_LIBRARIES)
  ENDIF(UMFPACK_INCLUDE_DIRS)

ENDIF(UMFPACK_LIBS_ONLY)

IF(UMFPACK_FOUND)
  IF(NOT UMFPACK_FIND_QUIETLY)
    IF(UMFPACK_LIBS_ONLY)
      MESSAGE(STATUS "Found UMFPACK: libraries at ${UMFPACK_LIBRARIES}")
    ELSE(UMFPACK_LIBS_ONLY)
      MESSAGE(STATUS "Found UMFPACK: headers at ${UMFPACK_INCLUDE_DIRS}, libraries at ${UMFPACK_LIBRARIES}")
    ENDIF(UMFPACK_LIBS_ONLY)
  ENDIF(NOT UMFPACK_FIND_QUIETLY)
ELSE(UMFPACK_FOUND)
  IF(UMFPACK_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "UMFPACK not found")
  ENDIF(UMFPACK_FIND_REQUIRED)
ENDIF(UMFPACK_FOUND)
