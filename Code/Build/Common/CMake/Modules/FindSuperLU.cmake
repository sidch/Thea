# Searches for an installation of the SuperLU library. On success, it sets the following variables:
#
#   SuperLU_FOUND         Set to true to indicate the library was found
#   SuperLU_INCLUDE_DIRS  Location of SuperLU header files (not set if SuperLU_LIBS_ONLY is true -- see below)
#   SuperLU_LIBRARIES     Location of SuperLU libraries (with full path)
#
# To specify an additional directory to search, set SuperLU_ROOT.
# To ignore looking for the headers (e.g. for ARPACK++, which has its own versions of the headers) set SuperLU_LIBS_ONLY
#
# Copyright (C) Siddhartha Chaudhuri, 2009
#

SET(SuperLU_FOUND FALSE)

IF(SuperLU_LIBS_ONLY)

  # First look for libsuperlu in user-provided root directory, then look in system locations
  FIND_LIBRARY(SuperLU_LIBRARIES NAMES superlu libsuperlu SuperLU libSuperLU PATHS "${SuperLU_ROOT}" "${SuperLU_ROOT}/lib"
               NO_DEFAULT_PATH)
  IF(NOT SuperLU_LIBRARIES)
    FIND_LIBRARY(SuperLU_LIBRARIES NAMES superlu libsuperlu SuperLU libSuperLU)
  ENDIF(NOT SuperLU_LIBRARIES)

  IF(SuperLU_LIBRARIES)
    SET(SuperLU_FOUND TRUE)
  ENDIF(SuperLU_LIBRARIES)

ELSEIF(SuperLU_LIBS_ONLY)

  # Look for the SuperLU header, first in the user-specified location and then in the system locations
  SET(SuperLU_INCLUDE_DOC "The directory containing the SuperLU header file slu_util.h")
  FIND_PATH(SuperLU_INCLUDE_DIRS NAMES slu_util.h PATHS ${SuperLU_ROOT} ${SuperLU_ROOT}/include DOC ${SuperLU_INCLUDE_DOC}
            NO_DEFAULT_PATH)
  IF(NOT SuperLU_INCLUDE_DIRS)  # now look in system locations
    FIND_PATH(SuperLU_INCLUDE_DIRS NAMES slu_util.h DOC ${SuperLU_INCLUDE_DOC})
  ENDIF(NOT SuperLU_INCLUDE_DIRS)

  IF(SuperLU_INCLUDE_DIRS)
    SET(SuperLU_LIBRARY_DIRS ${SuperLU_INCLUDE_DIRS})

    IF("${SuperLU_LIBRARY_DIRS}" MATCHES "/include$")
      # Strip off the trailing "/include" in the path.
      GET_FILENAME_COMPONENT(SuperLU_LIBRARY_DIRS ${SuperLU_LIBRARY_DIRS} PATH)
    ENDIF("${SuperLU_LIBRARY_DIRS}" MATCHES "/include$")

    IF(EXISTS "${SuperLU_LIBRARY_DIRS}/lib")
      SET(SuperLU_LIBRARY_DIRS ${SuperLU_LIBRARY_DIRS}/lib)
    ENDIF(EXISTS "${SuperLU_LIBRARY_DIRS}/lib")

    # Find SuperLU library
    FIND_LIBRARY(SuperLU_LIBRARIES NAMES superlu libsuperlu SuperLU libSuperLU PATHS ${SuperLU_LIBRARY_DIRS} NO_DEFAULT_PATH)

    IF(SuperLU_LIBRARIES)
      SET(SuperLU_FOUND TRUE)
    ENDIF(SuperLU_LIBRARIES)
  ENDIF(SuperLU_INCLUDE_DIRS)

ENDIF(SuperLU_LIBS_ONLY)

IF(SuperLU_FOUND)
  IF(NOT SuperLU_FIND_QUIETLY)
    IF(SuperLU_LIBS_ONLY)
      MESSAGE(STATUS "Found SuperLU: libraries at ${SuperLU_LIBRARIES}")
    ELSE(SuperLU_LIBS_ONLY)
      MESSAGE(STATUS "Found SuperLU: headers at ${SuperLU_INCLUDE_DIRS}, libraries at ${SuperLU_LIBRARIES}")
    ENDIF(SuperLU_LIBS_ONLY)
  ENDIF(NOT SuperLU_FIND_QUIETLY)
ELSE(SuperLU_FOUND)
  IF(SuperLU_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "SuperLU not found")
  ENDIF(SuperLU_FIND_REQUIRED)
ENDIF(SuperLU_FOUND)
