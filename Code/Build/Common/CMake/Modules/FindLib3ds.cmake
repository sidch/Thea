# Searches for an installation of the lib3ds library (http://www.lib3ds.org)
#
# Defines:
#
#   Lib3ds_FOUND          True if Lib3ds was found, else false
#   Lib3ds_LIBRARIES      Libraries to link
#   Lib3ds_INCLUDE_DIRS   The directories containing the header files
#   Lib3ds_VERSION_MAJOR  The major version of the library
#
# To specify an additional directory to search, set Lib3ds_ROOT.
#
# Author: Siddhartha Chaudhuri, 2009
#

SET(Lib3ds_FOUND FALSE)

IF(NOT Lib3ds_VERSION_MAJOR)
  SET(Lib3ds_VERSION_MAJOR 1 CACHE INTERNAL "Major version number of lib3ds")
ENDIF(NOT Lib3ds_VERSION_MAJOR)

IF(NOT Lib3ds_INCLUDE_DIRS)
  # Look for the Lib3ds header, first in the user-specified location and then in the system locations
  SET(Lib3ds_INCLUDE_DOC
      "The directory containing the lib3ds include files (pre-version 2: lib3ds/mesh.h; version 2 and later: lib3ds.h)")
  FIND_PATH(Lib3ds_INCLUDE_DIRS NAMES lib3ds.h PATHS ${Lib3ds_ROOT} ${Lib3ds_ROOT}/include DOC ${Lib3ds_INCLUDE_DOC}
            NO_DEFAULT_PATH)
  IF(NOT Lib3ds_INCLUDE_DIRS)  # now look in system locations
    FIND_PATH(Lib3ds_INCLUDE_DIRS NAMES lib3ds.h DOC ${Lib3ds_INCLUDE_DOC})
  ENDIF(NOT Lib3ds_INCLUDE_DIRS)

  IF(Lib3ds_INCLUDE_DIRS)
    # Version 2 onwards has unified header
    SET(Lib3ds_VERSION_MAJOR 2 CACHE INTERNAL "Major version number of lib3ds")
  ELSE(Lib3ds_INCLUDE_DIRS)
    FIND_PATH(Lib3ds_INCLUDE_DIRS NAMES lib3ds/types.h PATHS ${Lib3ds_ROOT} ${Lib3ds_ROOT}/include DOC ${Lib3ds_INCLUDE_DOC}
              NO_DEFAULT_PATH)
    IF(NOT Lib3ds_INCLUDE_DIRS)  # now look in system locations
      FIND_PATH(Lib3ds_INCLUDE_DIRS NAMES lib3ds/types.h DOC ${Lib3ds_INCLUDE_DOC})
    ENDIF(NOT Lib3ds_INCLUDE_DIRS)
  ENDIF(Lib3ds_INCLUDE_DIRS)
ENDIF(NOT Lib3ds_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(Lib3ds_INCLUDE_DIRS)
  SET(Lib3ds_LIBRARY_DIRS ${Lib3ds_INCLUDE_DIRS})
  IF("${Lib3ds_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" from the path
    GET_FILENAME_COMPONENT(Lib3ds_LIBRARY_DIRS ${Lib3ds_LIBRARY_DIRS} PATH)
  ENDIF("${Lib3ds_LIBRARY_DIRS}" MATCHES "/include$")

  FIND_LIBRARY(Lib3ds_DEBUG_LIBRARY NAMES 3ds_d lib3ds_d 3dsd lib3dsd PATH_SUFFIXES "" Debug
               PATHS ${Lib3ds_LIBRARY_DIRS} ${Lib3ds_LIBRARY_DIRS}/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Lib3ds_RELEASE_LIBRARY NAMES 3ds lib3ds PATH_SUFFIXES "" Release
               PATHS ${Lib3ds_LIBRARY_DIRS} ${Lib3ds_LIBRARY_DIRS}/lib NO_DEFAULT_PATH)

  SET(Lib3ds_LIBRARIES)
  IF(Lib3ds_DEBUG_LIBRARY AND Lib3ds_RELEASE_LIBRARY)
    SET(Lib3ds_LIBRARIES debug ${Lib3ds_DEBUG_LIBRARY} optimized ${Lib3ds_RELEASE_LIBRARY})
  ELSEIF(Lib3ds_DEBUG_LIBRARY)
    SET(Lib3ds_LIBRARIES ${Lib3ds_DEBUG_LIBRARY})
  ELSEIF(Lib3ds_RELEASE_LIBRARY)
    SET(Lib3ds_LIBRARIES ${Lib3ds_RELEASE_LIBRARY})
  ENDIF(Lib3ds_DEBUG_LIBRARY AND Lib3ds_RELEASE_LIBRARY)

  IF(Lib3ds_LIBRARIES)
    SET(Lib3ds_FOUND TRUE)
  ENDIF(Lib3ds_LIBRARIES)
ENDIF(Lib3ds_INCLUDE_DIRS)

IF(Lib3ds_FOUND)
  IF(NOT Lib3ds_FIND_QUIETLY)
    MESSAGE(STATUS "Found lib3ds: headers at ${Lib3ds_INCLUDE_DIRS}, libraries at ${Lib3ds_LIBRARIES}")
  ENDIF(NOT Lib3ds_FIND_QUIETLY)
ELSE(Lib3ds_FOUND)
  IF(Lib3ds_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "lib3ds not found")
  ENDIF(Lib3ds_FIND_REQUIRED)
ENDIF(Lib3ds_FOUND)
