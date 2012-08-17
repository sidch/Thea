# - Searches for an installation of the Cesium identifiers library
#
# Defines:
#
#   Cesium_FOUND         True if Cesium was found, else false
#   Cesium_LIBRARIES     Libraries to link
#   Cesium_INCLUDE_DIRS  The directories containing the header files
#   Cesium_CFLAGS        Extra compiler flags
#   Cesium_LDFLAGS       Extra linker flags
#
# To specify an additional directory to search, set Cesium_ROOT.
#
# Author: Ewen Cheslack-Postava, 2008
#

SET(Cesium_FOUND FALSE)
SET(Cesium_CFLAGS)
SET(Cesium_LDFLAGS)

# Look for a cesium header, first in the user-specified location and then in the system locations
SET(Cesium_INCLUDE_DOC "The directory containing the cesium include file pipeline.hpp")
FIND_PATH(Cesium_INCLUDE_DIRS NAMES cesium/pipeline.hpp PATHS ${cesium_ROOT} ${Cesium_ROOT}/include
          DOC ${Cesium_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT Cesium_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(Cesium_INCLUDE_DIRS NAMES cesium/pipeline.hpp DOC ${Cesium_INCLUDE_DOC})
ENDIF(NOT Cesium_INCLUDE_DIRS)

# Only look for the library file in the immediate neighborhood of the include directory
IF(Cesium_INCLUDE_DIRS)
  SET(Cesium_PROJECT_DIR ${Cesium_INCLUDE_DIRS})
  IF("${Cesium_PROJECT_DIR}" MATCHES "/include$")
    # Strip off the trailing "/include" from the path
    GET_FILENAME_COMPONENT(Cesium_PROJECT_DIR ${Cesium_PROJECT_DIR} PATH)
  ENDIF("${Cesium_PROJECT_DIR}" MATCHES "/include$")

  FIND_LIBRARY(Cesium_DEBUG_LIBRARY
               NAMES cesium_d cesiumd
               PATH_SUFFIXES "" Debug
               PATHS ${Cesium_PROJECT_DIR} ${Cesium_PROJECT_DIR}/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Cesium_RELEASE_LIBRARY
               NAMES cesium
               PATH_SUFFIXES "" Release
               PATHS ${Cesium_PROJECT_DIR} ${Cesium_PROJECT_DIR}/lib NO_DEFAULT_PATH)

  SET(Cesium_LIBRARIES)
  IF(Cesium_DEBUG_LIBRARY AND Cesium_RELEASE_LIBRARY)
    SET(Cesium_LIBRARIES debug ${Cesium_DEBUG_LIBRARY} optimized ${Cesium_RELEASE_LIBRARY})
  ELSEIF(Cesium_DEBUG_LIBRARY)
    SET(Cesium_LIBRARIES ${Cesium_DEBUG_LIBRARY})
  ELSEIF(Cesium_RELEASE_LIBRARY)
    SET(Cesium_LIBRARIES ${Cesium_RELEASE_LIBRARY})
  ENDIF(Cesium_DEBUG_LIBRARY AND Cesium_RELEASE_LIBRARY)

  IF(Cesium_LIBRARIES)
    SET(Cesium_FOUND TRUE)
  ENDIF(Cesium_LIBRARIES)
ENDIF(Cesium_INCLUDE_DIRS)

SET(Cesium_CFLAGS ${Cesium_CFLAGS}  CACHE STRING "Extra compiler flags required by cesium")
SET(Cesium_LDLAGS ${Cesium_LDFLAGS} CACHE STRING "Extra linker flags required by cesium")

IF(Cesium_FOUND)
  IF(NOT Cesium_FIND_QUIETLY)
    MESSAGE(STATUS "Found Cesium: headers at ${Cesium_INCLUDE_DIRS}, libraries at ${Cesium_LIBRARIES}")
  ENDIF(NOT Cesium_FIND_QUIETLY)
ELSE(Cesium_FOUND)
  IF(Cesium_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Cesium not found")
  ENDIF(Cesium_FIND_REQUIRED)
ENDIF(Cesium_FOUND)
