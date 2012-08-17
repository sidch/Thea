# - Searches for an installation of the Mutil utilities library
#
# Defines:
#
#   Mutil_FOUND         True if Mutil was found, else false
#   Mutil_LIBRARIES     The libraries to link (currently blank)
#   Mutil_INCLUDE_DIRS  The directories containing the header files
#   Mutil_CFLAGS        Extra compiler flags
#   Mutil_LDFLAGS       Extra linker flags
#
# To specify an additional directory to search, set Mutil_ROOT.
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(Mutil_FOUND FALSE)
SET(Mutil_CFLAGS)
SET(Mutil_LDFLAGS)

# Look for the Mutil header, first in the user-specified location and then in the system locations
SET(Mutil_INCLUDE_DOC "The directory containing the Mutil include file Mutil/Mutil.hpp")
FIND_PATH(Mutil_INCLUDE_DIRS NAMES Mutil/Mutil.hpp PATHS ${Mutil_ROOT} ${Mutil_ROOT}/include ${Mutil_ROOT}/Source
          DOC ${Mutil_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT Mutil_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(Mutil_INCLUDE_DIRS NAMES Mutil/Mutil.hpp DOC ${Mutil_INCLUDE_DOC})
ENDIF(NOT Mutil_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(Mutil_INCLUDE_DIRS)
  SET(Mutil_FOUND TRUE)

  # Maybe there is also a library file
  SET(Mutil_LIBRARY_DIRS ${Mutil_INCLUDE_DIRS})
  IF("${Mutil_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mutil_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(Mutil_LIBRARY_DIRS ${Mutil_LIBRARY_DIRS} PATH)
  ENDIF("${Mutil_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mutil_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(Mutil_DEBUG_LIBRARY
               NAMES Mutil_d Mutild
               PATH_SUFFIXES "" Debug
               PATHS ${Mutil_LIBRARY_DIRS} ${Mutil_LIBRARY_DIRS}/lib ${Mutil_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Mutil_RELEASE_LIBRARY
               NAMES Mutil
               PATH_SUFFIXES "" Release
               PATHS ${Mutil_LIBRARY_DIRS} ${Mutil_LIBRARY_DIRS}/lib ${Mutil_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  SET(Mutil_LIBRARIES)
  IF(Mutil_DEBUG_LIBRARY AND Mutil_RELEASE_LIBRARY)
    SET(Mutil_LIBRARIES debug ${Mutil_DEBUG_LIBRARY} optimized ${Mutil_RELEASE_LIBRARY})
  ELSEIF(Mutil_DEBUG_LIBRARY)
    SET(Mutil_LIBRARIES ${Mutil_DEBUG_LIBRARY})
  ELSEIF(Mutil_RELEASE_LIBRARY)
    SET(Mutil_LIBRARIES ${Mutil_RELEASE_LIBRARY})
  ENDIF(Mutil_DEBUG_LIBRARY AND Mutil_RELEASE_LIBRARY)

  IF(Mutil_LIBRARIES)
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      SET(Mutil_CFLAGS "-DMUTIL_DLL")
    ELSE(WIN32)
      # Flags for importing symbols from dynamically linked libraries
      IF("${Mutil_LIBRARIES}" MATCHES ".dylib$" OR "${Mutil_LIBRARIES}" MATCHES ".so$")
        SET(Mutil_CFLAGS "-DMUTIL_DLL")
      ENDIF("${Mutil_LIBRARIES}" MATCHES ".dylib$" OR "${Mutil_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)
  ENDIF(Mutil_LIBRARIES)
ENDIF(Mutil_INCLUDE_DIRS)

SET(Mutil_CFLAGS ${Mutil_CFLAGS}  CACHE STRING "Extra compiler flags required by Mutil")
SET(Mutil_LDLAGS ${Mutil_LDFLAGS} CACHE STRING "Extra linker flags required by Mutil")

IF(Mutil_FOUND)
  IF(NOT Mutil_FIND_QUIETLY)
    IF(Mutil_LIBRARIES)
      MESSAGE(STATUS "Found Mutil: headers at ${Mutil_INCLUDE_DIRS}, libraries at ${Mutil_LIBRARIES}")
    ELSE(Mutil_LIBRARIES)
      MESSAGE(STATUS "Found Mutil: headers at ${Mutil_INCLUDE_DIRS} (no compiled library found)")
    ENDIF(Mutil_LIBRARIES)
  ENDIF(NOT Mutil_FIND_QUIETLY)
ELSE(Mutil_FOUND)
  IF(Mutil_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Mutil not found")
  ENDIF(Mutil_FIND_REQUIRED)
ENDIF(Mutil_FOUND)
