# - Searches for an installation of the Mid identifiers library
#
# Defines:
#
#   Mid_FOUND         True if Mid was found, else false
#   Mid_LIBRARIES     Libraries to link
#   Mid_INCLUDE_DIRS  The directories containing the header files
#   Mid_CFLAGS        Extra compiler flags
#   Mid_LDFLAGS       Extra linker flags
#
# To specify an additional directory to search, set Mid_ROOT.
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(Mid_FOUND FALSE)
SET(Mid_CFLAGS)
SET(Mid_LDFLAGS)

# Look for the Mid header, first in the user-specified location and then in the system locations
SET(Mid_INCLUDE_DOC "The directory containing the Mid include file Mid/Mid.hpp")
FIND_PATH(Mid_INCLUDE_DIRS NAMES Mid/Mid.hpp PATHS ${Mid_ROOT} ${Mid_ROOT}/include ${Mid_ROOT}/Source
          DOC ${Mid_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT Mid_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(Mid_INCLUDE_DIRS NAMES Mid/Mid.hpp DOC ${Mid_INCLUDE_DOC})
ENDIF(NOT Mid_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(Mid_INCLUDE_DIRS)
  SET(Mid_LIBRARY_DIRS ${Mid_INCLUDE_DIRS})
  IF("${Mid_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mid_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(Mid_LIBRARY_DIRS ${Mid_LIBRARY_DIRS} PATH)
  ENDIF("${Mid_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mid_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(Mid_DEBUG_LIBRARY
               NAMES Mid_d Midd
               PATH_SUFFIXES "" Debug
               PATHS ${Mid_LIBRARY_DIRS} ${Mid_LIBRARY_DIRS}/lib ${Mid_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Mid_RELEASE_LIBRARY
               NAMES Mid
               PATH_SUFFIXES "" Release
               PATHS ${Mid_LIBRARY_DIRS} ${Mid_LIBRARY_DIRS}/lib ${Mid_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  SET(Mid_LIBRARIES)
  IF(Mid_DEBUG_LIBRARY AND Mid_RELEASE_LIBRARY)
    SET(Mid_LIBRARIES debug ${Mid_DEBUG_LIBRARY} optimized ${Mid_RELEASE_LIBRARY})
  ELSEIF(Mid_DEBUG_LIBRARY)
    SET(Mid_LIBRARIES ${Mid_DEBUG_LIBRARY})
  ELSEIF(Mid_RELEASE_LIBRARY)
    SET(Mid_LIBRARIES ${Mid_RELEASE_LIBRARY})
  ENDIF(Mid_DEBUG_LIBRARY AND Mid_RELEASE_LIBRARY)

  IF(Mid_LIBRARIES)
    SET(Mid_FOUND TRUE)

    # Flags for importing symbols from dynamically linked libraries
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      #SET(Mid_CFLAGS "-DMID_DLL") # using static libs on windows for now
    ELSE(WIN32)
      IF("${Mid_LIBRARIES}" MATCHES ".dylib$" OR "${Mid_LIBRARIES}" MATCHES ".so$")
        SET(Mid_CFLAGS "-DMID_DLL")
      ENDIF("${Mid_LIBRARIES}" MATCHES ".dylib$" OR "${Mid_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)
  ENDIF(Mid_LIBRARIES)
ENDIF(Mid_INCLUDE_DIRS)

SET(Mid_CFLAGS ${Mid_CFLAGS}  CACHE STRING "Extra compiler flags required by Mid")
SET(Mid_LDLAGS ${Mid_LDFLAGS} CACHE STRING "Extra linker flags required by Mid")

IF(Mid_FOUND)
  IF(NOT Mid_FIND_QUIETLY)
    MESSAGE(STATUS "Found Mid: headers at ${Mid_INCLUDE_DIRS}, libraries at ${Mid_LIBRARIES}")
  ENDIF(NOT Mid_FIND_QUIETLY)
ELSE(Mid_FOUND)
  IF(Mid_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Mid not found")
  ENDIF(Mid_FIND_REQUIRED)
ENDIF(Mid_FOUND)
