# - Searches for an installation of the SNL library
#
# Defines:
#
#   SNL_FOUND           True if SNL was found, else false
#   SNL_LIBRARIES       Libraries to link
#   SNL_INCLUDE_DIRS    The directories containing the header files
#   SNL_CFLAGS          Extra compiler flags
#
# To specify an additional directory to search, set SNL_ROOT.
#
# Author: Siddhartha Chaudhuri, 2012
#

SET(SNL_FOUND FALSE)
SET(SNL_CFLAGS )

# Look for the SNL header, first in the user-specified location and then in the system locations
SET(SNL_INCLUDE_DOC "The directory containing the SNL include file SNL/SNL.hpp")
FIND_PATH(SNL_INCLUDE_DIRS NAMES SNL/Common.hpp PATHS ${SNL_ROOT} ${SNL_ROOT}/include ${SNL_ROOT}/Source
          DOC ${SNL_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT SNL_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(SNL_INCLUDE_DIRS NAMES SNL/Common.hpp DOC ${SNL_INCLUDE_DOC})
ENDIF(NOT SNL_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(SNL_INCLUDE_DIRS)
  SET(SNL_LIBRARY_DIRS ${SNL_INCLUDE_DIRS})
  IF("${SNL_LIBRARY_DIRS}" MATCHES "/include$" OR "${SNL_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(SNL_LIBRARY_DIRS ${SNL_LIBRARY_DIRS} PATH)
  ENDIF("${SNL_LIBRARY_DIRS}" MATCHES "/include$" OR "${SNL_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(SNL_DEBUG_LIBRARY
               NAMES SNL_d SNLd
               PATH_SUFFIXES "" Debug
               PATHS ${SNL_LIBRARY_DIRS} ${SNL_LIBRARY_DIRS}/lib ${SNL_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(SNL_RELEASE_LIBRARY
               NAMES SNL
               PATH_SUFFIXES "" Release
               PATHS ${SNL_LIBRARY_DIRS} ${SNL_LIBRARY_DIRS}/lib ${SNL_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  SET(SNL_LIBRARIES)
  IF(SNL_DEBUG_LIBRARY AND SNL_RELEASE_LIBRARY)
    SET(SNL_LIBRARIES debug ${SNL_DEBUG_LIBRARY} optimized ${SNL_RELEASE_LIBRARY})
  ELSEIF(SNL_DEBUG_LIBRARY)
    SET(SNL_LIBRARIES ${SNL_DEBUG_LIBRARY})
  ELSEIF(SNL_RELEASE_LIBRARY)
    SET(SNL_LIBRARIES ${SNL_RELEASE_LIBRARY})
  ENDIF(SNL_DEBUG_LIBRARY AND SNL_RELEASE_LIBRARY)

  IF(SNL_LIBRARIES)
    SET(SNL_FOUND TRUE)

    # Flags for importing symbols from dynamically linked libraries
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      SET(SNL_CFLAGS "-DSNL_DLL -DSNL_DLL_IMPORTS")
    ELSE(WIN32)
      IF("${SNL_LIBRARIES}" MATCHES ".dylib$" OR "${SNL_LIBRARIES}" MATCHES ".so$")
        SET(SNL_CFLAGS "-DSNL_DLL -DSNL_DLL_IMPORTS")
      ENDIF("${SNL_LIBRARIES}" MATCHES ".dylib$" OR "${SNL_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)

    # Read extra flags to be used to build SNL
    SET(SNL_BUILD_FLAGS_FILE "${SNL_INCLUDE_DIRS}/SNL/BuildFlags.txt")
    IF(EXISTS "${SNL_BUILD_FLAGS_FILE}")
      FILE(READ "${SNL_BUILD_FLAGS_FILE}" SNL_BUILD_FLAGS)
      STRING(REGEX REPLACE "\n" " " SNL_BUILD_FLAGS "${SNL_BUILD_FLAGS}")
      SET(SNL_CFLAGS "${SNL_CFLAGS} ${SNL_BUILD_FLAGS}")
    ENDIF(EXISTS "${SNL_BUILD_FLAGS_FILE}")

  ENDIF(SNL_LIBRARIES)
ENDIF(SNL_INCLUDE_DIRS)

# Remove duplicate entries from lists, else the same dirs and flags can repeat many times
IF(SNL_LIBRARIES)
  LIST(REMOVE_DUPLICATES SNL_LIBRARIES)
ENDIF(SNL_LIBRARIES)

IF(SNL_INCLUDE_DIRS)
  LIST(REMOVE_DUPLICATES SNL_INCLUDE_DIRS)
ENDIF(SNL_INCLUDE_DIRS)

IF(SNL_CFLAGS)
  LIST(REMOVE_DUPLICATES SNL_CFLAGS)
ENDIF(SNL_CFLAGS)

SET(SNL_CFLAGS ${SNL_CFLAGS}  CACHE STRING "Extra compiler flags required by SNL")

IF(SNL_FOUND)
  IF(NOT SNL_FIND_QUIETLY)
    MESSAGE(STATUS "Found SNL: headers at ${SNL_INCLUDE_DIRS}, libraries at ${SNL_LIBRARIES}")
  ENDIF(NOT SNL_FIND_QUIETLY)
ELSE(SNL_FOUND)
  IF(SNL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "SNL not found")
  ENDIF(SNL_FIND_REQUIRED)
ENDIF(SNL_FOUND)
