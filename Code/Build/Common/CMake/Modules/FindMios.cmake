# - Searches for an installation of the Mios communications library
#
# Defines:
#
#   Mios_FOUND         True if Mios was found, else false
#   Mios_LIBRARIES     Libraries to link
#   Mios_INCLUDE_DIRS  The directories containing the header files
#   Mios_CFLAGS        Extra compiler flags
#   Mios_LDFLAGS       Extra linker flags
#
# To specify an additional directory to search, set Mios_ROOT.
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(Mios_FOUND FALSE)
SET(Mios_CFLAGS)
SET(Mios_LDFLAGS)

# Look for the Mios header, first in the user-specified location and then in the system locations
SET(Mios_INCLUDE_DOC "The directory containing the Mios include file Mios/Mios.hpp")
FIND_PATH(Mios_INCLUDE_DIRS NAMES Mios/Mios.hpp PATHS ${Mios_ROOT} ${Mios_ROOT}/include ${Mios_ROOT}/Source
          DOC ${Mios_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT Mios_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(Mios_INCLUDE_DIRS NAMES Mios/Mios.hpp DOC ${Mios_INCLUDE_DOC})
ENDIF(NOT Mios_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(Mios_INCLUDE_DIRS)
  SET(Mios_LIBRARY_DIRS ${Mios_INCLUDE_DIRS})
  IF("${Mios_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mios_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(Mios_LIBRARY_DIRS ${Mios_LIBRARY_DIRS} PATH)
  ENDIF("${Mios_LIBRARY_DIRS}" MATCHES "/include$" OR "${Mios_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(Mios_DEBUG_LIBRARY
               NAMES Mios_d Miosd
               PATH_SUFFIXES "" Debug
               PATHS ${Mios_LIBRARY_DIRS} ${Mios_LIBRARY_DIRS}/lib ${Mios_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Mios_RELEASE_LIBRARY
               NAMES Mios
               PATH_SUFFIXES "" Release
               PATHS ${Mios_LIBRARY_DIRS} ${Mios_LIBRARY_DIRS}/lib ${Mios_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  SET(Mios_LIBRARIES)
  IF(Mios_DEBUG_LIBRARY AND Mios_RELEASE_LIBRARY)
    SET(Mios_LIBRARIES debug ${Mios_DEBUG_LIBRARY} optimized ${Mios_RELEASE_LIBRARY})
  ELSEIF(Mios_DEBUG_LIBRARY)
    SET(Mios_LIBRARIES ${Mios_DEBUG_LIBRARY})
  ELSEIF(Mios_RELEASE_LIBRARY)
    SET(Mios_LIBRARIES ${Mios_RELEASE_LIBRARY})
  ENDIF(Mios_DEBUG_LIBRARY AND Mios_RELEASE_LIBRARY)

  IF(Mios_LIBRARIES)
    SET(Mios_FOUND TRUE)

    # Flags for importing symbols from dynamically linked libraries
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      #SET(Mios_CFLAGS "-DMIOS_DLL") # using static libs on windows for now
    ELSE(WIN32)
      IF("${Mios_LIBRARIES}" MATCHES ".dylib$" OR "${Mios_LIBRARIES}" MATCHES ".so$")
        SET(Mios_CFLAGS "-DMIOS_DLL")
      ENDIF("${Mios_LIBRARIES}" MATCHES ".dylib$" OR "${Mios_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)
  ENDIF(Mios_LIBRARIES)
ENDIF(Mios_INCLUDE_DIRS)

SET(Mios_CFLAGS ${Mios_CFLAGS}  CACHE STRING "Extra compiler flags required by Mios")
SET(Mios_LDLAGS ${Mios_LDFLAGS} CACHE STRING "Extra linker flags required by Mios")

IF(Mios_FOUND)
  IF(NOT Mios_FIND_QUIETLY)
    MESSAGE(STATUS "Found Mios: headers at ${Mios_INCLUDE_DIRS}, libraries at ${Mios_LIBRARIES}")
  ENDIF(NOT Mios_FIND_QUIETLY)
ELSE(Mios_FOUND)
  IF(Mios_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Mios not found")
  ENDIF(Mios_FIND_REQUIRED)
ENDIF(Mios_FOUND)
