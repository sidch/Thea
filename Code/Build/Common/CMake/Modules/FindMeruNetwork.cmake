# - Searches for an installation of the Meru networking library
#
# Defines:
#
#   MeruNetwork_FOUND         True if MeruNetwork was found, else false
#   MeruNetwork_LIBRARIES     Libraries to link
#   MeruNetwork_INCLUDE_DIRS  The directories containing the header files
#   MeruNetwork_CFLAGS        Extra compiler flags
#   MeruNetwork_LDFLAGS       Extra linker flags
#
# To specify an additional directory to search, set MeruNetwork_ROOT.
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(MeruNetwork_FOUND FALSE)
SET(MeruNetwork_CFLAGS)
SET(MeruNetwork_LDFLAGS)

# Look for the MeruNetwork header, first in the user-specified location and then in the system locations
SET(MeruNetwork_INCLUDE_DOC "The directory containing the MeruNetwork include file Meru/Network/Networking.hpp")
FIND_PATH(MeruNetwork_INCLUDE_DIRS NAMES Meru/Network/Networking.hpp
          PATHS ${MeruNetwork_ROOT} ${MeruNetwork_ROOT}/include ${MeruNetwork_ROOT}/Source
          DOC ${MeruNetwork_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT MeruNetwork_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(MeruNetwork_INCLUDE_DIRS NAMES Meru/Network/Networking.hpp DOC ${MeruNetwork_INCLUDE_DOC})
ENDIF(NOT MeruNetwork_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(MeruNetwork_INCLUDE_DIRS)
  SET(MeruNetwork_LIBRARY_DIRS ${MeruNetwork_INCLUDE_DIRS})
  IF("${MeruNetwork_LIBRARY_DIRS}" MATCHES "/include$" OR "${MeruNetwork_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(MeruNetwork_LIBRARY_DIRS ${MeruNetwork_LIBRARY_DIRS} PATH)
  ENDIF("${MeruNetwork_LIBRARY_DIRS}" MATCHES "/include$" OR "${MeruNetwork_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(MeruNetwork_DEBUG_LIBRARY
               NAMES MeruNetwork_d MeruNetworkd
               PATH_SUFFIXES "" Debug
               PATHS ${MeruNetwork_LIBRARY_DIRS} ${MeruNetwork_LIBRARY_DIRS}/lib ${MeruNetwork_LIBRARY_DIRS}/Build/lib
               NO_DEFAULT_PATH)

  FIND_LIBRARY(MeruNetwork_RELEASE_LIBRARY
               NAMES MeruNetwork
               PATH_SUFFIXES "" Release
               PATHS ${MeruNetwork_LIBRARY_DIRS} ${MeruNetwork_LIBRARY_DIRS}/lib ${MeruNetwork_LIBRARY_DIRS}/Build/lib
               NO_DEFAULT_PATH)

  SET(MeruNetwork_LIBRARIES)
  IF(MeruNetwork_DEBUG_LIBRARY AND MeruNetwork_RELEASE_LIBRARY)
    SET(MeruNetwork_LIBRARIES debug ${MeruNetwork_DEBUG_LIBRARY} optimized ${MeruNetwork_RELEASE_LIBRARY})
  ELSEIF(MeruNetwork_DEBUG_LIBRARY)
    SET(MeruNetwork_LIBRARIES ${MeruNetwork_DEBUG_LIBRARY})
  ELSEIF(MeruNetwork_RELEASE_LIBRARY)
    SET(MeruNetwork_LIBRARIES ${MeruNetwork_RELEASE_LIBRARY})
  ENDIF(MeruNetwork_DEBUG_LIBRARY AND MeruNetwork_RELEASE_LIBRARY)

  IF(MeruNetwork_LIBRARIES)
    SET(MeruNetwork_FOUND TRUE)

    # Flags for importing symbols from dynamically linked libraries
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      # SET(MeruNetwork_CFLAGS "-DMERU_NETWORK_DLL") # using static libs on windows for now
    ELSE(WIN32)
      IF("${MeruNetwork_LIBRARIES}" MATCHES ".dylib$" OR "${MeruNetwork_LIBRARIES}" MATCHES ".so$")
        SET(MeruNetwork_CFLAGS "-DMERU_NETWORK_DLL")
      ENDIF("${MeruNetwork_LIBRARIES}" MATCHES ".dylib$" OR "${MeruNetwork_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)
  ENDIF(MeruNetwork_LIBRARIES)
ENDIF(MeruNetwork_INCLUDE_DIRS)

SET(MeruNetwork_CFLAGS ${MeruNetwork_CFLAGS}  CACHE STRING "Extra compiler flags required by MeruNetwork")
SET(MeruNetwork_LDLAGS ${MeruNetwork_LDFLAGS} CACHE STRING "Extra linker flags required by MeruNetwork")

IF(MeruNetwork_FOUND)
  IF(NOT MeruNetwork_FIND_QUIETLY)
    MESSAGE(STATUS "Found MeruNetwork: headers at ${MeruNetwork_INCLUDE_DIRS}, libraries at ${MeruNetwork_LIBRARIES}")
  ENDIF(NOT MeruNetwork_FIND_QUIETLY)
ELSE(MeruNetwork_FOUND)
  IF(MeruNetwork_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "MeruNetwork not found")
  ENDIF(MeruNetwork_FIND_REQUIRED)
ENDIF(MeruNetwork_FOUND)
