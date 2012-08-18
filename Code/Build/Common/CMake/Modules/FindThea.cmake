# - Searches for an installation of the Thea library
#
# Defines:
#
#   Thea_FOUND           True if Thea was found, else false
#   Thea_LIBRARIES       Libraries to link
#   Thea_LIBRARY_DIRS    Additional directories for libraries. These do not necessarily correspond to Thea_LIBRARIES, and both
#                        variables must be passed to the linker.
#   Thea_INCLUDE_DIRS    The directories containing the header files
#   Thea_CFLAGS          Extra compiler flags
#   Thea_DEBUG_CFLAGS    Extra compiler flags to be used only in debug builds
#   Thea_RELEASE_CFLAGS  Extra compiler flags to be used only in release builds
#   Thea_LDFLAGS         Extra linker flags
#
# To specify an additional directory to search, set Thea_ROOT.
# To prevent automatically searching for all dependencies, set Thea_NO_DEPENDENCIES to true.
#
# Author: Siddhartha Chaudhuri, 2009
#
# Revisions:
#   - 2011-04-12: Locate dependencies automatically, without requiring the caller to do so separately. [SC]
#

SET(Thea_FOUND FALSE)
SET(Thea_LIBRARY_DIRS )
SET(Thea_CFLAGS )
SET(Thea_LDFLAGS )

# Look for the Thea header, first in the user-specified location and then in the system locations
SET(Thea_INCLUDE_DOC "The directory containing the Thea include file Thea/Thea.hpp")
FIND_PATH(Thea_INCLUDE_DIRS NAMES Thea/Common.hpp PATHS ${Thea_ROOT} ${Thea_ROOT}/include ${Thea_ROOT}/Source
          DOC ${Thea_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT Thea_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(Thea_INCLUDE_DIRS NAMES Thea/Common.hpp DOC ${Thea_INCLUDE_DOC})
ENDIF(NOT Thea_INCLUDE_DIRS)

# Only look for the library file in the immediate neighbourhood of the include directory
IF(Thea_INCLUDE_DIRS)
  SET(Thea_LIBRARY_DIRS ${Thea_INCLUDE_DIRS})
  IF("${Thea_LIBRARY_DIRS}" MATCHES "/include$" OR "${Thea_LIBRARY_DIRS}" MATCHES "/Source$")
    # Strip off the trailing "/include" or "/Source" from the path
    GET_FILENAME_COMPONENT(Thea_LIBRARY_DIRS ${Thea_LIBRARY_DIRS} PATH)
  ENDIF("${Thea_LIBRARY_DIRS}" MATCHES "/include$" OR "${Thea_LIBRARY_DIRS}" MATCHES "/Source$")

  FIND_LIBRARY(Thea_DEBUG_LIBRARY
               NAMES Thea_d Thead
               PATH_SUFFIXES "" Debug
               PATHS ${Thea_LIBRARY_DIRS} ${Thea_LIBRARY_DIRS}/lib ${Thea_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  FIND_LIBRARY(Thea_RELEASE_LIBRARY
               NAMES Thea
               PATH_SUFFIXES "" Release
               PATHS ${Thea_LIBRARY_DIRS} ${Thea_LIBRARY_DIRS}/lib ${Thea_LIBRARY_DIRS}/Build/lib NO_DEFAULT_PATH)

  SET(Thea_LIBRARIES)
  IF(Thea_DEBUG_LIBRARY AND Thea_RELEASE_LIBRARY)
    SET(Thea_LIBRARIES debug ${Thea_DEBUG_LIBRARY} optimized ${Thea_RELEASE_LIBRARY})
  ELSEIF(Thea_DEBUG_LIBRARY)
    SET(Thea_LIBRARIES ${Thea_DEBUG_LIBRARY})
  ELSEIF(Thea_RELEASE_LIBRARY)
    SET(Thea_LIBRARIES ${Thea_RELEASE_LIBRARY})
  ENDIF(Thea_DEBUG_LIBRARY AND Thea_RELEASE_LIBRARY)

  IF(Thea_LIBRARIES)
    SET(Thea_FOUND TRUE)

    # Flags for importing symbols from dynamically linked libraries
    IF(WIN32)
      # What's a good way of testing whether the .lib is static, or merely exports symbols from a DLL? For now, let's assume
      # it always exports (or hope that __declspec(dllimport) is a noop for static libraries)
      SET(Thea_CFLAGS "-DTHEA_DLL -DTHEA_DLL_IMPORTS")
    ELSE(WIN32)
      IF("${Thea_LIBRARIES}" MATCHES ".dylib$" OR "${Thea_LIBRARIES}" MATCHES ".so$")
        SET(Thea_CFLAGS "-DTHEA_DLL -DTHEA_DLL_IMPORTS")
      ENDIF("${Thea_LIBRARIES}" MATCHES ".dylib$" OR "${Thea_LIBRARIES}" MATCHES ".so$")
    ENDIF(WIN32)

    # Read extra flags to be used to build Thea
    SET(Thea_BUILD_FLAGS_FILE "${Thea_INCLUDE_DIRS}/Thea/BuildFlags.txt")
    IF(EXISTS "${Thea_BUILD_FLAGS_FILE}")
      FILE(READ "${Thea_BUILD_FLAGS_FILE}" Thea_BUILD_FLAGS)
      STRING(REGEX REPLACE "\n" " " Thea_BUILD_FLAGS "${Thea_BUILD_FLAGS}")
      SET(Thea_CFLAGS "${Thea_CFLAGS} ${Thea_BUILD_FLAGS}")
    ENDIF(EXISTS "${Thea_BUILD_FLAGS_FILE}")

  ENDIF(Thea_LIBRARIES)
ENDIF(Thea_INCLUDE_DIRS)

IF(NOT Thea_NO_DEPENDENCIES)

# Dependency: G3D
IF(Thea_FOUND)
  IF(EXISTS ${Thea_ROOT}/installed-g3d)
    SET(G3D_ROOT ${Thea_ROOT}/installed-g3d)
  ELSE()
    SET(G3D_ROOT ${Thea_ROOT})
  ENDIF()
  FIND_PACKAGE(G3D)
  IF(G3D_FOUND)
    SET(Thea_LIBRARIES ${Thea_LIBRARIES} ${G3D_LIBRARIES})
    SET(Thea_INCLUDE_DIRS ${Thea_INCLUDE_DIRS} ${G3D_INCLUDE_DIRS})
  ELSE(G3D_FOUND)
    MESSAGE(STATUS "Thea: G3D not found")
    SET(Thea_FOUND FALSE)
  ENDIF(G3D_FOUND)
ENDIF(Thea_FOUND)

# Dependency: Boost
IF(Thea_FOUND)
  SET(Boost_USE_STATIC_LIBS      ON)
  SET(Boost_USE_MULTITHREADED    ON)
  SET(Boost_USE_STATIC_RUNTIME  OFF)
  INCLUDE(BoostAdditionalVersions)
  IF(EXISTS ${Thea_ROOT}/installed-boost)
    SET(BOOST_ROOT ${Thea_ROOT}/installed-boost)
  ELSE()
    SET(BOOST_ROOT ${Thea_ROOT})
  ENDIF()
  FIND_PACKAGE(Boost)
  IF(Boost_FOUND)
    SET(Thea_INCLUDE_DIRS ${Thea_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    # We'll add the libraries below, after CGAL, which depends on Boost
  ELSE(Boost_FOUND)
    MESSAGE(STATUS "Thea: Boost not found")
    SET(Thea_FOUND FALSE)
  ENDIF(Boost_FOUND)
ENDIF(Thea_FOUND)

# Dependency: CGAL
IF(Thea_FOUND)
  IF(EXISTS ${Thea_ROOT}/installed-cgal)
    SET(CGAL_ROOT ${Thea_ROOT}/installed-cgal)
  ELSE()
    SET(CGAL_ROOT ${Thea_ROOT})
  ENDIF()
  FIND_PACKAGE(CGAL)
  IF(CGAL_FOUND)
    SET(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS} ${CGAL_3RD_PARTY_INCLUDE_DIRS})
    SET(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES})
    IF(NOT CGAL_LIBRARY)
      MESSAGE(STATUS "CGAL libraries will be auto-linked")
      SET(Thea_LIBRARY_DIRS ${Thea_LIBRARY_DIRS} ${CGAL_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
    ENDIF()
    # -O3 might cause problems with old versions of gcc 4 on OS X
    LIST(REMOVE_ITEM CGAL_RELEASE_CFLAGS "-O3")

    SET(Thea_LIBRARIES ${Thea_LIBRARIES} ${CGAL_LIBRARIES})
    SET(Thea_INCLUDE_DIRS ${Thea_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS})
    SET(Thea_CFLAGS ${Thea_CFLAGS} ${CGAL_3RD_PARTY_DEFINITIONS})
    SET(Thea_DEBUG_CFLAGS ${Thea_DEBUG_CFLAGS} ${CGAL_DEBUG_CFLAGS})
    SET(Thea_RELEASE_CFLAGS ${Thea_RELEASE_CFLAGS} ${CGAL_RELEASE_CFLAGS})
  ELSE(CGAL_FOUND)
    MESSAGE(STATUS "Thea: CGAL not found")
    SET(Thea_FOUND FALSE)
  ENDIF(CGAL_FOUND)
ENDIF(Thea_FOUND)

# CGAL depends on Boost but the block above also needs access to the Boost library dirs, so we add the libraries here, in the
# correct sequence
IF(Thea_FOUND)
  SET(Thea_LIBRARIES ${Thea_LIBRARIES} ${Boost_LIBRARIES})
ENDIF(Thea_FOUND)

# Dependency: Lib3ds
IF(Thea_FOUND)
  IF(EXISTS ${Thea_ROOT}/installed-lib3ds)
    SET(Lib3ds_ROOT ${Thea_ROOT}/installed-lib3ds)
  ELSE()
    SET(Lib3ds_ROOT ${Thea_ROOT})
  ENDIF()
  FIND_PACKAGE(Lib3ds)
  IF(Lib3ds_FOUND)
    SET(Thea_LIBRARIES ${Thea_LIBRARIES} ${Lib3ds_LIBRARIES})
    SET(Thea_INCLUDE_DIRS ${Thea_INCLUDE_DIRS} ${Lib3ds_INCLUDE_DIRS})
    SET(Thea_CFLAGS ${Thea_CFLAGS} -DTHEA_LIB3DS_VERSION_MAJOR=${Lib3ds_VERSION_MAJOR})
  ELSE(Lib3ds_FOUND)
    MESSAGE(STATUS "Thea: lib3ds not found")
    SET(Thea_FOUND FALSE)
  ENDIF(Lib3ds_FOUND)
ENDIF(Thea_FOUND)

# Dependency: CLUTO
IF(Thea_FOUND)
  IF(EXISTS ${Thea_ROOT}/installed-cluto)
    SET(CLUTO_ROOT ${Thea_ROOT}/installed-cluto)
  ELSE()
    SET(CLUTO_ROOT ${Thea_ROOT})
  ENDIF()
  FIND_PACKAGE(CLUTO)
  IF(CLUTO_FOUND)
    SET(Thea_LIBRARIES ${Thea_LIBRARIES} ${CLUTO_LIBRARIES})
    SET(Thea_INCLUDE_DIRS ${Thea_INCLUDE_DIRS} ${CLUTO_INCLUDE_DIRS})
    SET(Thea_CFLAGS ${Thea_CFLAGS} -DTHEA_ENABLE_CLUTO)
  ELSE(CLUTO_FOUND)
    MESSAGE(STATUS "Thea: CLUTO not found")  # this is not a fatal error
  ENDIF(CLUTO_FOUND)
ENDIF(Thea_FOUND)

ENDIF(NOT Thea_NO_DEPENDENCIES)

# Remove duplicate entries from lists, else the same dirs and flags can repeat many times

# Don't remove duplicates from Thea_LIBRARIES -- the list includes repetitions of "debug" and "optimized"

IF(Thea_LIBRARY_DIRS)
  LIST(REMOVE_DUPLICATES Thea_LIBRARY_DIRS)
ENDIF(Thea_LIBRARY_DIRS)

IF(Thea_INCLUDE_DIRS)
  LIST(REMOVE_DUPLICATES Thea_INCLUDE_DIRS)
ENDIF(Thea_INCLUDE_DIRS)

IF(Thea_CFLAGS)
  LIST(REMOVE_DUPLICATES Thea_CFLAGS)
ENDIF(Thea_CFLAGS)

IF(Thea_DEBUG_CFLAGS)
  LIST(REMOVE_DUPLICATES Thea_DEBUG_CFLAGS)
ENDIF(Thea_DEBUG_CFLAGS)

IF(Thea_RELEASE_CFLAGS)
  LIST(REMOVE_DUPLICATES Thea_RELEASE_CFLAGS)
ENDIF(Thea_RELEASE_CFLAGS)

IF(Thea_LDFLAGS)
  LIST(REMOVE_DUPLICATES Thea_LDFLAGS)
ENDIF(Thea_LDFLAGS)

SET(Thea_LIBRARY_DIRS ${Thea_LIBRARY_DIRS} CACHE STRING "Additional directories for libraries required by Thea")
SET(Thea_CFLAGS ${Thea_CFLAGS}  CACHE STRING "Extra compiler flags required by Thea")
SET(Thea_LDFLAGS ${Thea_LDFLAGS} CACHE STRING "Extra linker flags required by Thea")

IF(Thea_FOUND)
  IF(NOT Thea_FIND_QUIETLY)
    MESSAGE(STATUS "Found Thea: headers at ${Thea_INCLUDE_DIRS}, libraries at ${Thea_LIBRARIES}")
  ENDIF(NOT Thea_FIND_QUIETLY)
ELSE(Thea_FOUND)
  IF(Thea_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Thea not found")
  ENDIF(Thea_FIND_REQUIRED)
ENDIF(Thea_FOUND)
