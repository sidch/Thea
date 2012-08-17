# Find the Offscreen Mesa (OSMesa) headers and libraries.
#
# Once done this will define
#  OSMesa_FOUND             - system has OSMesa
#  OSMesa_INCLUDE_DIRS      - the OSMesa include directory (containing GL/osmesa.h)
#  OSMesa_CFLAGS            - extra compiler flags
#  OSMesa_LIBRARIES         - link to these to use OSMesa
#  OSMesa_LIBRARY_DIRS      - the directories containing the libraries
#  OSMesa_LDFLAGS           - extra linker flags
#  OSMesa_GL_LIBRARIES      - an OpenGL library (libGL.so or equivalent), if found in the same folder, plus dependencies. This
#                             may be needed to build apps with OSMesa. If OSMesa is located via pkg-config, any libGL
#                             dependencies are automatically pulled into OSMesa_LIBRARIES and this variable is blank. Else, it
#                             is difficult to determine whether or not libGL is really needed (or is even compatible), so it is
#                             detected separately and returned via this variable so that the caller can decide. If Mesa is
#                             installed cleanly into its own folder, it is probably safe to always add this variable to the link
#                             path -- it is a noop if no library was found.
#
# For GLU, use the following flags:
#  OSMesa_GLU_FOUND         - system has GLU (probably) compatible with OSMesa
#  OSMesa_GLU_INCLUDE_DIRS  - the GLU include directory (containing GL/glu.h)
#  OSMesa_GLU_CFLAGS        - extra compiler flags
#  OSMesa_GLU_LIBRARIES     - link to these to use GLU
#  OSMesa_GLU_LIBRARY_DIRS  - the directories containing the libraries
#  OSMesa_GLU_LDFLAGS       - extra linker flags for GLU
#
# To specify an additional directory to search, set OSMesa_ROOT.
#
# Author: Siddhartha Chaudhuri, 2009. Distributed under the BSD license.
#

SET(OSMesa_INCLUDE_PATH_DESCRIPTION
    "directory containing the OSMesa include file. E.g /usr/local/Mesa-7.3/include/ or c:\\Mesa-7.3\\include")
SET(OSMesa_DIR_MESSAGE "Set the OSMesa_INCLUDE_DIR CMake cache entry to the ${OSMesa_INCLUDE_PATH_DESCRIPTION}")
SET(OSMesa_DIR_SEARCH ${OSMesa_ROOT} ${OSMesa_ROOT}/include)

# Assume we didn't find anything
SET(OSMesa_FOUND FALSE)
SET(OSMesa_INCLUDE_DIRS )
SET(OSMesa_CFLAGS "")
SET(OSMesa_LIBRARY_DIRS )
SET(OSMesa_LDFLAGS "")

SET(OSMesa_GLU_FOUND FALSE)
SET(OSMesa_GLU_INCLUDE_DIRS )
SET(OSMesa_GLU_CFLAGS "")
SET(OSMesa_GLU_LIBRARY_DIRS )
SET(OSMesa_GLU_LDFLAGS "")

#
# Look for an installation, first in the user-specified locations and then in the system locations
#

FIND_PATH(OSMesa_INCLUDE_DIR NAMES GL/osmesa.h PATHS ${OSMesa_DIR_SEARCH} DOC "The ${OSMesa_INCLUDE_PATH_DESCRIPTION}"
	  NO_DEFAULT_PATH)

IF(NOT OSMesa_INCLUDE_DIR)  # now look in system locations
  FIND_PATH(OSMesa_INCLUDE_DIR NAMES GL/osmesa.h DOC "The ${OSMesa_INCLUDE_PATH_DESCRIPTION}")
ENDIF(NOT OSMesa_INCLUDE_DIR)

# Now try to get the include and library paths
IF(OSMesa_INCLUDE_DIR AND EXISTS "${OSMesa_INCLUDE_DIR}")

  SET(OSMesa_INCLUDE_DIRS ${OSMesa_INCLUDE_DIR})
  SET(OSMesa_GLU_INCLUDE_DIRS ${OSMesa_INCLUDE_DIR})

  # Look for the OSMesa library path.
  SET(OSMesa_LIBRARY_DIR ${OSMesa_INCLUDE_DIR})

  IF(EXISTS ${OSMesa_LIBRARY_DIR}/lib)
    SET(OSMesa_LIBRARY_DIR ${OSMesa_LIBRARY_DIR}/lib)
  ELSE(EXISTS ${OSMesa_LIBRARY_DIR}/lib)
    IF("${OSMesa_LIBRARY_DIR}" MATCHES "/include$")
      # Strip off the trailing "/include" in the path.
      GET_FILENAME_COMPONENT(OSMesa_LIBRARY_DIR ${OSMesa_LIBRARY_DIR} PATH)
    ENDIF("${OSMesa_LIBRARY_DIR}" MATCHES "/include$")

    IF(EXISTS "${OSMesa_LIBRARY_DIR}/lib")
      SET (OSMesa_LIBRARY_DIR ${OSMesa_LIBRARY_DIR}/lib)
    ENDIF(EXISTS "${OSMesa_LIBRARY_DIR}/lib")
  ENDIF(EXISTS ${OSMesa_LIBRARY_DIR}/lib)

  IF(OSMesa_LIBRARY_DIR AND EXISTS "${OSMesa_LIBRARY_DIR}")
    SET(OSMesa_LIBRARY_DIRS ${OSMesa_LIBRARY_DIR})
    SET(OSMesa_GLU_LIBRARY_DIRS ${OSMesa_LIBRARY_DIR})
  ENDIF(OSMesa_LIBRARY_DIR AND EXISTS "${OSMesa_LIBRARY_DIR}")
ENDIF(OSMesa_INCLUDE_DIR AND EXISTS "${OSMesa_INCLUDE_DIR}")

# Look for the actual library files
IF(OSMesa_LIBRARY_DIRS)

  # Try to use pkg-config
  SET(PKG_CONFIG_FAILED TRUE)
  SET(OSMesa_PC_DIR ${OSMesa_LIBRARY_DIR}/pkgconfig)

  IF(EXISTS ${OSMesa_PC_DIR}/osmesa.pc)
    FIND_PACKAGE(PkgConfig)
    IF(PKG_CONFIG_FOUND)
      SET(ENV{PKG_CONFIG_PATH} ${OSMesa_LIBRARY_DIR}/pkgconfig)

      # First try to find the main library OSMesa
      PKG_CHECK_MODULES(OSMesa_PC_OSMesa osmesa)
      IF(OSMesa_PC_OSMesa_FOUND)
	SET(PKG_CONFIG_FAILED FALSE)

	SET(OSMesa_LIBRARIES ${OSMesa_PC_OSMesa_LIBRARIES})
	SET(OSMesa_LIBRARY_DIRS ${OSMesa_PC_OSMesa_LIBRARY_DIRS})
	SET(OSMesa_INCLUDE_DIRS ${OSMesa_PC_OSMesa_INCLUDE_DIRS})
	SET(OSMesa_LDFLAGS ${OSMesa_PC_OSMesa_LDFLAGS_OTHER})
	SET(OSMesa_CFLAGS ${OSMesa_PC_OSMesa_CFLAGS_OTHER})
	SET(OSMesa_GL_LIBRARIES )  # libGL is pulled in via pkg-config into OSMesa_LIBRARIES if it is required

	# Finally try to find libGLU
	IF(EXISTS ${OSMesa_PC_DIR}/glu.pc)  # only a libGLU from the same installation is valid
	  PKG_CHECK_MODULES(OSMesa_PC_GLU glu)
	ENDIF(EXISTS ${OSMesa_PC_DIR}/glu.pc)

	IF(OSMesa_PC_GLU_FOUND)
	  SET(OSMesa_GLU_LIBRARIES ${OSMesa_PC_GLU_LIBRARIES})
	  SET(OSMesa_GLU_LIBRARY_DIRS ${OSMesa_PC_GLU_LIBRARY_DIRS})
	  SET(OSMesa_GLU_INCLUDE_DIRS ${OSMesa_PC_GLU_INCLUDE_DIRS})
	  SET(OSMesa_GLU_LDFLAGS ${OSMesa_PC_GLU_LDFLAGS_OTHER})
	  SET(OSMesa_GLU_CFLAGS ${OSMesa_PC_GLU_CFLAGS_OTHER})
	ELSE(OSMesa_PC_GLU_FOUND)
	  SET(OSMesa_GLU_LIBRARIES )
	ENDIF(OSMesa_PC_GLU_FOUND)

      ENDIF(OSMesa_PC_OSMesa_FOUND)
    ENDIF(PKG_CONFIG_FOUND)
  ENDIF(EXISTS ${OSMesa_PC_DIR}/osmesa.pc)

  # If we couldn't use pkg-config, try to locate the libs manually
  IF(PKG_CONFIG_FAILED)
    FIND_LIBRARY(OSMesa_LIBRARIES NAMES OSMesa osmesa32 PATHS ${OSMesa_LIBRARY_DIRS} NO_DEFAULT_PATH)
    FIND_LIBRARY(OSMesa_GL_LIBRARIES NAMES GL opengl opengl32 PATHS ${OSMesa_LIBRARY_DIRS} NO_DEFAULT_PATH)
    FIND_LIBRARY(OSMesa_GLU_LIBRARIES NAMES GLU glu32 PATHS ${OSMesa_LIBRARY_DIRS} NO_DEFAULT_PATH)
  ENDIF(PKG_CONFIG_FAILED)
ENDIF(OSMesa_LIBRARY_DIRS)

# If libGL was not found, set the variable to null so that it can be safely passed to the linker as a noop
IF(NOT OSMesa_GL_LIBRARIES)
  SET(OSMesa_GL_LIBRARIES "")
ENDIF(NOT OSMesa_GL_LIBRARIES)

# If the libraries were found, set the appropriate flags
IF(OSMesa_LIBRARIES)
  SET(OSMesa_FOUND TRUE)

  # Try to fix a bug in Linux that omits libdl
  IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    SET(OSMesa_LIBRARIES ${OSMesa_LIBRARIES} -ldl)
  ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")

  IF(OSMesa_GLU_LIBRARIES)
    SET(OSMesa_GLU_FOUND TRUE)
  ENDIF(OSMesa_GLU_LIBRARIES)
ENDIF(OSMesa_LIBRARIES)

IF(OSMesa_FOUND)
  IF(NOT OSMesa_FIND_QUIETLY)
    MESSAGE(STATUS "Found OSMesa: headers at ${OSMesa_INCLUDE_DIRS}, libraries at ${OSMesa_LIBRARY_DIRS}")

    IF(NOT OSMesa_GLU_FOUND)
      MESSAGE(STATUS "OSMesa-compatible GLU not found")
    ENDIF(NOT OSMesa_GLU_FOUND)
  ENDIF(NOT OSMesa_FIND_QUIETLY)
ELSE(OSMesa_FOUND)
  IF(OSMesa_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "OSMesa was not found. ${OSMesa_DIR_MESSAGE}")
  ELSE(OSMesa_FIND_REQUIRED)
    IF(NOT OSMesa_FIND_QUIETLY)
      MESSAGE(STATUS "OSMesa was not found. ${OSMesa_DIR_MESSAGE}")
    ENDIF(NOT OSMesa_FIND_QUIETLY)
  ENDIF(OSMesa_FIND_REQUIRED)
ENDIF(OSMesa_FOUND)
