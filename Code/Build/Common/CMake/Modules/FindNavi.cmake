# Searches for a Navi installation
#
# Defines:
#
#   Navi_FOUND         True if Navi was found, else false
#   Navi_LIBRARIES     Libraries to link
#   Navi_INCLUDE_DIRS  The directories containing the header files
#   Navi_CFLAGS        Additional compiler flags
#
# Currently, this code requires you to manually set the root directory of the Navi installation, as the value of the variable
# Navi_ROOT.
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(Navi_FOUND FALSE)

IF(Navi_ROOT AND EXISTS ${Navi_ROOT})
  IF(WIN32)  # Windows

    SET(Navi_INCLUDE_DIRS ${Navi_ROOT}/Dependencies/win32/llmozlib/include ${Navi_ROOT}/Include)
    SET(Navi_LIBRARIES
        debug ${Navi_ROOT}/Lib/Navi_d_DLL.lib
        debug ${Navi_ROOT}/Dependencies/win32/llmozlib/lib/llmozlib_d.lib
        optimized ${Navi_ROOT}/Lib/Navi_DLL.lib
        optimized ${Navi_ROOT}/Dependencies/win32/llmozlib/lib/llmozlib.lib)
    SET(Navi_CFLAGS)
    SET(Navi_FOUND TRUE)

  ELSEIF(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")  # OS X

    SET(Navi_INCLUDE_DIRS
        ${Navi_ROOT}/Dependencies/win32/llmozlib/include
        ${Navi_ROOT}/Dependencies/all/utf8-cpp
        ${Navi_ROOT}/Include)
    FIND_LIBRARY(Navi_LIBRARIES NAMES Navi PATHS ${Navi_ROOT}/Source NO_DEFAULT_PATH)
    SET(Navi_CFLAGS)

    IF(Navi_LIBRARIES)
      SET(Navi_FOUND TRUE)
    ENDIF(Navi_LIBRARIES)

  ELSE(WIN32)  # Linux etc

    SET(Navi_INCLUDE_DIRS
        ${Navi_ROOT}/NaviDemo/Headers
        ${Navi_ROOT}/Dependencies/win32/llmozlib/include
        ${Navi_ROOT}/Dependencies/all/utf8-cpp
        ${Navi_ROOT}/Include)
    FIND_LIBRARY(Navi_LIBRARIES NAMES Navi PATHS ${Navi_ROOT}/Dependencies/linux/alt/usr/lib NO_DEFAULT_PATH)

    FIND_PACKAGE(PkgConfig)
    IF(NOT PKG_CONFIG_FOUND)
      MESSAGE("Could not find pkg-config (to get Navi build flags)")
    ELSE(NOT PKG_CONFIG_FOUND)
      SET(PKG_CONFIG_PATH_BACKUP $ENV{PKG_CONFIG_PATH})  # save current PKG_CONFIG_PATH
      SET(ENV{PKG_CONFIG_PATH} ${Navi_ROOT}/Dependencies/linux/alt/usr/lib/pkgconfig)

      EXEC_PROGRAM(${PKG_CONFIG_EXECUTABLE}
                   ARGS --cflags OIS OGRE cairo webkit-1.0 libcurl
                   OUTPUT_VARIABLE Navi_CFLAGS)

      EXEC_PROGRAM(${PKG_CONFIG_EXECUTABLE}
                   ARGS --libs --static webkit-1.0 libcurl
                   OUTPUT_VARIABLE Navi_DEPS_1)
      EXEC_PROGRAM(${PKG_CONFIG_EXECUTABLE}
                   ARGS --libs gthread-2.0 OIS OGRE cairo sqlite3 libxslt
                   OUTPUT_VARIABLE Navi_DEPS_2)

      SET(ENV{PKG_CONFIG_PATH} ${PKG_CONFIG_PATH_BACKUP})  # restore saved PKG_CONFIG_PATH

      # The WebKit static lib installation is screwed up because the pkg-config file does not have the full list
      # of dependencies. We'll correct this by parsing the .la file created by Libtool.
      SET(WebKit_LIBTOOL_ARCHIVE ${Navi_ROOT}/Dependencies/linux/alt/usr/lib/libwebkit-1.0.la)
      IF(EXISTS ${WebKit_LIBTOOL_ARCHIVE})
        INCLUDE(GetLibtoolizedDependencies)
        GET_LIBTOOLIZED_DEPENDENCIES(${WebKit_LIBTOOL_ARCHIVE} WebKit_MAIN_ARCHIVE WebKit_DEPS)
      ENDIF(EXISTS ${WebKit_LIBTOOL_ARCHIVE})

      SET(Navi_DEPS ${Navi_DEPS_1} ${WebKit_DEPS} ${Navi_DEPS_2} "-ldl")
      STRING(REGEX REPLACE ";" " " Navi_DEPS "${Navi_DEPS}")  # make it all one long space-separated line
      SET(Navi_LIBRARIES ${Navi_LIBRARIES} ${Navi_DEPS})

      SET(Navi_FOUND TRUE)
    ENDIF(NOT PKG_CONFIG_FOUND)

  ENDIF(WIN32)
ENDIF(Navi_ROOT AND EXISTS ${Navi_ROOT})

IF(Navi_FOUND)
  IF(NOT Navi_FIND_QUIETLY)
    MESSAGE(STATUS "Found Navi: headers at ${Navi_INCLUDE_DIRS}, libraries at ${Navi_LIBRARIES}")
  ENDIF(NOT Navi_FIND_QUIETLY)
ELSE(Navi_FOUND)
  IF(Navi_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Navi not found")
  ENDIF(Navi_FIND_REQUIRED)
ENDIF(Navi_FOUND)
