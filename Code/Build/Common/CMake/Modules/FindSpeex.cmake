# - Try to find the Speex library
# Once done this will define
#
#  SPEEX_FOUND           - system has Speex
#  SPEEX_INCLUDE_DIR     - the Speex include directory
#  SPEEX_LIBRARY         - The Speex library
#
# To specify an additional directory to search for Speex, set Speex_ROOT.
#
# Copyright (c) 2006, Richard Laerkaeng, <richard@goteborg.utfors.se>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
# (Copied and slightly modified from the KDE project)
#

INCLUDE(CheckLibraryExists)

FIND_PATH(SPEEX_INCLUDE_DIR speex/speex.h PATHS ${Speex_ROOT} ${Speex_ROOT}/include NO_DEFAULT_PATH)
IF(NOT SPEEX_INCLUDE_DIR)
  FIND_PATH(SPEEX_INCLUDE_DIR speex/speex.h)
ENDIF(NOT SPEEX_INCLUDE_DIR)


FIND_LIBRARY(SPEEX_DEBUG_LIBRARY NAMES speexd speex_d libspeexd libspeex_d PATH_SUFFIXES "" Debug
             PATHS ${Speex_ROOT} ${Speex_ROOT}/lib NO_DEFAULT_PATH)
IF(NOT SPEEX_DEBUG_LIBRARY)  # now look in system locations
  FIND_LIBRARY(SPEEX_DEBUG_LIBRARY NAMES speexd speex_d libspeexd libspeex_d PATH_SUFFIXES "" Debug)
ENDIF(NOT SPEEX_DEBUG_LIBRARY)


FIND_LIBRARY(SPEEX_RELEASE_LIBRARY NAMES speex libspeex PATH_SUFFIXES "" Release
             PATHS ${Speex_ROOT} ${Speex_ROOT}/lib NO_DEFAULT_PATH)
IF(NOT SPEEX_RELEASE_LIBRARY)  # now look in system locations
  FIND_LIBRARY(SPEEX_RELEASE_LIBRARY NAMES speex libspeex PATH_SUFFIXES "" Release)
ENDIF(NOT SPEEX_RELEASE_LIBRARY)


SET(SPEEX_LIBRARY)
IF(SPEEX_DEBUG_LIBRARY AND SPEEX_RELEASE_LIBRARY)
  SET(SPEEX_LIBRARY debug ${SPEEX_DEBUG_LIBRARY} optimized ${SPEEX_RELEASE_LIBRARY})
ELSEIF(SPEEX_DEBUG_LIBRARY)
  SET(SPEEX_LIBRARY ${SPEEX_DEBUG_LIBRARY})
ELSEIF(SPEEX_RELEASE_LIBRARY)
  SET(SPEEX_LIBRARY ${SPEEX_RELEASE_LIBRARY})
ENDIF(SPEEX_DEBUG_LIBRARY AND SPEEX_RELEASE_LIBRARY)


IF(SPEEX_INCLUDE_DIR AND SPEEX_LIBRARY)
   SET(SPEEX_FOUND TRUE)
ELSE(SPEEX_INCLUDE_DIR AND SPEEX_LIBRARY)
   SET(SPEEX_FOUND FALSE)
ENDIF(SPEEX_INCLUDE_DIR AND SPEEX_LIBRARY)


IF(SPEEX_FOUND)
   IF(NOT Speex_FIND_QUIETLY)
      MESSAGE(STATUS "Found Speex: headers at ${SPEEX_INCLUDE_DIR}, libraries at ${SPEEX_LIBRARY}")
   ENDIF(NOT Speex_FIND_QUIETLY)
ELSE(SPEEX_FOUND)
   IF(Speex_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Speex library")
   ELSE(Speex_FIND_REQUIRED)
      MESSAGE(STATUS_ERROR "Could not find Speex library")
   ENDIF(Speex_FIND_REQUIRED)
ENDIF(SPEEX_FOUND)
