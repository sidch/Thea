# - Searches for an installation of the OpenAL library.
#
# This replaces the version of FindOpenAL.cmake supplied with CMake to comply with Meru requirements.
#
# Defines:
#
#   OPENAL_FOUND          True if OpenAL was found, else false
#   OPENAL_LIBRARY        The location of the OpenAL library (might include support libraries)
#   OPENAL_INCLUDE_DIR    The directory containing the OpenAL header file (al.h)
#   ALUT_FOUND            True if ALUT was found, else false
#   ALUT_LIBRARY          The location of the ALUT library (might include support libraries)
#   ALUT_INCLUDE_DIR      The directory containing the ALUT header file (alut.h)
#
# To specify an additional directory to search for OpenAL, set OPENAL_ROOT.
# To specify an additional directory to search for ALUT, set ALUT_ROOT.
# To skip checking for ALUT, set ALUT_SKIP to true.
#
# On OS X, to avoid linking to a framework, set OPENAL_NO_FRAMEWORK to true ("cmake -DOPENAL_NO_FRAMEWORK=TRUE ."):
#
# Author: Siddhartha Chaudhuri, 2008
#

SET(OPENAL_FOUND FALSE)

# Only look for a framework in OS X if the flags are all false
IF(NOT OPENAL_NO_FRAMEWORK)

  IF(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")  # OS X

    INCLUDE(CMakeFindFrameworks)
    CMAKE_FIND_FRAMEWORKS(OpenAL)
    IF(OpenAL_FRAMEWORKS)
      LIST(GET OpenAL_FRAMEWORKS 0 OPENAL_LIBRARY)
      SET(OPENAL_INCLUDE_DIR ${OPENAL_LIBRARY}/Headers)
      SET(OPENAL_FOUND TRUE)
    ENDIF(OpenAL_FRAMEWORKS)

  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")

ENDIF(NOT OPENAL_NO_FRAMEWORK)


# Look for a non-framework install
IF(NOT OPENAL_FOUND)

  SET(OPENAL_INCLUDE_DOC "The directory containing the OpenAL include file al.h")
  FIND_PATH(OPENAL_INCLUDE_DIR NAMES al.h PATH_SUFFIXES "" "AL" "OpenAL"
            PATHS ${OPENAL_ROOT} ${OPENAL_ROOT}/include DOC ${OPENAL_INCLUDE_DOC} NO_DEFAULT_PATH)
  IF(NOT OPENAL_INCLUDE_DIR)  # now look in system locations
    FIND_PATH(OPENAL_INCLUDE_DIR NAMES al.h PATH_SUFFIXES "" "AL" "OpenAL" DOC ${OPENAL_INCLUDE_DOC})
  ENDIF(NOT OPENAL_INCLUDE_DIR)

  INCLUDE(CheckFunctionExists)
  MACRO(OPENAL_CHECK_LIBRARY variable)  # checks whether we have located all the required libraries
    SET(${variable})

    SET(TMP_CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS})
    SET(TMP_CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS})
    SET(TMP_CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES})
    SET(TMP_CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES})

    SET(CMAKE_REQUIRED_FLAGS)
    SET(CMAKE_REQUIRED_DEFINITIONS)
    SET(CMAKE_REQUIRED_INCLUDES)
    SET(CMAKE_REQUIRED_LIBRARIES ${OPENAL_LIBRARY})

    CHECK_FUNCTION_EXISTS(alcGetCurrentContext ${variable})

    SET(CMAKE_REQUIRED_FLAGS        ${TMP_CMAKE_REQUIRED_FLAGS})
    SET(CMAKE_REQUIRED_DEFINITIONS  ${TMP_CMAKE_REQUIRED_DEFINITIONS})
    SET(CMAKE_REQUIRED_INCLUDES     ${TMP_CMAKE_REQUIRED_INCLUDES})
    SET(CMAKE_REQUIRED_LIBRARIES    ${TMP_CMAKE_REQUIRED_LIBRARIES})
  ENDMACRO(OPENAL_CHECK_LIBRARY)

  IF(OPENAL_INCLUDE_DIR)
    # First look for the OpenAL library
    SET(OPENAL_LIBRARY_DOC "The directory containing the OpenAL library file")
    FIND_LIBRARY(OPENAL_LIBRARY NAMES openal al OpenAL32 PATHS ${OPENAL_ROOT} ${OPENAL_ROOT}/lib ${OPENAL_ROOT}/libs
                 DOC ${OPENAL_LIBRARY_DOC} NO_DEFAULT_PATH)
    IF(NOT OPENAL_LIBRARY)  # now look in system locations
      FIND_LIBRARY(OPENAL_LIBRARY NAMES openal al OpenAL32 DOC ${OPENAL_LIBRARY_DOC})
    ENDIF(NOT OPENAL_LIBRARY)

    IF(OPENAL_LIBRARY)
      # Check if this is sufficient, or whether we need to link in SDL
      OPENAL_CHECK_LIBRARY(OPENAL_LIBRARY_OK)

      IF(OPENAL_LIBRARY_OK)
        SET(OPENAL_FOUND TRUE)
      ELSE(OPENAL_LIBRARY_OK)
        # We need to additionally link in SDL
        FIND_PACKAGE(SDL)

        IF(SDL_FOUND)
          SET(OPENAL_LIBRARY ${OPENAL_LIBRARY} ${SDL_LIBRARY})
          OPENAL_CHECK_LIBRARY(OPENAL_SDL_LIBRARY_OK)  # CHECK_FUNCTION_EXISTS caches result, so must use different variable

          IF(OPENAL_SDL_LIBRARY_OK)
            SET(OPENAL_FOUND TRUE)
          ENDIF(OPENAL_SDL_LIBRARY_OK)
        ENDIF(SDL_FOUND)

      ENDIF(OPENAL_LIBRARY_OK)
    ENDIF(OPENAL_LIBRARY)
  ENDIF(OPENAL_INCLUDE_DIR)

ENDIF(NOT OPENAL_FOUND)


# Now look for ALUT
SET(ALUT_FOUND FALSE)

IF(NOT ALUT_SKIP)
  SET(ALUT_INCLUDE_DOC "The directory containing the ALUT include file alut.h")
  FIND_PATH(ALUT_INCLUDE_DIR NAMES alut.h PATH_SUFFIXES "" "AL" "OpenAL" "ALUT" "alut" "ALut"
            PATHS ${ALUT_ROOT} ${ALUT_ROOT}/include DOC ${ALUT_INCLUDE_DOC} NO_DEFAULT_PATH)
  IF(NOT ALUT_INCLUDE_DIR)  # now look in system locations
    FIND_PATH(ALUT_INCLUDE_DIR NAMES alut.h PATH_SUFFIXES "" "AL" "OpenAL" "ALUT" "alut" "ALut" DOC ${ALUT_INCLUDE_DOC})
  ENDIF(NOT ALUT_INCLUDE_DIR)

  IF(ALUT_INCLUDE_DIR)
    SET(ALUT_LIBRARY_DOC "The directory containing the ALUT library file")
    FIND_LIBRARY(ALUT_LIBRARY NAMES alut ALut ALut32 PATHS ${ALUT_ROOT} ${ALUT_ROOT}/lib ${ALUT_ROOT}/libs
                 DOC ${ALUT_LIBRARY_DOC} NO_DEFAULT_PATH)
    IF(NOT ALUT_LIBRARY)  # now look in system locations
      FIND_LIBRARY(ALUT_LIBRARY NAMES alut ALut ALut32 DOC ${ALUT_LIBRARY_DOC})
    ENDIF(NOT ALUT_LIBRARY)

    IF(ALUT_LIBRARY)
      SET(ALUT_FOUND TRUE)
    ENDIF(ALUT_LIBRARY)
  ENDIF(ALUT_INCLUDE_DIR)
ENDIF(NOT ALUT_SKIP)


IF(OPENAL_FOUND)
  IF(NOT OPENAL_FIND_QUIETLY)
    MESSAGE(STATUS "Found OpenAL: headers at ${OPENAL_INCLUDE_DIR}, libraries at ${OPENAL_LIBRARY}")
  ENDIF(NOT OPENAL_FIND_QUIETLY)
ELSE(OPENAL_FOUND)
  IF(OPENAL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "OpenAL not found")
  ENDIF(OPENAL_FIND_REQUIRED)
ENDIF(OPENAL_FOUND)

IF(NOT ALUT_SKIP)
  IF(ALUT_FOUND)
    IF(NOT OPENAL_FIND_QUIETLY)
      MESSAGE(STATUS "Found ALUT: headers at ${ALUT_INCLUDE_DIR}, libraries at ${ALUT_LIBRARY}")
    ENDIF(NOT OPENAL_FIND_QUIETLY)
  ENDIF(ALUT_FOUND)
ENDIF(NOT ALUT_SKIP)
