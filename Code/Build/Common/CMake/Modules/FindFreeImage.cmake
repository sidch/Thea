# Searches for an installation of the FreeImage library. On success, it sets the following variables:
#
#   FreeImage_FOUND              Set to true to indicate the library was found
#   FreeImage_INCLUDE_DIRS       The directory containing the main FreeImage header file FreeImage.h
#   FreeImage_LIBRARY_DIRS       The directory containing the FreeImage libraries (for info only, should not be required)
#   FreeImage_LIBRARIES          The libraries needed to use FreeImage (with full paths)
#   FreeImage_VERSION_MAJOR      The major version of the installation
#   FreeImage_VERSION_MINOR      The minor version of the installation
#   FreeImage_VERSION_PATCH      The patch/subminor version of the installation
#
# To specify an additional directory to search, set FreeImage_ROOT.
# To search for the C++ wrapper FreeImagePlus instead, set FreeImage_LANGUAGE to "C++".
#
# Copyright (C) Siddhartha Chaudhuri, 2008
#

IF(FreeImage_LANGUAGE MATCHES "[Cc][+][+]")
  SET(FreeImage_LANGUAGE_INTERNAL "C++")
  SET(FreeImage_HEADER_NAME FreeImagePlus.h)
  IF(WIN32)
    SET(FreeImage_DEBUG_LIBRARY_NAME FreeImagePlusd)
    SET(FreeImage_RELEASE_LIBRARY_NAME FreeImagePlus)
  ELSE(WIN32)
    SET(FreeImage_LIBRARY_NAME freeimageplus)
  ENDIF(WIN32)
ELSE(FreeImage_LANGUAGE MATCHES "[Cc][+][+]")
  SET(FreeImage_LANGUAGE_INTERNAL "C")
  SET(FreeImage_HEADER_NAME FreeImage.h)
  IF(WIN32)
    SET(FreeImage_DEBUG_LIBRARY_NAME FreeImaged)
    SET(FreeImage_RELEASE_LIBRARY_NAME FreeImage)
  ELSE(WIN32)
    SET(FreeImage_LIBRARY_NAME freeimage)
  ENDIF(WIN32)
ENDIF(FreeImage_LANGUAGE MATCHES "[Cc][+][+]")

# Look for the FreeImage header, first in the user-specified location and then in the system locations
SET(FreeImage_INCLUDE_DOC "The directory containing the FreeImage header file FreeImage.h")
FIND_PATH(FreeImage_INCLUDE_DIRS NAMES ${FreeImage_HEADER_NAME} PATHS ${FreeImage_ROOT} ${FreeImage_ROOT}/include
          DOC ${FreeImage_INCLUDE_DOC} NO_DEFAULT_PATH)
IF(NOT FreeImage_INCLUDE_DIRS)  # now look in system locations
  FIND_PATH(FreeImage_INCLUDE_DIRS NAMES ${FreeImage_HEADER_NAME} DOC ${FreeImage_INCLUDE_DOC})
ENDIF(NOT FreeImage_INCLUDE_DIRS)

SET(FreeImage_FOUND FALSE)

IF(FreeImage_INCLUDE_DIRS)
  SET(FreeImage_LIBRARY_DIRS ${FreeImage_INCLUDE_DIRS})

  IF("${FreeImage_LIBRARY_DIRS}" MATCHES "/include$")
    # Strip off the trailing "/include" in the path.
    GET_FILENAME_COMPONENT(FreeImage_LIBRARY_DIRS ${FreeImage_LIBRARY_DIRS} PATH)
  ENDIF("${FreeImage_LIBRARY_DIRS}" MATCHES "/include$")

  IF(EXISTS "${FreeImage_LIBRARY_DIRS}/lib")
    SET(FreeImage_LIBRARY_DIRS ${FreeImage_LIBRARY_DIRS}/lib)
  ENDIF(EXISTS "${FreeImage_LIBRARY_DIRS}/lib")

  IF(WIN32)
    FILE(GLOB FreeImage_ALL_DEBUG_LIBS   ${FreeImage_LIBRARY_DIRS}/${FreeImage_DEBUG_LIBRARY_NAME}.lib)
    FILE(GLOB FreeImage_ALL_RELEASE_LIBS ${FreeImage_LIBRARY_DIRS}/${FreeImage_RELEASE_LIBRARY_NAME}.lib)

    IF(FreeImage_ALL_DEBUG_LIBS)
      LIST(GET FreeImage_ALL_DEBUG_LIBS 0 FreeImage_DEBUG_LIBRARY)
    ELSE(FreeImage_ALL_DEBUG_LIBS)
      SET(FreeImage_DEBUG_LIBRARY )
    ENDIF(FreeImage_ALL_DEBUG_LIBS)

    IF(FreeImage_ALL_RELEASE_LIBS)
      LIST(GET FreeImage_ALL_RELEASE_LIBS 0 FreeImage_RELEASE_LIBRARY)
    ELSE(FreeImage_ALL_RELEASE_LIBS)
      SET(FreeImage_RELEASE_LIBRARY )
    ENDIF(FreeImage_ALL_RELEASE_LIBS)

    SET(FreeImage_LIBRARIES)
    IF(FreeImage_DEBUG_LIBRARY AND FreeImage_RELEASE_LIBRARY)
      SET(FreeImage_LIBRARIES debug ${FreeImage_DEBUG_LIBRARY} optimized ${FreeImage_RELEASE_LIBRARY})
    ELSEIF(FreeImage_DEBUG_LIBRARY)
      SET(FreeImage_LIBRARIES ${FreeImage_DEBUG_LIBRARY})
    ELSEIF(FreeImage_RELEASE_LIBRARY)
      SET(FreeImage_LIBRARIES ${FreeImage_RELEASE_LIBRARY})
    ENDIF(FreeImage_DEBUG_LIBRARY AND FreeImage_RELEASE_LIBRARY)
  ELSE(WIN32)
    FILE(GLOB FreeImage_ALL_LIBS ${FreeImage_LIBRARY_DIRS}/lib${FreeImage_LIBRARY_NAME}*.so
                                 ${FreeImage_LIBRARY_DIRS}/lib${FreeImage_LIBRARY_NAME}*.so.*
                                 ${FreeImage_LIBRARY_DIRS}/lib${FreeImage_LIBRARY_NAME}*.dylib
                                 ${FreeImage_LIBRARY_DIRS}/lib${FreeImage_LIBRARY_NAME}*.dylib.*
                                 ${FreeImage_LIBRARY_DIRS}/lib${FreeImage_LIBRARY_NAME}*.a)
    IF(FreeImage_ALL_LIBS)
      LIST(GET FreeImage_ALL_LIBS 0 FreeImage_LIBRARIES)
    ELSE(FreeImage_ALL_LIBS)
      MESSAGE(STATUS "Couldn't find FreeImage")
      SET(FreeImage_LIBRARIES )
    ENDIF(FreeImage_ALL_LIBS)
  ENDIF(WIN32)

  SET(FreeImage_LIBRARIES ${FreeImage_LIBRARIES} CACHE FILEPATH "The location of the FreeImage library")

  # The library version can be determined by parsing the file FreeImage.h (NOT FreeImagePlus.h even when searching for the C++
  # wrapper -- FreeImage.h should still be present in the same directory since it is included by FreeImagePlus.h)
  IF(FreeImage_LIBRARIES)
    SET(FreeImage_VERSION_MAJOR)
    SET(FreeImage_VERSION_MINOR)
    SET(FreeImage_VERSION_PATCH)
    FILE(READ "${FreeImage_INCLUDE_DIRS}/FreeImage.h" _FreeImage_H_CONTENTS)

    STRING(REGEX MATCH "#define[ \t]+FREEIMAGE_MAJOR_VERSION[ \t]+[0-9]+" FreeImage_VERSION_MAJOR "${_FreeImage_H_CONTENTS}")
    STRING(REGEX MATCH "[0-9]+$" FreeImage_VERSION_MAJOR ${FreeImage_VERSION_MAJOR})

    STRING(REGEX MATCH "#define[ \t]+FREEIMAGE_MINOR_VERSION[ \t]+[0-9]+" FreeImage_VERSION_MINOR "${_FreeImage_H_CONTENTS}")
    STRING(REGEX MATCH "[0-9]+$" FreeImage_VERSION_MINOR ${FreeImage_VERSION_MINOR})

    STRING(REGEX MATCH "#define[ \t]+FREEIMAGE_RELEASE_SERIAL[ \t]+[0-9]+" FreeImage_VERSION_PATCH
           "${_FreeImage_H_CONTENTS}")
    STRING(REGEX MATCH "[0-9]+$" FreeImage_VERSION_PATCH ${FreeImage_VERSION_PATCH})

    SET(FreeImage_VERSION_MAJOR ${FreeImage_VERSION_MAJOR} CACHE INTERNAL "The major version number for the FreeImage library")
    SET(FreeImage_VERSION_MINOR ${FreeImage_VERSION_MINOR} CACHE INTERNAL "The minor version number for the FreeImage library")
    SET(FreeImage_VERSION_PATCH ${FreeImage_VERSION_PATCH}
        CACHE INTERNAL "The patch/subminor version number for the FreeImage library")

    SET(FreeImage_FOUND TRUE)
  ENDIF(FreeImage_LIBRARIES)
ENDIF(FreeImage_INCLUDE_DIRS)

IF(FreeImage_FOUND)
  IF(NOT FreeImage_FIND_QUIETLY)
    MESSAGE(STATUS "Found FreeImage (language ${FreeImage_LANGUAGE_INTERNAL}):"
                   " headers at ${FreeImage_INCLUDE_DIRS}, libraries at ${FreeImage_LIBRARIES}")
  ENDIF(NOT FreeImage_FIND_QUIETLY)
ELSE(FreeImage_FOUND)
  IF(FreeImage_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "FreeImage (language ${FreeImage_LANGUAGE_INTERNAL}) not found")
  ENDIF(FreeImage_FIND_REQUIRED)
ENDIF(FreeImage_FOUND)
