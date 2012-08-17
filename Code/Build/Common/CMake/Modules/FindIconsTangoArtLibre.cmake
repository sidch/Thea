# Searches for Tango Art Libre icons
#
# Defines:
#
#   IconsTangoArtLibre_FOUND  True if the Tango Art Libre icons were found, else false
#   IconsTangoArtLibre_DIR    Directory containing the Tango Art Libre icons (i.e. the directory containing 16x16)
#
# To specify an additional directory to search, set IconsTangoArtLibre_ROOT. This should be the directory containing
# TangoArtLibre/16x16 or tango-art-libre/16x16.
#
# Copyright (C) Siddhartha Chaudhuri, 2010
#

# Look for the 16x16 subdirectory, first in the user-specified location and then in the system locations
SET(IconsTangoArtLibre_DOC "The directory containing the Tango Art Libre icons subdirectory 16x16")
# ... though we actually look for {TangoArtLibre, tango-art-libre}/16x16

FIND_PATH(IconsTangoArtLibre_DIR NAMES TangoArtLibre/16x16 tango-art-libre/16x16
          PATHS ${IconsTangoArtLibre_ROOT}
          PATH_SUFFIXES share share/icons
          DOC ${IconsTangoArtLibre_DOC} NO_DEFAULT_PATH)
IF(NOT IconsTangoArtLibre_DIR)  # now look in system locations
  FIND_PATH(IconsTangoArtLibre_DIR NAMES TangoArtLibre/16x16 tango-art-libre/16x16
            PATH_SUFFIXES share share/icons
            DOC ${IconsTangoArtLibre_DOC})
ENDIF(NOT IconsTangoArtLibre_DIR)

IF(IconsTangoArtLibre_DIR)
  SET(IconsTangoArtLibre_FOUND TRUE)

  # Add TangoArtLibre or tango-art-libre as appropriate to the final path
  IF(EXISTS ${IconsTangoArtLibre_DIR}/TangoArtLibre)
    SET(IconsTangoArtLibre_DIR ${IconsTangoArtLibre_DIR}/TangoArtLibre CACHE PATH ${IconsTangoArtLibre_DOC} FORCE)
  ELSEIF(EXISTS ${IconsTangoArtLibre_DIR}/tango-art-libre)
    SET(IconsTangoArtLibre_DIR ${IconsTangoArtLibre_DIR}/tango-art-libre CACHE PATH ${IconsTangoArtLibre_DOC} FORCE)
  ENDIF(EXISTS ${IconsTangoArtLibre_DIR}/TangoArtLibre)

ELSE(IconsTangoArtLibre_DIR)
  SET(IconsTangoArtLibre_FOUND FALSE)
ENDIF(IconsTangoArtLibre_DIR)

IF(IconsTangoArtLibre_FOUND)
  IF(NOT IconsTangoArtLibre_FIND_QUIETLY)
    MESSAGE(STATUS "Found Tango Art Libre icons at ${IconsTangoArtLibre_DIR}")
  ENDIF(NOT IconsTangoArtLibre_FIND_QUIETLY)
ELSE(IconsTangoArtLibre_FOUND)
  IF(IconsTangoArtLibre_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Tango Art Libre icons not found")
  ENDIF(IconsTangoArtLibre_FIND_REQUIRED)
ENDIF(IconsTangoArtLibre_FOUND)
