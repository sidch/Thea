# Searches for Tango icons
#
# Defines:
#
#   IconsTango_FOUND  True if the Tango icons were found, else false
#   IconsTango_DIR    Directory containing the Tango icons (i.e. the directory containing 16x16)
#
# To specify an additional directory to search, set IconsTango_ROOT. This should be the directory containing Tango/16x16 (or
# tango/16x16).
#
# Copyright (C) Siddhartha Chaudhuri, 2010
#

# Look for the 16x16 subdirectory, first in the user-specified location and then in the system locations
SET(IconsTango_DOC "The directory containing the Tango icons subdirectory 16x16")  # though we actually look for Tango/16x16
FIND_PATH(IconsTango_DIR NAMES Tango/16x16 tango/16x16
          PATHS ${IconsTango_ROOT}
          PATH_SUFFIXES share share/icons
          DOC ${IconsTango_DOC} NO_DEFAULT_PATH)
IF(NOT IconsTango_DIR)  # now look in system locations
  FIND_PATH(IconsTango_DIR NAMES Tango/16x16 tango/16x16
            PATH_SUFFIXES share share/icons
            DOC ${IconsTango_DOC})
ENDIF(NOT IconsTango_DIR)

IF(IconsTango_DIR)
  SET(IconsTango_FOUND TRUE)

  # Add Tango or tango as appropriate to the final path
  IF(EXISTS ${IconsTango_DIR}/Tango)
    SET(IconsTango_DIR ${IconsTango_DIR}/Tango CACHE PATH ${IconsTango_DOC} FORCE)
  ELSEIF(EXISTS ${IconsTango_DIR}/tango)
    SET(IconsTango_DIR ${IconsTango_DIR}/tango CACHE PATH ${IconsTango_DOC} FORCE)
  ENDIF(EXISTS ${IconsTango_DIR}/Tango)

ELSE(IconsTango_DIR)
  SET(IconsTango_FOUND FALSE)
ENDIF(IconsTango_DIR)

IF(IconsTango_FOUND)
  IF(NOT IconsTango_FIND_QUIETLY)
    MESSAGE(STATUS "Found Tango icons at ${IconsTango_DIR}")
  ENDIF(NOT IconsTango_FIND_QUIETLY)
ELSE(IconsTango_FOUND)
  IF(IconsTango_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Tango icons not found")
  ENDIF(IconsTango_FIND_REQUIRED)
ENDIF(IconsTango_FOUND)
