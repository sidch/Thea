# OS X requires a Mach-O dynamic library to have a baked "install name", that is used by other modules to link to it. Depending
# on how the library is built, the install name is not always an absolute path, nor necessarily the same as the name of the
# library file itself. This macro takes as input the name of a target, and a list of libraries that it links to (the output of
# FIND_PACKAGE or FIND_LIBRARY calls), and generates a set of custom, post-build commands that, for each linked dylib, changes
# the name the target uses to refer to it with a fully-qualified (absolute) version of the library's own install name. This
# helps ensure that the target can be used from any location while still being able to locate the linked dynamic libraries.
#
# Note that this script does NOT handle the case when a linked library itself refers to another library using a non-absolute
# name (Boost is a notorious example). To avoid such issues, it is recommended to use a static library instead of a shared one
# in a non-standard location. Alternatively, set DYLD_LIBRARY_PATH to include these non-standard locations when running the
# program (not recommended).
#
# Copyright (C) Siddhartha Chaudhuri, 2009.
#

MACRO(OSX_FIX_DYLIB_REFERENCES target libraries)

  IF(APPLE)
    IF(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
      GET_TARGET_PROPERTY(OFIN_${target}_Output ${target} DEBUG_LOCATION)
      IF(NOT DEBUG_LOCATION)  # variable name changed in 2.6
        GET_TARGET_PROPERTY(OFIN_${target}_Output ${target} LOCATION_DEBUG)
      ENDIF()
    ELSE()
      GET_TARGET_PROPERTY(OFIN_${target}_Output ${target} LOCATION)
    ENDIF()

    FOREACH(OFIN_${target}_Library ${libraries})
      IF(${OFIN_${target}_Library} MATCHES ".*[.]dylib")
        # Resolve symlinks and get absolute path
        GET_FILENAME_COMPONENT(OFIN_${target}_LibraryAbsolute ${OFIN_${target}_Library} ABSOLUTE)

        # Get the baked install name of the library
        EXECUTE_PROCESS(COMMAND otool -D ${OFIN_${target}_LibraryAbsolute}
                        OUTPUT_VARIABLE OFIN_${target}_LibraryInstallNameOutput
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        STRING(REGEX REPLACE "[\r\n]" " " OFIN_${target}_LibraryInstallNameOutput ${OFIN_${target}_LibraryInstallNameOutput})
        SEPARATE_ARGUMENTS(OFIN_${target}_LibraryInstallNameOutput)
        LIST(GET OFIN_${target}_LibraryInstallNameOutput 1 OFIN_${target}_LibraryInstallName)

        # Replace the filename with the install name
        IF(${OFIN_${target}_LibraryInstallName} MATCHES ".*/.*")  # install name has path baked in
          GET_FILENAME_COMPONENT(OFIN_${target}_LibraryAbsolute ${OFIN_${target}_LibraryInstallName} ABSOLUTE)
        ELSE()  # install name is just unqualified filename
          GET_FILENAME_COMPONENT(OFIN_${target}_LibraryAbsolutePath ${OFIN_${target}_LibraryAbsolute} PATH)
          SET(OFIN_${target}_LibraryAbsolute "${OFIN_${target}_LibraryAbsolutePath}/${OFIN_${target}_LibraryInstallName}")
        ENDIF()

        # Replace the unqualified filename, if it appears, with the absolute location
        GET_FILENAME_COMPONENT(OFIN_${target}_LibraryFilename ${OFIN_${target}_LibraryAbsolute} NAME)
        ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD
                           COMMAND install_name_tool
                           ARGS -change
                                ${OFIN_${target}_LibraryFilename}
                                ${OFIN_${target}_LibraryAbsolute}
                                ${OFIN_${target}_Output})
      ENDIF()
    ENDFOREACH(OFIN_${target}_Library)
  ENDIF()

ENDMACRO(OSX_FIX_DYLIB_REFERENCES)
