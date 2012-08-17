# - Get all the dependencies of a library from its accompanying libtool archive (.la file)
# Any .la dependency is recursively searched.
#
# Example: GET_LIBTOOLIZED_DEPENDENCIES("libhello.la" HELLO_MAIN_ARCHIVE HELLO_DEPENDENCIES)
#   HELLO_MAIN_ARCHIVE is now libhello.a
#   HELLO_DEPENDENCIES is now -lX11;-lc;/usr/local/lib/libzip.a;/usr/lib/libcurl.so
#
# Author: Siddhartha Chaudhuri, 2008
#

MACRO(GET_LIBTOOLIZED_DEPENDENCIES libtool_archive main_archive dependencies)
  _GET_LIBTOOLIZED_DEPENDENCIES_NON_RECURSIVE(${libtool_archive} ${main_archive} ${dependencies})

  SET(_HAS_UNRESOLVED_LA TRUE)
  WHILE(_HAS_UNRESOLVED_LA)
    SET(_HAS_UNRESOLVED_LA FALSE)

    IF(${dependencies})
      LIST(LENGTH ${dependencies} _LA_NUM_DEPS)
      SET(_LA_DEP_INDEX 0)

      # Cycle through the dependencies, recursing if we find another libtoolized archive
      WHILE(${_LA_DEP_INDEX} LESS ${_LA_NUM_DEPS})
        LIST(GET ${dependencies} ${_LA_DEP_INDEX} _LA_DEP)  # get the "i'th" dependency

        IF("${_LA_DEP}" MATCHES ".la$")                                      # the dependency is itself a libtoolized archive
          # Get the dependencies of the dependency
          _GET_LIBTOOLIZED_DEPENDENCIES_NON_RECURSIVE(${_LA_DEP} _LA_DEP_MAIN_ARCHIVE _LA_DEP_DEPS)
          SET(_LA_DEP_DEPS ${_LA_DEP_MAIN_ARCHIVE} ${_LA_DEP_DEPS})

          LIST(REMOVE_AT ${dependencies} ${_LA_DEP_INDEX})                   # get rid of the .la entry
          IF(_LA_DEP_DEPS)                                                   # replace it with the new dependencies
            LIST(LENGTH ${dependencies} _LA_NUM_DEPS)                        # - recompute the number of dependencies
            IF(${_LA_DEP_INDEX} EQUAL ${_LA_NUM_DEPS})                       # - we removed the last element
              LIST(APPEND ${dependencies} ${_LA_DEP_DEPS})                   # - append the new dependencies
            ELSE(${_LA_DEP_INDEX} EQUAL ${_LA_NUM_DEPS})
              LIST(INSERT ${dependencies} ${_LA_DEP_INDEX} ${_LA_DEP_DEPS})  # - insert  the new dependencies
            ENDIF(${_LA_DEP_INDEX} EQUAL ${_LA_NUM_DEPS})
          ENDIF(_LA_DEP_DEPS)
          SET(_HAS_UNRESOLVED_LA TRUE)                                       # maybe we have new .la dependencies
          SET(_LA_DEP_INDEX ${_LA_NUM_DEPS})                                 # break the loop
        ELSE("${_LA_DEP}" MATCHES ".la$")
          MATH(EXPR _LA_DEP_INDEX "${_LA_DEP_INDEX} + 1")                    # move to the next dependency
        ENDIF("${_LA_DEP}" MATCHES ".la$")
      ENDWHILE(${_LA_DEP_INDEX} LESS ${_LA_NUM_DEPS})
    ENDIF(${dependencies})
  ENDWHILE(_HAS_UNRESOLVED_LA)
ENDMACRO(GET_LIBTOOLIZED_DEPENDENCIES)


MACRO(_GET_LIBTOOLIZED_DEPENDENCIES_NON_RECURSIVE nr_libtool_archive nr_main_archive nr_dependencies)
  SET(${nr_main_archive})
  SET(${nr_dependencies})
  FILE(READ ${nr_libtool_archive} _NR_LA_CONTENTS)

  # Get the main archive name
  STRING(REGEX REPLACE ".*old_library[ \t]*=[ \t]*'([^']*).*" "\\1" ${nr_main_archive} "${_NR_LA_CONTENTS}")
  IF(NOT ${nr_main_archive})
    # The library is shared, so look for the dlname entry
    STRING(REGEX REPLACE ".*dlname[ \t]*=[ \t]*'([^']*).*" "\\1" ${nr_main_archive} "${_NR_LA_CONTENTS}")
  ENDIF(NOT ${nr_main_archive})

  IF(${nr_main_archive})
    IF(NOT "${${nr_main_archive}}" MATCHES "/")  # relative filename was specified (probably a needless check)
      GET_FILENAME_COMPONENT(_NR_LA_PATH "${nr_libtool_archive}" ABSOLUTE)
      GET_FILENAME_COMPONENT(_NR_LA_PATH "${_NR_LA_PATH}" PATH)

      IF(_NR_LA_PATH)
        SET(${nr_main_archive} "${_NR_LA_PATH}/${${nr_main_archive}}")
      ENDIF(_NR_LA_PATH)
    ENDIF(NOT "${${nr_main_archive}}" MATCHES "/")
  ENDIF(${nr_main_archive})

  # Get the list of dependencies
  STRING(REGEX REPLACE ".*dependency_libs[ \t]*=[ \t]*'([^']*).*" "\\1" ${nr_dependencies} "${_NR_LA_CONTENTS}")
  SEPARATE_ARGUMENTS(${nr_dependencies})
ENDMACRO(_GET_LIBTOOLIZED_DEPENDENCIES_NON_RECURSIVE)
