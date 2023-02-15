#===============================================================================================================================
#
# Utility code for finding libraries that Thea depends on
#
# Defines:
#
#   Thea_DEPS_FOUND           True if dependencies were found, else false
#   Thea_DEPS_LIBRARIES       Libraries to link
#   Thea_DEPS_LIBRARY_DIRS    Additional directories for libraries. These do not necessarily correspond to Thea_DEPS_LIBRARIES,
#                             and both variables must be passed to the linker.
#   Thea_DEPS_INCLUDE_DIRS    The directories containing the header files
#   Thea_DEPS_CFLAGS          Extra compiler flags
#   Thea_DEPS_DEBUG_CFLAGS    Extra compiler flags to be used only in debug builds
#   Thea_DEPS_RELEASE_CFLAGS  Extra compiler flags to be used only in release builds
#   Thea_DEPS_LDFLAGS         Extra linker flags
#
# To specify an additional directory to search, set THEA_DEPS_ROOT in the parent script or on the command line.
#
# To specify which libraries to search for, set Thea_FIND_<PackageName> to TRUE in the parent script, or just set Thea_FIND_ALL
# to TRUE to search for all supported dependencies. Note: the Thea_FIND_<PackageName> flags will be set to FALSE when this
# script finishes!
#
# To suppress searching/linking optional dependencies from the command line, set WITH_<PACKAGENAME> (all uppercase) to false. To
# omit all optional dependencies, set LITE to true.
#
# Copyright (C) Siddhartha Chaudhuri, 2011
#
#===============================================================================================================================

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Required unless explicitly omitted
IF(NOT DEFINED Thea_FIND_Eigen3)
  SET(Thea_FIND_Eigen3 TRUE)
ENDIF()

IF(Thea_FIND_ALL)
  SET(Thea_FIND_CGAL       TRUE)
  SET(Thea_FIND_CLUTO      TRUE)
  SET(Thea_FIND_FreeImage  TRUE)
  SET(Thea_FIND_Lib3ds     TRUE)
  # SET(Thea_FIND_HDF5       TRUE)
ELSEIF(Thea_LITE)
  SET(Thea_FIND_CGAL       FALSE)
  SET(Thea_FIND_CLUTO      FALSE)
  SET(Thea_FIND_FreeImage  FALSE)
  SET(Thea_FIND_Lib3ds     FALSE)
  # SET(Thea_FIND_HDF5       FALSE)
ENDIF()

# Optional libraries are enabled by default, unless LITE is true
IF(LITE)
  SET(WITH_OPTIONAL FALSE)
ELSE()
  SET(WITH_OPTIONAL TRUE)
ENDIF()

IF(NOT DEFINED WITH_CGAL)
  SET(WITH_CGAL ${WITH_OPTIONAL})
ENDIF()

IF(NOT DEFINED WITH_CLUTO)
  SET(WITH_CLUTO ${WITH_OPTIONAL})
ENDIF()

IF(NOT DEFINED WITH_FREEIMAGE)
  SET(WITH_FREEIMAGE ${WITH_OPTIONAL})
ENDIF()

IF(NOT DEFINED WITH_LIB3DS)
  SET(WITH_LIB3DS ${WITH_OPTIONAL})
ENDIF()

SET(Thea_DEPS_INCLUDE_DIRS   )
SET(Thea_DEPS_LIBRARIES      )
SET(Thea_DEPS_LIBRARY_DIRS   )
SET(Thea_DEPS_CFLAGS         )  # common flags for all builds
SET(Thea_DEPS_DEBUG_CFLAGS   )  # extra flags for debug builds
SET(Thea_DEPS_RELEASE_CFLAGS )  # extra flags for release builds
SET(Thea_DEPS_LDFLAGS        )

# Dependency: Eigen3 (required)
IF(Thea_FIND_Eigen3)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-eigen3)
    SET(EIGEN3_ROOT ${THEA_DEPS_ROOT}/installed-eigen3)
  ELSE()
    SET(EIGEN3_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  FIND_PACKAGE(Eigen3 REQUIRED)
  SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
ENDIF()

# Dependency: CGAL (optional)
IF(WITH_CGAL AND Thea_FIND_CGAL)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-cgal)
    SET(CGAL_ROOT ${THEA_DEPS_ROOT}/installed-cgal)
  ELSE()
    SET(CGAL_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  FIND_PACKAGE(CGAL)  # not a required package
  IF(CGAL_FOUND)
    SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS})
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CGAL=1")
    SET(Thea_DEPS_DEBUG_CFLAGS "${Thea_DEPS_DEBUG_CFLAGS} ${CGAL_DEBUG_CFLAGS}")
    SET(Thea_DEPS_RELEASE_CFLAGS "${Thea_DEPS_RELEASE_CFLAGS} ${CGAL_RELEASE_CFLAGS}")

    IF(CGAL_LIBRARY)
      SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${CGAL_LIBRARY})
      SET(Thea_DEPS_LIBRARY_DIRS ${Thea_DEPS_LIBRARY_DIRS} ${CGAL_LIBRARY_DIRS})
    ENDIF()

    # CGAL appends the directory containing its own CMake modules to the module search path. We shouldn't need it after this
    # point, so let's drop everything on the module path other than the first component.
    LIST(GET CMAKE_MODULE_PATH 0 CMAKE_MODULE_PATH)

  ELSE()  # this is not a fatal error
    MESSAGE(STATUS "CGAL not found: library will be built without CGAL-dependent geometry processing components")
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CGAL=0")
  ENDIF()
ELSE()
  SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CGAL=0")
ENDIF()

# Dependency: CLUTO (optional)
IF(WITH_CLUTO AND Thea_FIND_CLUTO)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-cluto)
    SET(CLUTO_ROOT ${THEA_DEPS_ROOT}/installed-cluto)
  ELSE()
    SET(CLUTO_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  FIND_PACKAGE(CLUTO)  # not a required package
  IF(CLUTO_FOUND)
    SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${CLUTO_LIBRARIES})
    SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${CLUTO_INCLUDE_DIRS})
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CLUTO=1")
  ELSE()  # this is not a fatal error
    MESSAGE(STATUS "CLUTO not found: library will be built without CLUTO clustering support")
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CLUTO=0")
  ENDIF()
ELSE()
  SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_CLUTO=0")
ENDIF()

# Dependency: FreeImage (optional)
IF(WITH_FREEIMAGE AND Thea_FIND_FreeImage)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-freeimage)
    SET(FreeImage_ROOT ${THEA_DEPS_ROOT}/installed-freeimage)
  ELSE()
    SET(FreeImage_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  SET(FreeImage_LANGUAGE "C++")
  FIND_PACKAGE(FreeImage)

  IF(FreeImage_FOUND)
    SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${FreeImage_LIBRARIES})
    SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${FreeImage_INCLUDE_DIRS})
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_FREEIMAGE=1")
  ELSE()  # this is not a fatal error
    MESSAGE(STATUS "FreeImage not found: library will be built without FreeImage image codecs")
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_FREEIMAGE=0")
  ENDIF()
ELSE()
  SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_FREEIMAGE=0")
ENDIF()

# Dependency: Lib3ds (optional)
IF(WITH_LIB3DS AND Thea_FIND_Lib3ds)
  IF(EXISTS ${THEA_DEPS_ROOT}/installed-lib3ds)
    SET(Lib3ds_ROOT ${THEA_DEPS_ROOT}/installed-lib3ds)
  ELSE()
    SET(Lib3ds_ROOT ${THEA_DEPS_ROOT})
  ENDIF()
  FIND_PACKAGE(Lib3ds)

  IF(Lib3ds_FOUND)
    SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${Lib3ds_LIBRARIES})
    SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${Lib3ds_INCLUDE_DIRS})
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_LIB3DS=1")
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_LIB3DS_VERSION_MAJOR=${Lib3ds_VERSION_MAJOR}")
  ELSE()  # this is not a fatal error
    MESSAGE(STATUS "Lib3ds not found: library will be built without 3DS mesh file support")
    SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_LIB3DS=0")
  ENDIF()
ELSE()
  SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} -DTHEA_ENABLE_LIB3DS=0")
ENDIF()

# # Dependency: HDF5 (optional)
# IF(WITH_HDF5 AND Thea_FIND_HDF5)
#   IF(EXISTS ${THEA_DEPS_ROOT}/installed-hdf5)
#     SET(HDF5_ROOT ${THEA_DEPS_ROOT}/installed-hdf5)
#   ELSE()
#     SET(HDF5_ROOT ${THEA_DEPS_ROOT})
#   ENDIF()
#
#   # Work around a FindHDF5 bug (?) that causes system paths to be ignored if HDF5_ROOT is specified
#   FIND_PACKAGE(HDF5 QUIET)
#   IF(NOT HDF5_FOUND)
#     SET(HDF5_ROOT )
#     FIND_PACKAGE(HDF5)
#   ENDIF()
#
#   IF(HDF5_FOUND)
#     SET(Thea_DEPS_INCLUDE_DIRS ${Thea_DEPS_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
#     SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} ${HDF5_DEFINITIONS} -DTHEA_ENABLE_HDF5=1")
#     SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${HDF5_LIBRARIES})
#     SET(Thea_DEPS_LIBRARY_DIRS ${Thea_DEPS_LIBRARY_DIRS} ${HDF5_LIBRARY_DIRS})
#   ELSE()  # this is not a fatal error
#     MESSAGE(STATUS "HDF5 not found: library will be built without HDF5 matrix codecs")
#     SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} ${HDF5_DEFINITIONS} -DTHEA_ENABLE_HDF5=0")
#   ENDIF()
# ELSE()
#     SET(Thea_DEPS_CFLAGS "${Thea_DEPS_CFLAGS} ${HDF5_DEFINITIONS} -DTHEA_ENABLE_HDF5=0")
# ENDIF()

# Additional platform-specific libraries
FIND_PACKAGE(Threads REQUIRED)
SET(Thea_DEPS_PLATFORM_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(Thea_DEPS_PLATFORM_LIBRARIES ${Thea_DEPS_PLATFORM_LIBRARIES} "-framework Carbon")
ENDIF()

SET(Thea_DEPS_PLATFORM_LIBRARIES ${Thea_DEPS_PLATFORM_LIBRARIES} ${CMAKE_DL_LIBS})  # for loading plugins with DynLib
SET(Thea_DEPS_LIBRARIES ${Thea_DEPS_LIBRARIES} ${Thea_DEPS_PLATFORM_LIBRARIES})

# Unset parameters
SET(Thea_FIND_Eigen3     FALSE)
SET(Thea_FIND_CGAL       FALSE)
SET(Thea_FIND_CLUTO      FALSE)
SET(Thea_FIND_FreeImage  FALSE)
SET(Thea_FIND_Lib3ds     FALSE)
# SET(Thea_FIND_HDF5       FALSE)
SET(Thea_FIND_ALL        FALSE)
