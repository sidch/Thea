#===============================================================================================================================
#
# Utility code for finding libraries that Thea depends on
#
# Copyright (C) Siddhartha Chaudhuri, 2011
#
#===============================================================================================================================

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

IF(Thea_FIND_ALL)
  SET(Thea_FIND_Boost      TRUE)
  SET(Thea_FIND_FreeImage  TRUE)
  SET(Thea_FIND_CGAL       TRUE)
  SET(Thea_FIND_Lib3ds     TRUE)
  SET(Thea_FIND_CLUTO      TRUE)
ENDIF()

# Dependency: Boost
IF(Thea_FIND_Boost OR Thea_FIND_CGAL)
  SET(Boost_USE_STATIC_LIBS      ON)
  SET(Boost_USE_MULTITHREADED    ON)
  SET(Boost_USE_STATIC_RUNTIME  OFF)
  INCLUDE(BoostAdditionalVersions)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-boost)
    SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-boost)
  ELSE()
    SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  FIND_PACKAGE(Boost COMPONENTS filesystem system thread REQUIRED)
ENDIF()

# Dependency: FreeImage
IF(Thea_FIND_FreeImage)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-freeimage)
    SET(FreeImage_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-freeimage)
  ELSE()
    SET(FreeImage_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  SET(FreeImage_LANGUAGE "C++")
  FIND_PACKAGE(FreeImage REQUIRED)
ENDIF()

# Dependency: CGAL
IF(Thea_FIND_CGAL)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-cgal)
    SET(CGAL_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-cgal)
  ELSE()
    SET(CGAL_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  FIND_PACKAGE(CGAL REQUIRED)
  ADD_DEFINITIONS(${CGAL_3RD_PARTY_DEFINITIONS})
  SET(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS} ${CGAL_3RD_PARTY_INCLUDE_DIRS})
  SET(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES})
  IF(NOT CGAL_LIBRARY)
    MESSAGE(STATUS "CGAL libraries will be auto-linked")
    LINK_DIRECTORIES(${CGAL_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
  ENDIF()
  # -O3 might cause problems with old versions of gcc 4 on OS X
  LIST(REMOVE_ITEM CGAL_RELEASE_CFLAGS "-O3")
  ADD_DEFINITIONS(-DTHEA_ENABLE_CGAL)
ENDIF()

# Dependency: Lib3ds
IF(Thea_FIND_Lib3ds)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-lib3ds)
    SET(Lib3ds_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-lib3ds)
  ELSE()
    SET(Lib3ds_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  FIND_PACKAGE(Lib3ds REQUIRED)
  ADD_DEFINITIONS(-DTHEA_LIB3DS_VERSION_MAJOR=${Lib3ds_VERSION_MAJOR})
ENDIF()

# Dependency: CLUTO
IF(Thea_FIND_CLUTO)
  IF(APPLE)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 4)  # We don't have Cluto built for 64-bit OS X, i.e. Snow Leopard (10.6)
      SET(_Thea_FIND_CLUTO TRUE)
    ELSE()
      SET(_Thea_FIND_CLUTO FALSE)
    ENDIF()
  ELSE()
    SET(_Thea_FIND_CLUTO TRUE)
  ENDIF()

  IF(_Thea_FIND_CLUTO)
    IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-cluto)
      SET(CLUTO_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-cluto)
    ELSE()
      SET(CLUTO_ROOT ${THEA_INSTALLATIONS_ROOT})
    ENDIF()
    FIND_PACKAGE(CLUTO REQUIRED)
  ELSE()
    MESSAGE(STATUS "NOTE: CLUTO not available for this system: Thea will be built without CLUTO clustering support")
  ENDIF()
ENDIF()

# Additional platform-specific libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(PLATFORM_LIBRARIES "-framework Carbon")
ENDIF()

SET(PLATFORM_LIBRARIES ${PLATFORM_LIBRARIES} ${CMAKE_DL_LIBS})  # for loading plugins with DynLib

# Unset parameters
SET(Thea_FIND_Boost      FALSE)
SET(Thea_FIND_G3D        FALSE)
SET(Thea_FIND_FreeImage  FALSE)
SET(Thea_FIND_CGAL       FALSE)
SET(Thea_FIND_Lib3ds     FALSE)
SET(Thea_FIND_CLUTO      FALSE)
SET(Thea_FIND_ALL        FALSE)
