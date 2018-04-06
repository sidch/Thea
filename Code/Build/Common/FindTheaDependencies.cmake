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
  SET(Thea_FIND_Lib3ds     TRUE)
  SET(Thea_FIND_CLUTO      TRUE)
ENDIF()

# Dependency: Boost
IF(Thea_FIND_Boost)
  SET(Boost_USE_MULTITHREADED    ON)
  # SET(Boost_USE_STATIC_LIBS      ON)
  # SET(Boost_USE_STATIC_RUNTIME  OFF)
  INCLUDE(BoostAdditionalVersions)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-boost)
    SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-boost)
  ELSE()
    SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  IF(NOT Thea_FIND_Boost_COMPONENTS)
    SET(Thea_FIND_Boost_COMPONENTS filesystem system thread)
  ENDIF(NOT Thea_FIND_Boost_COMPONENTS)
  FIND_PACKAGE(Boost COMPONENTS ${Thea_FIND_Boost_COMPONENTS} REQUIRED)
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

# Dependency: CLUTO (optional)
IF(Thea_FIND_CLUTO)
  IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-cluto)
    SET(CLUTO_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-cluto)
  ELSE()
    SET(CLUTO_ROOT ${THEA_INSTALLATIONS_ROOT})
  ENDIF()
  FIND_PACKAGE(CLUTO)  # not a required package
  IF(CLUTO_FOUND)
    ADD_DEFINITIONS(-DTHEA_ENABLE_CLUTO)
  ENDIF()
ENDIF()

# Additional platform-specific libraries
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(PLATFORM_LIBRARIES "-framework Carbon")
ENDIF()

SET(PLATFORM_LIBRARIES ${PLATFORM_LIBRARIES} ${CMAKE_DL_LIBS})  # for loading plugins with DynLib

# Unset parameters
SET(Thea_FIND_Boost      FALSE)
SET(Thea_FIND_Boost_COMPONENTS )
SET(Thea_FIND_FreeImage  FALSE)
SET(Thea_FIND_Lib3ds     FALSE)
SET(Thea_FIND_CLUTO      FALSE)
SET(Thea_FIND_ALL        FALSE)
