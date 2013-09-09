#===============================================================================================================================
#
# Build script for the Thea tools.
#
# Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaTools)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# See cmake --help-policy CMP0003 for details on this one
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

# If you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Postfix for debug builds
SET(CMAKE_DEBUG_POSTFIX "d")

# Project root path
GET_FILENAME_COMPONENT(ProjectRoot ../.. ABSOLUTE)

# Path for build products
SET(OutputRoot ${ProjectRoot}/Build/Output)

# Path to put executables in
SET(EXECUTABLE_OUTPUT_PATH ${OutputRoot}/bin)

# Path to put libraries in
SET(LIBRARY_OUTPUT_PATH ${OutputRoot}/lib)

# Path for customized CMake modules
IF(NOT CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH ${ProjectRoot}/Build/Common/CMake/Modules)
ENDIF()
GET_FILENAME_COMPONENT(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ABSOLUTE)

# Path to root folder for source code
SET(SourceRoot ${ProjectRoot}/Source)

# Locate Thea dependencies
SET(Thea_FIND_Lib3ds     TRUE)
SET(Thea_FIND_FreeImage  TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Dependency: Boost
SET(Boost_USE_STATIC_LIBS      ON)
SET(Boost_USE_MULTITHREADED    ON)
SET(Boost_USE_STATIC_RUNTIME  OFF)
INCLUDE(BoostAdditionalVersions)
IF(EXISTS ${THEA_INSTALLATIONS_ROOT}/installed-boost)
  SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT}/installed-boost)
ELSE()
  SET(BOOST_ROOT ${THEA_INSTALLATIONS_ROOT})
ENDIF()
FIND_PACKAGE(Boost COMPONENTS program_options thread REQUIRED)

# Dependency: OpenGL
FIND_PACKAGE(OpenGL)
SET(OpenGL_FOUND ${OPENGL_FOUND})
SET(OpenGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
SET(OpenGL_CFLAGS )
SET(OpenGL_LIBRARIES    ${OPENGL_LIBRARIES})
SET(OpenGL_LIBRARY_DIRS )
SET(OpenGL_LDFLAGS )
IF(APPLE)
  ADD_DEFINITIONS(-DTHEA_GL_USE_FRAMEWORK)
ENDIF()
MESSAGE(STATUS "Found OpenGL: headers at ${OpenGL_INCLUDE_DIRS}, libraries at ${OpenGL_LIBRARIES}")

# Dependency: Qt4
FIND_PACKAGE(Qt4)
SET(QT_USE_QTOPENGL 1)
INCLUDE(${QT_USE_FILE})

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(${CGAL_DEBUG_CFLAGS} -Wall -g2 -fno-strict-aliasing -frounding-math)
  ELSE()
    ADD_DEFINITIONS(${CGAL_RELEASE_CFLAGS} -Wall -O2 -fno-strict-aliasing -DNDEBUG)
  ENDIF()
ELSEIF(MSVC)
  IF(NOT CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    ADD_DEFINITIONS(/O2)
  ENDIF()
ENDIF()

# Shared library flags
ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)

# Include directories
INCLUDE_DIRECTORIES(BEFORE
                    ${CMAKE_CURRENT_BINARY_DIR}  # QT4_WRAP_UI generates ui_*.h files here
                    ${SourceRoot}
                    ${Boost_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${OutputRoot}/lib)

#===========================================================
# Browse3D
#===========================================================

IF(QT4_FOUND AND OpenGL_FOUND)
  # Run Qt .ui files through uic -- produces ui_*.h files in current directory
  QT4_WRAP_UI(Browse3D_QtUicSources
              ${SourceRoot}/Tools/Browse3D/MainWindow.ui)

  # Run Qt .hpp files defining QObjects through moc -- produces moc_*.cxx files in source directory
  QT4_WRAP_CPP(Browse3D_QtMocSources
               ${SourceRoot}/Tools/Browse3D/App.hpp
               ${SourceRoot}/Tools/Browse3D/MainWindow.hpp
               ${SourceRoot}/Tools/Browse3D/Model.hpp
               ${SourceRoot}/Tools/Browse3D/ModelDisplay.hpp)

  # Source file lists
  FILE(GLOB Browse3DSources
       ${SourceRoot}/Tools/Browse3D/*.cpp)

  SET(Browse3DSources ${Browse3D_QtUicSources} ${Browse3D_QtMocSources} ${Browse3DSources})

  # Libraries to link to
  SET(Browse3DLibraries
      Thea
      ${QT_LIBRARIES}
      ${Lib3ds_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${Boost_LIBRARIES}
      ${PLATFORM_LIBRARIES})

  # Build products
  ADD_EXECUTABLE(Browse3D ${Browse3DSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(Browse3D ${Browse3DLibraries})

  # Linker flags
  IF(Thea_LDFLAGS)
    SET_TARGET_PROPERTIES(Browse3D PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
  ENDIF()

  # Fix library install names on OS X
  IF(APPLE)
    INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
    OSX_FIX_DYLIB_REFERENCES(Browse3D "${Browse3DLibraries}")
  ENDIF()

  SET(TheaToolsDependencies ${TheaToolsDependencies} Browse3D)

ENDIF()

#===========================================================
# MeshConv
#===========================================================

# Source file lists
FILE(GLOB MeshConvSources
     ${SourceRoot}/Tools/MeshConv/*.cpp)

# Libraries to link to
SET(MeshConvLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshConv ${MeshConvSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshConv ${MeshConvLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshConv PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshConv "${MeshConvLibraries}")
ENDIF()

SET(TheaToolsDependencies ${TheaToolsDependencies} MeshConv)

#===========================================================
# MeshFeatures
#===========================================================

# Source file lists
FILE(GLOB MeshFeaturesSources
     ${SourceRoot}/Tools/MeshFeatures/*.cpp)

# Libraries to link to
SET(MeshFeaturesLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshFeatures ${MeshFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshFeatures ${MeshFeaturesLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshFeatures PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshFeatures "${MeshFeaturesLibraries}")
ENDIF()

SET(TheaToolsDependencies ${TheaToolsDependencies} MeshFeatures)

#===========================================================
# MeshFix
#===========================================================

# Source file lists
FILE(GLOB MeshFixSources
     ${SourceRoot}/Tools/MeshFix/*.cpp)

# Libraries to link to
SET(MeshFixLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshFix ${MeshFixSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshFix ${MeshFixLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshFix PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshFix "${MeshFixLibraries}")
ENDIF()

SET(TheaToolsDependencies ${TheaToolsDependencies} MeshFix)

#===========================================================
# Register
#===========================================================

# Source file lists
FILE(GLOB RegisterSources
     ${SourceRoot}/Tools/Register/*.cpp)

# Libraries to link to
SET(RegisterLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(Register ${RegisterSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(Register ${RegisterLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(Register PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(Register "${RegisterLibraries}")
ENDIF()

SET(TheaToolsDependencies ${TheaToolsDependencies} Register)

#===========================================================
# SegmentSDF
#===========================================================

# Source file lists
FILE(GLOB SegmentSDFSources
     ${SourceRoot}/Tools/SegmentSDF/*.cpp)

# Libraries to link to
SET(SegmentSDFLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(SegmentSDF ${SegmentSDFSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(SegmentSDF ${SegmentSDFLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(SegmentSDF PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(SegmentSDF "${SegmentSDFLibraries}")
ENDIF()

SET(TheaToolsDependencies ${TheaToolsDependencies} SegmentSDF)

#===========================================================
# Target for all tools
#===========================================================

ADD_CUSTOM_TARGET(TheaTools DEPENDS ${TheaToolsDependencies})
