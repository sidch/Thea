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
SET(Thea_FIND_CGAL       TRUE)
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

# Dependency: OpenCL
FIND_PACKAGE(OpenCL)

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX)

  STRING(REPLACE ";" " " EXTRA_DEBUG_CFLAGS "${CGAL_DEBUG_CFLAGS}")
  STRING(REPLACE ";" " " EXTRA_RELEASE_CFLAGS "${CGAL_RELEASE_CFLAGS}")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fno-strict-aliasing")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS} -g2")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -O2")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS} -DNDEBUG -g2 -O2")

ELSEIF(MSVC)
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ENDIF()

# Shared library flags
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE
                    ${CMAKE_CURRENT_BINARY_DIR}  # QT4_WRAP_UI generates ui_*.h files here
                    ${SourceRoot}
                    ${Boost_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS}
                    ${CGAL_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${OutputRoot}/lib)

# List of all tools
SET(TheaAllTools )

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

  SET(TheaAllTools ${TheaAllTools} Browse3D)

ENDIF()

#===========================================================
# GIP
#===========================================================

IF(OPENCL_FOUND)
  # Source file lists
  FILE(GLOB GIPSources
       ${SourceRoot}/Tools/GIP/*.cpp)

  # Libraries to link to
  SET(GIPLibraries
      Thea
      ${OPENCL_LIBRARIES}
      ${PLATFORM_LIBRARIES})

  # Build products
  ADD_EXECUTABLE(gip ${GIPSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(gip ${GIPLibraries})

  # Linker flags
  IF(Thea_LDFLAGS)
    SET_TARGET_PROPERTIES(gip PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
  ENDIF()

  # Fix library install names on OS X
  IF(APPLE)
    INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
    OSX_FIX_DYLIB_REFERENCES(GIP "${GIPLibraries}")
  ENDIF()

  SET(TheaAllTools ${TheaAllTools} gip)

ENDIF()

#===========================================================
# ImageFeatures
#===========================================================

# Source file lists
FILE(GLOB ImageFeaturesSources
     ${SourceRoot}/Tools/ImageFeatures/*.cpp)

# Libraries to link to
SET(ImageFeaturesLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(ImageFeatures ${ImageFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(ImageFeatures ${ImageFeaturesLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(ImageFeatures PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(ImageFeatures "${ImageFeaturesLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} ImageFeatures)

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

SET(TheaAllTools ${TheaAllTools} MeshConv)

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

SET(TheaAllTools ${TheaAllTools} MeshFeatures)

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

SET(TheaAllTools ${TheaAllTools} MeshFix)

#===========================================================
# MeshSample
#===========================================================

# Source file lists
FILE(GLOB MeshSampleSources
     ${SourceRoot}/Tools/MeshSample/*.cpp)

# Libraries to link to
SET(MeshSampleLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshSample ${MeshSampleSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshSample ${MeshSampleLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshSample PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshSample "${MeshSampleLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} MeshSample)

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

SET(TheaAllTools ${TheaAllTools} Register)

#===========================================================
# SampleGraph
#===========================================================

# Source file lists
FILE(GLOB SampleGraphSources
     ${SourceRoot}/Tools/SampleGraph/*.cpp)

# Libraries to link to
SET(SampleGraphLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(SampleGraph ${SampleGraphSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(SampleGraph ${SampleGraphLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(SampleGraph PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(SampleGraph "${SampleGraphLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} SampleGraph)

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

SET(TheaAllTools ${TheaAllTools} SegmentSDF)

#===========================================================
# Target for all tools
#===========================================================

ADD_CUSTOM_TARGET(TheaTools DEPENDS ${TheaAllTools})
