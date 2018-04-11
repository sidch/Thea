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
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

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
SET(Thea_FIND_ALL               FALSE)
SET(Thea_FIND_Boost             TRUE)
SET(Thea_FIND_Boost_COMPONENTS  filesystem program_options system thread)
SET(Thea_FIND_CLUTO             TRUE)
SET(Thea_FIND_FreeImage         TRUE)
SET(Thea_FIND_Lib3ds            TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Dependency: OpenGL
IF(THEA_GL_OSMESA)
  FIND_PACKAGE(OSMesa)

  IF(OSMesa_FOUND)
    SET(OPENGL_FOUND TRUE)
  ENDIF()

  IF(OSMesa_GLU_FOUND)
    SET(OPENGL_GLU_FOUND TRUE)
  ENDIF()

  # Replace GL libraries with OSMesa ones, if available
  IF(OSMesa_FOUND)
    SET(OPENGL_INCLUDE_DIRS ${OSMesa_INCLUDE_DIRS} ${OSMesa_GLU_INCLUDE_DIRS})
    SET(OPENGL_CFLAGS ${OSMesa_CFLAGS} ${OSMesa_GLU_CFLAGS})
    SET(OPENGL_LIBRARIES ${OSMesa_GLU_LIBRARIES} ${OSMesa_LIBRARIES})
    SET(OPENGL_LIBRARY_DIRS ${OSMesa_LIBRARY_DIRS} ${OSMesa_GLU_LIBRARY_DIRS})
    SET(OPENGL_LDFLAGS ${OSMesa_LDFLAGS} ${OSMesa_GLU_LDFLAGS})
  ENDIF()

  ADD_DEFINITIONS(-DTHEA_GL_OSMESA=1)
ELSE()
  FIND_PACKAGE(OpenGL)
ENDIF()

# Dependency: wxWidgets
FIND_PACKAGE(wxWidgets COMPONENTS core base gl)
IF(wxWidgets_FOUND)
  INCLUDE(${wxWidgets_USE_FILE})
ENDIF()

# Dependency: OpenCL
FIND_PACKAGE(OpenCL)

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

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
                    ${SourceRoot}
                    ${OPENGL_INCLUDE_DIRS}
                    ${Boost_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS})

IF(OPENCL_FOUND)
  INCLUDE_DIRECTORIES(BEFORE ${OPENCL_INCLUDE_DIRS})
ENDIF()

# Link directories
LINK_DIRECTORIES(${OutputRoot}/lib
                 ${OPENGL_LIBRARY_DIRS})

# List of all tools
SET(TheaAllTools )

#===========================================================
# Browse3D
#===========================================================

IF(wxWidgets_FOUND AND OPENGL_FOUND)
  # Source file lists
  FILE(GLOB Browse3DSources
       ${SourceRoot}/Tools/Browse3D/*.cpp)

  # Libraries to link to
  SET(Browse3DLibraries
      Thea
      ${wxWidgets_LIBRARIES}
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
# GeoSphere
#===========================================================

# Source file lists
FILE(GLOB GeoSphereSources
     ${SourceRoot}/Tools/GeoSphere/*.cpp)

# Libraries to link to
SET(GeoSphereLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(GeoSphere ${GeoSphereSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(GeoSphere ${GeoSphereLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(GeoSphere PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(GeoSphere "${GeoSphereLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} GeoSphere)

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
# MeshGlobalFeatures
#===========================================================

# Source file lists
FILE(GLOB MeshGlobalFeaturesSources
     ${SourceRoot}/Tools/MeshGlobalFeatures/*.cpp)

# Libraries to link to
SET(MeshGlobalFeaturesLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshGlobalFeatures ${MeshGlobalFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshGlobalFeatures ${MeshGlobalFeaturesLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshGlobalFeatures PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshGlobalFeatures "${MeshGlobalFeaturesLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} MeshGlobalFeatures)

#===========================================================
# MeshLocalFeatures
#===========================================================

# Source file lists
FILE(GLOB MeshLocalFeaturesSources
     ${SourceRoot}/Tools/MeshLocalFeatures/*.cpp)

# Libraries to link to
SET(MeshLocalFeaturesLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshLocalFeatures ${MeshLocalFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshLocalFeatures ${MeshLocalFeaturesLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(MeshLocalFeatures PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(MeshLocalFeatures "${MeshLocalFeaturesLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} MeshLocalFeatures)

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
# RenderShape
#===========================================================

# Source file lists
FILE(GLOB RenderShapeSources
     ${SourceRoot}/Tools/RenderShape/*.cpp)

# Libraries to link to
SET(RenderShapeLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(RenderShape ${RenderShapeSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(RenderShape ${RenderShapeLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(RenderShape PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(RenderShape "${RenderShapeLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} RenderShape)

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

IF(CLUTO_FOUND)

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

ENDIF(CLUTO_FOUND)

#===========================================================
# ShapeAlign
#===========================================================

# Source file lists
FILE(GLOB ShapeAlignSources
     ${SourceRoot}/Tools/ShapeAlign/*.cpp)

# Libraries to link to
SET(ShapeAlignLibraries
    Thea
    ${Lib3ds_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(ShapeAlign ${ShapeAlignSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(ShapeAlign ${ShapeAlignLibraries})

# Linker flags
IF(Thea_LDFLAGS)
  SET_TARGET_PROPERTIES(ShapeAlign PROPERTIES LINK_FLAGS ${Thea_LDFLAGS})
ENDIF()

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(ShapeAlign "${ShapeAlignLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} ShapeAlign)

#===========================================================
# Target for all tools
#===========================================================

ADD_CUSTOM_TARGET(tools DEPENDS ${TheaAllTools})
