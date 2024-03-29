#===============================================================================================================================
#
# Build script for the Thea tools.
#
# Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

PROJECT(TheaTools CXX C)

# Shared library flags
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Dependency: OpenGL
IF(WITH_PLUGIN_GL)
  IF(WITH_OSMESA)
    SET(OPENGL_FOUND )
    SET(OPENGL_GLU_FOUND )

    IF(EXISTS ${THEA_DEPS_ROOT}/installed-osmesa)
      SET(OSMesa_ROOT ${THEA_DEPS_ROOT}/installed-osmesa)
    ELSE()
      SET(OSMesa_ROOT ${THEA_DEPS_ROOT})
    ENDIF()
    FIND_PACKAGE(OSMesa)

    IF(OSMesa_GLU_FOUND)
      SET(OPENGL_GLU_FOUND TRUE)
    ENDIF()

    IF(OSMesa_FOUND)
      SET(OPENGL_FOUND TRUE)

      # Replace GL libraries with OSMesa ones, if available
      SET(OPENGL_INCLUDE_DIRS ${OSMesa_INCLUDE_DIRS} ${OSMesa_GLU_INCLUDE_DIRS})
      SET(OPENGL_CFLAGS "${OSMesa_CFLAGS} ${OSMesa_GLU_CFLAGS}")
      SET(OPENGL_LIBRARIES ${OSMesa_GLU_LIBRARIES} ${OSMesa_LIBRARIES})
      SET(OPENGL_LIBRARY_DIRS ${OSMesa_LIBRARY_DIRS} ${OSMesa_GLU_LIBRARY_DIRS})
      SET(OPENGL_LDFLAGS "${OSMesa_LDFLAGS} ${OSMesa_GLU_LDFLAGS}")
      ADD_DEFINITIONS(-DTHEA_GL_OSMESA=1 -DGLEW_OSMESA=1)
    ENDIF()
  ELSE()
    FIND_PACKAGE(OpenGL)
  ENDIF()

  IF(OPENGL_FOUND)
    ADD_DEFINITIONS(-DGLEW_STATIC)

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENGL_CFLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENGL_CFLAGS}")
  ENDIF()

  # Dependency: wxWidgets
  FIND_PACKAGE(wxWidgets COMPONENTS core base gl)
  IF(wxWidgets_FOUND)
    INCLUDE(${wxWidgets_USE_FILE})
  ENDIF()
ENDIF()

# Dependency: OpenCL
FIND_PACKAGE(OpenCL)
IF(OPENCL_FOUND AND APPLE)
  ADD_DEFINITIONS(-DGL_SILENCE_DEPRECATION -DCL_SILENCE_DEPRECATION)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE
                    ${SourceRoot}
                    ${OPENGL_INCLUDE_DIRS}
                    ${OPENCL_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${OutputRoot}/lib
                 ${OPENGL_LIBRARY_DIRS})

# List of all tools
SET(TheaAllTools )

#===========================================================
# Browse3D
#===========================================================

IF(WITH_PLUGIN_GL AND wxWidgets_FOUND AND OPENGL_FOUND AND NOT WITH_OSMESA)
  # Source file lists
  FILE(GLOB Browse3DSources
       ${SourceRoot}/Tools/Browse3D/*.cpp
       ${SourceRoot}/Plugins/GL/glew.c)

  # Libraries to link to
  SET(Browse3DLibraries
      Thea
      ${wxWidgets_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${Thea_DEPS_LIBRARIES})

  # Additional platform-specific libraries
  IF(WIN32)
    SET(Browse3DLibraries ${Browse3DLibraries} "version.lib")
  ENDIF()

  # Build products
  ADD_EXECUTABLE(Browse3D ${Browse3DSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(Browse3D ${Browse3DLibraries})

  # Linker flags
  SET_TARGET_PROPERTIES(Browse3D
                          PROPERTIES
                            LINK_FLAGS "${Thea_DEPS_LDFLAGS} ${OPENGL_LDFLAGS}"
                            INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(GeoSphere ${GeoSphereSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(GeoSphere ${GeoSphereLibraries})

# Linker flags
SET_TARGET_PROPERTIES(GeoSphere
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
      ${Thea_DEPS_LIBRARIES})

  # Build products
  ADD_EXECUTABLE(gip ${GIPSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(gip ${GIPLibraries})

  # Linker flags
  SET_TARGET_PROPERTIES(gip
                          PROPERTIES
                            LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                            INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(ImageFeatures ${ImageFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(ImageFeatures ${ImageFeaturesLibraries})

# Linker flags
SET_TARGET_PROPERTIES(ImageFeatures
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshConv ${MeshConvSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshConv ${MeshConvLibraries})

# Linker flags
SET_TARGET_PROPERTIES(MeshConv
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshGlobalFeatures ${MeshGlobalFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshGlobalFeatures ${MeshGlobalFeaturesLibraries})

# Linker flags
SET_TARGET_PROPERTIES(MeshGlobalFeatures
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshLocalFeatures ${MeshLocalFeaturesSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshLocalFeatures ${MeshLocalFeaturesLibraries})

# Linker flags
SET_TARGET_PROPERTIES(MeshLocalFeatures
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshFix ${MeshFixSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshFix ${MeshFixLibraries})

# Linker flags
SET_TARGET_PROPERTIES(MeshFix
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(MeshSample ${MeshSampleSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(MeshSample ${MeshSampleLibraries})

# Linker flags
SET_TARGET_PROPERTIES(MeshSample
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(Register ${RegisterSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(Register ${RegisterLibraries})

# Linker flags
SET_TARGET_PROPERTIES(Register
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(Register "${RegisterLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} Register)

#===========================================================
# RenderShape
#===========================================================

IF(WITH_PLUGIN_GL AND OPENGL_FOUND)
  # Source file lists
  FILE(GLOB RenderShapeSources
       ${SourceRoot}/Tools/RenderShape/*.cpp
       ${SourceRoot}/Plugins/GL/glew.c)

  # Libraries to link to
  SET(RenderShapeLibraries
      Thea
      ${OPENGL_LIBRARIES}
      ${Thea_DEPS_LIBRARIES})

  # Additional platform-specific libraries
  IF(WIN32)
    SET(RenderShapeLibraries ${RenderShapeLibraries} "version.lib")
  ENDIF()

  # Build products
  ADD_EXECUTABLE(RenderShape ${RenderShapeSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(RenderShape ${RenderShapeLibraries})

  # Linker flags
  SET_TARGET_PROPERTIES(RenderShape
                          PROPERTIES
                            LINK_FLAGS "${Thea_DEPS_LDFLAGS} ${OPENGL_LDFLAGS}"
                            INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Fix library install names on OS X
  IF(APPLE)
    INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
    OSX_FIX_DYLIB_REFERENCES(RenderShape "${RenderShapeLibraries}")
  ENDIF()

  SET(TheaAllTools ${TheaAllTools} RenderShape)
ENDIF()

#===========================================================
# SampleGraph
#===========================================================

# Source file lists
FILE(GLOB SampleGraphSources
     ${SourceRoot}/Tools/SampleGraph/*.cpp)

# Libraries to link to
SET(SampleGraphLibraries
    Thea
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(SampleGraph ${SampleGraphSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(SampleGraph ${SampleGraphLibraries})

# Linker flags
SET_TARGET_PROPERTIES(SampleGraph
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(SampleGraph "${SampleGraphLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} SampleGraph)

#===========================================================
# SegmentSDF
#===========================================================

IF(WITH_CLUTO AND CLUTO_FOUND)
  # Source file lists
  FILE(GLOB SegmentSDFSources
       ${SourceRoot}/Tools/SegmentSDF/*.cpp)

  # Libraries to link to
  SET(SegmentSDFLibraries
      Thea
      ${Thea_DEPS_LIBRARIES})

  # Build products
  ADD_EXECUTABLE(SegmentSDF ${SegmentSDFSources})

  # Additional libraries to be linked
  TARGET_LINK_LIBRARIES(SegmentSDF ${SegmentSDFLibraries})

  # Linker flags
  SET_TARGET_PROPERTIES(SegmentSDF
                          PROPERTIES
                            LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                            INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Fix library install names on OS X
  IF(APPLE)
    INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
    OSX_FIX_DYLIB_REFERENCES(SegmentSDF "${SegmentSDFLibraries}")
  ENDIF()

  SET(TheaAllTools ${TheaAllTools} SegmentSDF)
ENDIF()

#===========================================================
# ShapeAlign
#===========================================================

# Source file lists
FILE(GLOB ShapeAlignSources
     ${SourceRoot}/Tools/ShapeAlign/*.cpp)

# Libraries to link to
SET(ShapeAlignLibraries
    Thea
    ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(ShapeAlign ${ShapeAlignSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(ShapeAlign ${ShapeAlignLibraries})

# Linker flags
SET_TARGET_PROPERTIES(ShapeAlign
                        PROPERTIES
                          LINK_FLAGS "${Thea_DEPS_LDFLAGS}"
                          INSTALL_RPATH_USE_LINK_PATH TRUE)

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(ShapeAlign "${ShapeAlignLibraries}")
ENDIF()

SET(TheaAllTools ${TheaAllTools} ShapeAlign)

#===========================================================
# Target for all tools, and installation
#===========================================================

ADD_CUSTOM_TARGET(tools DEPENDS ${TheaAllTools})

FOREACH(tgt ${TheaAllTools})
  INSTALL(TARGETS ${tgt} DESTINATION share/Thea/Build/Output/bin)
ENDFOREACH()

install_symlink(${CMAKE_INSTALL_PREFIX}/share/Thea/Build/Output/bin
                ${CMAKE_INSTALL_PREFIX}/bin/Thea)

INSTALL(DIRECTORY ${Thea_PROJECT_ROOT}/Resources DESTINATION share/Thea)
