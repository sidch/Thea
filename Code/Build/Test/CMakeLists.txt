#===============================================================================================================================
#
# Build script for the Thea tests.
#
# Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaTests)

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# See cmake --help-policy CMP0003 for details on this one
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

# See cmake --help-policy CMP0074 for details on this one
IF(POLICY CMP0074)
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF(POLICY CMP0074)

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

# Path to folder with installations of the dependencies
IF(NOT THEA_INSTALLATIONS_ROOT)
  SET(THEA_INSTALLATIONS_ROOT ${CMAKE_INSTALL_PREFIX})
ENDIF()
SET(THEA_INSTALLATIONS_ROOT ${THEA_INSTALLATIONS_ROOT} CACHE PATH "Path to folder with installations of dependencies")

# Locate dependencies
SET(Thea_FIND_ALL TRUE)
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Dependency: OpenGL
FIND_PACKAGE(OpenGL REQUIRED)
IF(NOT OPENGL_GLU_FOUND)
  MESSAGE(FATAL_ERROR "GLU was not found")
ENDIF()

# Dependency: GLUT
FIND_PACKAGE(GLUT REQUIRED)
IF(NOT GLUT_FOUND)
  MESSAGE(FATAL_ERROR "GLUT not found")
ENDIF()
MESSAGE(STATUS "Found GLUT: headers at ${GLUT_INCLUDE_DIR}, libraries at ${GLUT_LIBRARIES}")

# Definitions, compiler switches etc.
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")

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
                    ${SourceRoot}
                    ${Boost_INCLUDE_DIRS}
                    ${GLUT_INCLUDE_DIR}
                    ${OPENGL_INCLUDE_DIRS}
                    ${CGAL_INCLUDE_DIRS}
                    ${Lib3ds_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${OutputRoot}/lib)

#===========================================================
# TestARPACK
#===========================================================

IF(TARGET TheaPluginARPACK)

# Source file lists
SET(TheaTestARPACKSources
      ${SourceRoot}/Test/TestARPACK.cpp)

# Libraries to link to
SET(TheaTestARPACKLibraries
      Thea
      TheaPluginARPACK
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestARPACK ${TheaTestARPACKSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestARPACK ${TheaTestARPACKLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestARPACK "${TheaTestARPACKLibraries}")
ENDIF()

ENDIF(TARGET TheaPluginARPACK)

#===========================================================
# TestBagOfWords
#===========================================================

# Source file lists
SET(TheaTestBagOfWordsSources
      ${SourceRoot}/Test/TestBagOfWords.cpp)

# Libraries to link to
SET(TheaTestBagOfWordsLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestBagOfWords ${TheaTestBagOfWordsSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestBagOfWords ${TheaTestBagOfWordsLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestBagOfWords "${TheaTestBagOfWordsLibraries}")
ENDIF()

#===========================================================
# TestCSPARSE
#===========================================================

# Source file lists
SET(TheaTestCSPARSESources
      ${SourceRoot}/Test/TestCSPARSE.cpp)

# Libraries to link to
SET(TheaTestCSPARSELibraries
      Thea
      TheaPluginCSPARSE
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestCSPARSE ${TheaTestCSPARSESources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestCSPARSE ${TheaTestCSPARSELibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestCSPARSE "${TheaTestCSPARSELibraries}")
ENDIF()

#===========================================================
# TestDisplayMesh
#===========================================================

# Source file lists
SET(TheaTestDisplayMeshSources
      ${SourceRoot}/Test/TestDisplayMesh.cpp)

# Libraries to link to
SET(TheaTestDisplayMeshLibraries
      Thea
      TheaPluginGL
      ${GLUT_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestDisplayMesh ${TheaTestDisplayMeshSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestDisplayMesh ${TheaTestDisplayMeshLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestDisplayMesh "${TheaTestDisplayMeshLibraries}")
ENDIF()

#===========================================================
# TestGL
#===========================================================

# Source file lists
SET(TheaTestGLSources
      ${SourceRoot}/Test/TestGL.cpp)

# Libraries to link to
SET(TheaTestGLLibraries
      Thea
      TheaPluginGL
      ${GLUT_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestGL ${TheaTestGLSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestGL ${TheaTestGLLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestGL "${TheaTestGLLibraries}")
ENDIF()

#===========================================================
# TestHoughForest
#===========================================================

# Source file lists
SET(TheaTestHoughForestSources
      ${SourceRoot}/Test/TestHoughForest.cpp)

# Libraries to link to
SET(TheaTestHoughForestLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestHoughForest ${TheaTestHoughForestSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestHoughForest ${TheaTestHoughForestLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestHoughForest "${TheaTestHoughForestLibraries}")
ENDIF()

#===========================================================
# TestJointBoost
#===========================================================

# Source file lists
SET(TheaTestJointBoostSources
      ${SourceRoot}/Test/TestJointBoost.cpp)

# Libraries to link to
SET(TheaTestJointBoostLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestJointBoost ${TheaTestJointBoostSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestJointBoost ${TheaTestJointBoostLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestJointBoost "${TheaTestJointBoostLibraries}")
ENDIF()

#===========================================================
# TestKDTree3
#===========================================================

# Source file lists
SET(TheaTestKDTree3Sources
      ${SourceRoot}/Test/TestKDTree3.cpp)

# Libraries to link to
SET(TheaTestKDTree3Libraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestKDTree3 ${TheaTestKDTree3Sources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestKDTree3 ${TheaTestKDTree3Libraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestKDTree3 "${TheaTestKDTree3Libraries}")
ENDIF()

#===========================================================
# TestLaplaceBeltrami
#===========================================================

# Source file lists
SET(TheaTestLaplaceBeltramiSources
      ${SourceRoot}/Test/TestLaplaceBeltrami.cpp)

# Libraries to link to
SET(TheaTestLaplaceBeltramiLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestLaplaceBeltrami ${TheaTestLaplaceBeltramiSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestLaplaceBeltrami ${TheaTestLaplaceBeltramiLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestLaplaceBeltrami "${TheaTestLaplaceBeltramiLibraries}")
ENDIF()

#===========================================================
# TestMath
#===========================================================

# Source file lists
SET(TheaTestMathSources
      ${SourceRoot}/Test/TestMath.cpp)

# Libraries to link to
SET(TheaTestMathLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMath ${TheaTestMathSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMath ${TheaTestMathLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestMath "${TheaTestMathLibraries}")
ENDIF()

#===========================================================
# TestMetrics
#===========================================================

# Source file lists
SET(TheaTestMetricsSources
      ${SourceRoot}/Test/TestMetrics.cpp)

# Libraries to link to
SET(TheaTestMetricsLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMetrics ${TheaTestMetricsSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMetrics ${TheaTestMetricsLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestMetrics "${TheaTestMetricsLibraries}")
ENDIF()

#===========================================================
# TestMesh
#===========================================================

# Source file lists
SET(TheaTestMeshSources
      ${SourceRoot}/Test/TestMesh.cpp)

# Libraries to link to
SET(TheaTestMeshLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMesh ${TheaTestMeshSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMesh ${TheaTestMeshLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestMesh "${TheaTestMeshLibraries}")
ENDIF()

#===========================================================
# TestOPTPP
#===========================================================

# Source file lists
SET(TheaTestOPTPPSources
      ${SourceRoot}/Test/TestOPTPP.cpp)

# Libraries to link to
SET(TheaTestOPTPPLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestOPTPP ${TheaTestOPTPPSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestOPTPP ${TheaTestOPTPPLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestOPTPP "${TheaTestOPTPPLibraries}")
ENDIF()

#===========================================================
# TestPCA
#===========================================================

# Source file lists
SET(TheaTestPCASources
      ${SourceRoot}/Test/TestPCA.cpp)

# Libraries to link to
SET(TheaTestPCALibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestPCA ${TheaTestPCASources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestPCA ${TheaTestPCALibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestPCA "${TheaTestPCALibraries}")
ENDIF()

#===========================================================
# TestPyramidMatch
#===========================================================

# Source file lists
SET(TheaTestPyramidMatchSources
      ${SourceRoot}/Test/TestPyramidMatch.cpp)

# Libraries to link to
SET(TheaTestPyramidMatchLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestPyramidMatch ${TheaTestPyramidMatchSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestPyramidMatch ${TheaTestPyramidMatchLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestPyramidMatch "${TheaTestPyramidMatchLibraries}")
ENDIF()

#===========================================================
# TestZernike
#===========================================================

# Source file lists
SET(TheaTestZernikeSources
      ${SourceRoot}/Test/TestZernike.cpp)

# Libraries to link to
SET(TheaTestZernikeLibraries
      Thea
      ${PLATFORM_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestZernike ${TheaTestZernikeSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestZernike ${TheaTestZernikeLibraries})

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestZernike "${TheaTestZernikeLibraries}")
ENDIF()

#===========================================================
# Target for all tests
#===========================================================

SET(TheaTestsDependencies
    TheaTestBagOfWords
    TheaTestCSPARSE
    TheaTestDisplayMesh
    TheaTestGL
    TheaTestJointBoost
    TheaTestKDTree3
    TheaTestMath
    TheaTestMesh
    TheaTestMetrics
    TheaTestOPTPP
    TheaTestPCA
    TheaTestPyramidMatch
    TheaTestZernike)

IF(TARGET TheaTestARPACK)
  SET(TheaTestsDependencies ${TheaTestsDependencies} TheaTestARPACK)
ENDIF()

FOREACH(target ${TheaTestsDependencies})
  ADD_TEST(NAME ${target} COMMAND ${target})
ENDFOREACH()
