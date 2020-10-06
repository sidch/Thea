#===============================================================================================================================
#
# Build script for the Thea tests.
#
# Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
#
#===============================================================================================================================

PROJECT(TheaTests CXX C)

# Setup build
GET_FILENAME_COMPONENT(Thea_PROJECT_ROOT ../.. ABSOLUTE)
SET(Thea_FIND_ALL TRUE)
INCLUDE(${Thea_PROJECT_ROOT}/Build/Common/BuildSetup.cmake NO_POLICY_SCOPE)

# Dependency: OpenGL
IF(WITH_OSMESA)
  SET(OPENGL_FOUND FALSE)
ELSE()
  FIND_PACKAGE(OpenGL)
  IF(NOT OPENGL_GLU_FOUND)
    SET(OPENGL_FOUND FALSE)
  ENDIF()
ENDIF()

# Dependency: GLUT
IF(OPENGL_FOUND)
  FIND_PACKAGE(GLUT)
  MESSAGE(STATUS "Found GLUT: headers at ${GLUT_INCLUDE_DIR}, libraries at ${GLUT_LIBRARIES}")
ENDIF()

IF(OPENGL_FOUND)
  INCLUDE_DIRECTORIES(BEFORE
                      ${GLUT_INCLUDE_DIR}
                      ${OPENGL_INCLUDE_DIRS})

  IF(APPLE)
    ADD_DEFINITIONS(-DGL_SILENCE_DEPRECATION -DCL_SILENCE_DEPRECATION)
  ENDIF()
ENDIF()

# Shared library flags
IF(THEA_DLL)
  ADD_DEFINITIONS(-DTHEA_DLL -DTHEA_DLL_IMPORTS)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE ${SourceRoot})

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestARPACK ${TheaTestARPACKSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestARPACK ${TheaTestARPACKLibraries})
SET_TARGET_PROPERTIES(TheaTestARPACK PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestBagOfWords ${TheaTestBagOfWordsSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestBagOfWords ${TheaTestBagOfWordsLibraries})
SET_TARGET_PROPERTIES(TheaTestBagOfWords PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestBagOfWords "${TheaTestBagOfWordsLibraries}")
ENDIF()

#===========================================================
# TestBezier
#===========================================================

# Source file lists
SET(TheaTestBezierSources
      ${SourceRoot}/Test/TestBezier.cpp)

# Libraries to link to
SET(TheaTestBezierLibraries
      Thea
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestBezier ${TheaTestBezierSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestBezier ${TheaTestBezierLibraries})
SET_TARGET_PROPERTIES(TheaTestBezier PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestBezier "${TheaTestBezierLibraries}")
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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestCSPARSE ${TheaTestCSPARSESources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestCSPARSE ${TheaTestCSPARSELibraries})
SET_TARGET_PROPERTIES(TheaTestCSPARSE PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestCSPARSE "${TheaTestCSPARSELibraries}")
ENDIF()

#===========================================================
# TestDisplayMesh
#===========================================================

IF(OPENGL_FOUND AND GLUT_FOUND)

# Source file lists
SET(TheaTestDisplayMeshSources
      ${SourceRoot}/Test/TestDisplayMesh.cpp)

# Libraries to link to
SET(TheaTestDisplayMeshLibraries
      Thea
      TheaPluginGL
      ${GLUT_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestDisplayMesh ${TheaTestDisplayMeshSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestDisplayMesh ${TheaTestDisplayMeshLibraries})
SET_TARGET_PROPERTIES(TheaTestDisplayMesh PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestDisplayMesh "${TheaTestDisplayMeshLibraries}")
ENDIF()

ENDIF()

#===========================================================
# TestGL
#===========================================================

IF(OPENGL_FOUND AND GLUT_FOUND)

# Source file lists
SET(TheaTestGLSources
      ${SourceRoot}/Test/TestGL.cpp)

# Libraries to link to
SET(TheaTestGLLibraries
      Thea
      TheaPluginGL
      ${GLUT_LIBRARIES}
      ${OPENGL_LIBRARIES}
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestGL ${TheaTestGLSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestGL ${TheaTestGLLibraries})
SET_TARGET_PROPERTIES(TheaTestGL PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestGL "${TheaTestGLLibraries}")
ENDIF()

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestHoughForest ${TheaTestHoughForestSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestHoughForest ${TheaTestHoughForestLibraries})
SET_TARGET_PROPERTIES(TheaTestHoughForest PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestHoughForest "${TheaTestHoughForestLibraries}")
ENDIF()

#===========================================================
# TestIterators
#===========================================================

# Source file lists
SET(TheaTestIteratorsSources
      ${SourceRoot}/Test/TestIterators.cpp)

# Libraries to link to
SET(TheaTestIteratorsLibraries
      Thea
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestIterators ${TheaTestIteratorsSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestIterators ${TheaTestIteratorsLibraries})
SET_TARGET_PROPERTIES(TheaTestIterators PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestIterators "${TheaTestIteratorsLibraries}")
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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestJointBoost ${TheaTestJointBoostSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestJointBoost ${TheaTestJointBoostLibraries})
SET_TARGET_PROPERTIES(TheaTestJointBoost PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestKDTree3 ${TheaTestKDTree3Sources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestKDTree3 ${TheaTestKDTree3Libraries})
SET_TARGET_PROPERTIES(TheaTestKDTree3 PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestLaplaceBeltrami ${TheaTestLaplaceBeltramiSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestLaplaceBeltrami ${TheaTestLaplaceBeltramiLibraries})
SET_TARGET_PROPERTIES(TheaTestLaplaceBeltrami PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMath ${TheaTestMathSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMath ${TheaTestMathLibraries})
SET_TARGET_PROPERTIES(TheaTestMath PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMetrics ${TheaTestMetricsSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMetrics ${TheaTestMetricsLibraries})
SET_TARGET_PROPERTIES(TheaTestMetrics PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestMesh ${TheaTestMeshSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestMesh ${TheaTestMeshLibraries})
SET_TARGET_PROPERTIES(TheaTestMesh PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestMesh "${TheaTestMeshLibraries}")
ENDIF()

#===========================================================
# TestOPTPP
#===========================================================

IF(WITH_OPTPP)

# Source file lists
SET(TheaTestOPTPPSources
      ${SourceRoot}/Test/TestOPTPP.cpp)

# Libraries to link to
SET(TheaTestOPTPPLibraries
      Thea
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestOPTPP ${TheaTestOPTPPSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestOPTPP ${TheaTestOPTPPLibraries})
SET_TARGET_PROPERTIES(TheaTestOPTPP PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestOPTPP "${TheaTestOPTPPLibraries}")
ENDIF()

ENDIF()

#===========================================================
# TestParam
#===========================================================

# Source file lists
SET(TheaTestParamSources
      ${SourceRoot}/Test/TestParam.cpp)

# Libraries to link to
SET(TheaTestParamLibraries
      Thea
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestParam ${TheaTestParamSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestParam ${TheaTestParamLibraries})
SET_TARGET_PROPERTIES(TheaTestParam PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

# Fix library install names on OS X
IF(APPLE)
  INCLUDE(${CMAKE_MODULE_PATH}/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(TheaTestParam "${TheaTestParamLibraries}")
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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestPCA ${TheaTestPCASources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestPCA ${TheaTestPCALibraries})
SET_TARGET_PROPERTIES(TheaTestPCA PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestPyramidMatch ${TheaTestPyramidMatchSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestPyramidMatch ${TheaTestPyramidMatchLibraries})
SET_TARGET_PROPERTIES(TheaTestPyramidMatch PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
      ${Thea_DEPS_LIBRARIES})

# Build products
ADD_EXECUTABLE(TheaTestZernike ${TheaTestZernikeSources})

# Additional libraries to be linked
TARGET_LINK_LIBRARIES(TheaTestZernike ${TheaTestZernikeLibraries})
SET_TARGET_PROPERTIES(TheaTestZernike PROPERTIES LINK_FLAGS "${Thea_DEPS_LDFLAGS}")

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
